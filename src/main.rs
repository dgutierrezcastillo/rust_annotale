use anyhow::{Context, Result};
use bio::io::fasta;
use bio::alphabets::dna::revcomp;
use bio::seq_analysis::orf::{Finder, Orf};
use hmmer_pure_rs::alphabet::{Alphabet, AlphabetType};
use hmmer_pure_rs::bg::Bg;
use hmmer_pure_rs::hmmfile;
use hmmer_pure_rs::profile::{profile_config, P7_LOCAL};
use hmmer_pure_rs::sequence::Sequence as HmmSequence;
use hmmer_pure_rs::{Hmm, Pipeline, Profile, TopHits, OProfile};
use std::path::Path;
use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;


#[derive(Parser, Debug)]
#[command(author, version, about = "A Rust implementation of AnnoTALE")]
struct Args {
    #[arg(short, long)]
    fasta: String,

    #[arg(long)]
    hmm_dir: String,

    #[arg(short, long, default_value_t = 100.0)]
    threshold: f32,
}

// --- Translation Logic ---

struct Translator {
    table: HashMap<&'static [u8; 3], u8>,
}

impl Translator {
    fn new() -> Self {
        let mut table = HashMap::new();
        let codes = [
            (b"TTT", b'F'), (b"TTC", b'F'), (b"TTA", b'L'), (b"TTG", b'L'),
            (b"TCT", b'S'), (b"TCC", b'S'), (b"TCA", b'S'), (b"TCG", b'S'),
            (b"TAT", b'Y'), (b"TAC", b'Y'), (b"TAA", b'*'), (b"TAG", b'*'),
            (b"TGT", b'C'), (b"TGC", b'C'), (b"TGA", b'*'), (b"TGG", b'W'),
            (b"CTT", b'L'), (b"CTC", b'L'), (b"CTA", b'L'), (b"CTG", b'L'),
            (b"CCT", b'P'), (b"CCC", b'P'), (b"CCA", b'P'), (b"CCG", b'P'),
            (b"CAT", b'H'), (b"CAC", b'H'), (b"CAA", b'Q'), (b"CAG", b'Q'),
            (b"CGT", b'R'), (b"CGC", b'R'), (b"CGA", b'R'), (b"CGG", b'R'),
            (b"ATT", b'I'), (b"ATC", b'I'), (b"ATA", b'I'), (b"ATG", b'M'),
            (b"ACT", b'T'), (b"ACC", b'T'), (b"ACA", b'T'), (b"ACG", b'T'),
            (b"AAT", b'N'), (b"AAC", b'N'), (b"AAA", b'K'), (b"AAG", b'K'),
            (b"AGT", b'S'), (b"AGC", b'S'), (b"AGA", b'R'), (b"AGG", b'R'),
            (b"GTT", b'V'), (b"GTC", b'V'), (b"GTA", b'V'), (b"GTG", b'V'),
            (b"GCT", b'A'), (b"GCC", b'A'), (b"GCA", b'A'), (b"GCG", b'A'),
            (b"GAT", b'D'), (b"GAC", b'D'), (b"GAA", b'E'), (b"GAG", b'E'),
            (b"GGT", b'G'), (b"GGC", b'G'), (b"GGA", b'G'), (b"GGG", b'G'),
        ];
        for (codon, aa) in codes.iter() {
            table.insert(*codon, *aa);
        }
        Translator { table }
    }

    fn translate(&self, seq: &[u8]) -> Vec<u8> {
        let mut protein = Vec::new();
        for chunk in seq.chunks_exact(3) {
            let upper_chunk = [
                chunk[0].to_ascii_uppercase(),
                chunk[1].to_ascii_uppercase(),
                chunk[2].to_ascii_uppercase(),
            ];
            protein.push(*self.table.get(&upper_chunk).unwrap_or(&b'X'));
        }
        protein
    }
}

// --- TALE Finder Logic ---

#[derive(Debug, Clone)]
struct TALERegion {
    strand: char, 
    start: usize, // Genomic start
    end: usize,   // Genomic end
    score: f32,
    cds_start: usize, 
    cds_end: usize,   
    is_pseudo: bool,
    rvds: String,
}

struct TALEFinder {
    repeats_hmm: Hmm,
    abc: Alphabet,
    bg: Bg,
    threshold: f32,
    translator: Translator,
}

impl TALEFinder {
    fn new(hmm_dir: &str, threshold: f32) -> Result<Self> {
        let repeats_path = Path::new(hmm_dir).join("repeats.hmm");
        let repeats_hmms = hmmfile::read_hmm_file(&repeats_path)
            .with_context(|| format!("Failed to read {}", repeats_path.display()))?;

        let abc = Alphabet::new(AlphabetType::Dna);
        let bg = Bg::new(&abc);

        Ok(Self {
            repeats_hmm: repeats_hmms[0].clone(),
            abc,
            bg,
            threshold,
            translator: Translator::new(),
        })
    }

    fn scan_sequence(&self, record_id: &str, sequence: &[u8]) -> Vec<TALERegion> {
        let mut results = Vec::new();

        // 1. Scan Forward Strand
        results.extend(self.process_strand(record_id, sequence, '+'));

        // 2. Scan Reverse Strand
        let rev_seq = revcomp(sequence);
        results.extend(self.process_strand(record_id, &rev_seq, '-'));

        results
    }

    fn process_strand(&self, id: &str, sequence: &[u8], strand: char) -> Vec<TALERegion> {
        let mut raw_matches = self.run_hmmer_raw(id, sequence);
        if raw_matches.is_empty() { return Vec::new(); }

        raw_matches.sort_by_key(|m| m.0);
        let mut clusters = Vec::new();
        if !raw_matches.is_empty() {
            let mut current_cluster = Vec::new();
            current_cluster.push(raw_matches[0]);
            let mut c_end = raw_matches[0].1;
            let mut c_score = raw_matches[0].2;

            for i in 1..raw_matches.len() {
                let m = raw_matches[i];
                if m.0 < c_end + 500 {
                    c_end = std::cmp::max(c_end, m.1);
                    c_score += m.2;
                    current_cluster.push(m);
                } else {
                    clusters.push((current_cluster.clone(), c_score));
                    current_cluster.clear();
                    current_cluster.push(m);
                    c_end = m.1;
                    c_score = m.2;
                }
            }
            clusters.push((current_cluster, c_score));
        }

        let mut final_tales = Vec::new();
        for (domains, c_score) in clusters {
            let c_start = domains[0].0;
            let c_end = domains.last().unwrap().1;
            
            // Refine boundaries
            let buffer = 1200;
            let search_start = c_start.saturating_sub(buffer);
            let search_end = std::cmp::min(sequence.len(), c_end + buffer);
            
            // Refine CDS within search window
            let (cds_rel_start, cds_rel_end, is_pseudo) = self.refine_cds(sequence, search_start, search_end);
            let final_cds_start = search_start + cds_rel_start;
            let final_cds_end = search_start + cds_rel_end;

            // Extract RVDs
            let rvds = if !is_pseudo && (final_cds_end - final_cds_start) > 100 {
                self.extract_rvds(sequence, final_cds_start, &domains)
            } else {
                "N/A".to_string()
            };

            let (actual_start, actual_end) = if strand == '+' {
                (final_cds_start, final_cds_end)
            } else {
                let l = sequence.len();
                (l - final_cds_end, l - final_cds_start)
            };

            final_tales.push(TALERegion {
                strand,
                start: actual_start,
                end: actual_end,
                score: c_score,
                cds_start: actual_start, 
                cds_end: actual_end,
                is_pseudo,
                rvds,
            });
        }

        final_tales
    }

    fn run_hmmer_raw(&self, id: &str, sequence: &[u8]) -> Vec<(usize, usize, f32)> {
        let mut matches = Vec::new();
        let l_val = sequence.len();
        let mut gm = Profile::new(self.repeats_hmm.m, &self.abc);
        profile_config(&self.repeats_hmm, &self.bg, &mut gm, l_val as i32, P7_LOCAL);
        let mut om = OProfile::convert(&gm);

        let mut pli = Pipeline::new();
        pli.new_model(&gm);
        let mut th = TopHits::new();

        let dsq = self.abc.digitize(sequence);
        let sq = HmmSequence {
            name: id.to_string(),
            acc: String::new(),
            desc: String::new(),
            dsq,
            n: l_val,
            l: l_val,
        };

        pli.run(&mut gm, &mut om, &self.bg, &self.repeats_hmm, &sq, &mut th);

        for hit in &th.hits {
            for domain in &hit.dcl {
                if domain.bitscore > 10.0 {
                    matches.push((domain.iali as usize, domain.jali as usize, domain.bitscore));
                }
            }
        }
        matches
    }

    fn refine_cds(&self, sequence: &[u8], start: usize, end: usize) -> (usize, usize, bool) {
        let target_seq = &sequence[start..end];
        let start_codons = vec![b"ATG"];
        let stop_codons = vec![b"TGA", b"TAG", b"TAA"];
        let finder = Finder::new(start_codons, stop_codons, 300); // Min ORF 300bp

        let mut max_len = 0;
        let mut best_orf: Option<Orf> = None;

        for orf in finder.find_all(target_seq) {
            let len = orf.end - orf.start;
            if len > max_len {
                max_len = len;
                best_orf = Some(orf);
            }
        }

        if let Some(orf) = best_orf {
            // TALE CDS are typically > 1.5kb. If found ORF is < 1kb but we found repeats, it might be pseudo.
            let is_pseudo = (orf.end - orf.start) < 900;
            (orf.start, orf.end, is_pseudo)
        } else {
            (0, 0, true)
        }
    }

    fn extract_rvds(&self, sequence: &[u8], cds_start: usize, domains: &[(usize, usize, f32)]) -> String {
        let mut rvd_str = String::new();
        
        for m in domains {
            let domain_start = m.0;
            if domain_start < cds_start { continue; }
            
            // Find the offset from the CDS start
            let offset = domain_start - cds_start;
            
            // Align to the nearest codon
            let frame_shift = offset % 3;
            let aligned_start = domain_start - frame_shift;
            
            // The RVD is at amino acids 12 and 13 of the repeat.
            // That is 12 * 3 = 36 bp from the aligned start of the repeat.
            let rvd_dna_start = aligned_start + 36;
            let rvd_dna_end = rvd_dna_start + 6;
            
            if rvd_dna_end <= sequence.len() {
                let rvd_dna = &sequence[rvd_dna_start..rvd_dna_end];
                let rvd_aa = self.translator.translate(rvd_dna);
                if rvd_aa.len() >= 2 {
                    if !rvd_str.is_empty() { rvd_str.push('-'); }
                    rvd_str.push(rvd_aa[0] as char);
                    rvd_str.push(rvd_aa[1] as char);
                }
            }
        }
        
        if rvd_str.is_empty() {
            "N/A".to_string()
        } else {
            rvd_str
        }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    println!("Initializing TALEFinder with HMMs from {}...", args.hmm_dir);
    let finder = TALEFinder::new(&args.hmm_dir, args.threshold)?;

    let reader = fasta::Reader::from_file(&args.fasta)?;

    println!("Scanning {} for TALE effectors...", args.fasta);

    // Read all records to parallelize
    let mut records = Vec::new();
    for result in reader.records() {
        records.push(result?);
    }

    records.par_iter().for_each(|record| {
        let id = record.id();
        let seq = record.seq();
        
        println!("Processing sequence: {} (length: {})", id, seq.len());
        let mut matches = finder.scan_sequence(id, seq);
        
        if !matches.is_empty() {
            println!("\nFound {} potential TALE effectors in {}", matches.len(), id);
            matches.sort_by_key(|m| m.start);
            
            println!("{:<5} | {:<8} | {:<2} | {:<10} | {:<10} | {:<8} | {:<10}", 
                     "No.", "Type", "St", "Start", "End", "Score", "RVDs");
            println!("{:-<5}-|-{:-<8}-|-{:-<2}-|-{:-<10}-|-{:-<10}-|-{:-<8}-|-{:-<10}", 
                     "", "", "", "", "", "", "");

            for (i, region) in matches.iter().enumerate() {
                let status = if region.is_pseudo { "PSEUDO" } else { "CDS" };
                println!("{:5} | {:<8} | {:<2} | {:10} | {:10} | {:8.1} | {}", 
                    i + 1, status, region.strand, region.start, region.end, region.score, region.rvds);
            }
        } else {
            println!("No TALE effectors found in {}", id);
        }
    });

    Ok(())
}
