use anyhow::{Context, Result};
use bio::io::fasta;
use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "TALE Frameshift and Truncation Scanner", long_about = None)]
struct Args {
    /// Input FASTA file (can be DNA sequences or RVD sequences)
    #[arg(short = 'i', long = "input", required = true)]
    input: PathBuf,

    /// Output file path (optional, writes to stdout if not provided)
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,
}

// Minimal translation map for TALE repeats
fn translate(dna: &[u8]) -> Vec<u8> {
    let mut aa = Vec::new();
    for chunk in dna.chunks_exact(3) {
        let a = match chunk {
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
            b"TGT" | b"TGC" => b'C',
            b"GAT" | b"GAC" => b'D',
            b"GAA" | b"GAG" => b'E',
            b"TTT" | b"TTC" => b'F',
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
            b"CAT" | b"CAC" => b'H',
            b"ATT" | b"ATC" | b"ATA" => b'I',
            b"AAA" | b"AAG" => b'K',
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
            b"ATG" => b'M',
            b"AAT" | b"AAC" => b'N',
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
            b"CAA" | b"CAG" => b'Q',
            b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
            b"TGG" => b'W',
            b"TAT" | b"TAC" => b'Y',
            b"TAA" | b"TAG" | b"TGA" => b'*',
            _ => b'X',
        };
        aa.push(a);
    }
    aa
}

fn dna_to_rvds(seq: &[u8]) -> Vec<String> {
    let mut rvds = Vec::new();
    let aa_seq = translate(seq);
    let mut start_idx = None;
    for i in 0..aa_seq.len().saturating_sub(3) {
        if aa_seq[i] == b'L' && aa_seq[i+1] == b'T' && aa_seq[i+2] == b'P' {
            start_idx = Some(i);
            break;
        }
    }
    if let Some(start) = start_idx {
        let mut curr = start * 3;
        while curr + 102 <= seq.len() {
            let repeat_dna = &seq[curr..curr+102];
            let repeat_aa = translate(repeat_dna);
            if repeat_aa.len() >= 14 {
                let rvd = format!("{}{}", repeat_aa[12] as char, repeat_aa[13] as char);
                rvds.push(rvd);
            }
            curr += 102;
        }
    }
    rvds
}

fn is_dna(seq: &[u8]) -> bool {
    if seq.contains(&b'-') || seq.contains(&b',') {
        return false;
    }
    seq.iter().all(|&b| match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' | b'N' | b'U' => true,
        _ => false,
    })
}

fn parse_rvd_sequence(seq_str: &str) -> Vec<String> {
    let cleaned: String = seq_str.chars().filter(|c| !c.is_whitespace()).collect();
    if cleaned.contains('-') {
        cleaned.split('-').map(|s| s.to_string()).collect()
    } else if cleaned.contains(',') {
        cleaned.split(',').map(|s| s.to_string()).collect()
    } else {
        // Chunk into pairs of 2 characters
        let chars: Vec<char> = cleaned.chars().collect();
        let mut rvds = Vec::new();
        let mut i = 0;
        while i < chars.len() {
            if i + 1 < chars.len() {
                rvds.push(format!("{}{}", chars[i], chars[i+1]));
                i += 2;
            } else {
                rvds.push(chars[i].to_string());
                i += 1;
            }
        }
        rvds
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let reader = fasta::Reader::from_file(&args.input)
        .with_context(|| format!("Failed to read input FASTA file {:?}", args.input))?;

    let mut tales: Vec<(String, Vec<String>)> = Vec::new();

    for record in reader.records() {
        let rec = record?;
        let id = rec.id().to_string();
        let seq = rec.seq();
        if is_dna(seq) {
            let rvds = dna_to_rvds(seq);
            tales.push((id, rvds));
        } else {
            let seq_str = String::from_utf8_lossy(seq);
            let rvds = parse_rvd_sequence(&seq_str);
            tales.push((id, rvds));
        }
    }

    // Set up output target
    let mut writer: Box<dyn Write> = match args.output {
        Some(ref path) => {
            let file = File::create(path)
                .with_context(|| format!("Failed to create output file {:?}", path))?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(std::io::stdout()),
    };

    // First scan: Internal Size/Frameshifts
    // Compares all pairs (i, j) where i < j
    for i in 0..tales.len() {
        for j in (i + 1)..tales.len() {
            let (id1, rvds1) = &tales[i];
            let (id2, rvds2) = &tales[j];

            if rvds1.len() != rvds2.len() && rvds1.len() >= 4 && rvds2.len() >= 4 {
                let prefix_match = rvds1[0..4] == rvds2[0..4];
                let suffix_match = rvds1[rvds1.len() - 4..] == rvds2[rvds2.len() - 4..];

                if prefix_match && suffix_match {
                    writeln!(writer, "{}", id1)?;
                    writeln!(writer, "{}", rvds1.join("-"))?;
                    writeln!(writer, "{}", id2)?;
                    writeln!(writer, "{}", rvds2.join("-"))?;
                    writeln!(writer, "+++++++++++++++++++++++++++++++++++++++++++++++++")?;
                }
            }
        }
    }

    writeln!(writer, "#######################################################")?;

    // Second scan: Truncations
    // Compares all pairs (i, j) where i < j
    for i in 0..tales.len() {
        for j in (i + 1)..tales.len() {
            let (id1, rvds1) = &tales[i];
            let (id2, rvds2) = &tales[j];

            let (long_id, long_rvds, short_id, short_rvds) = if rvds1.len() > rvds2.len() {
                (id1, rvds1, id2, rvds2)
            } else if rvds1.len() < rvds2.len() {
                (id2, rvds2, id1, rvds1)
            } else {
                continue;
            };

            if short_rvds.is_empty() {
                continue;
            }

            let prefix_match = long_rvds[0..short_rvds.len()] == *short_rvds;
            let suffix_match = long_rvds[long_rvds.len() - short_rvds.len()..] == *short_rvds;

            if prefix_match || suffix_match {
                writeln!(writer, "{}", long_id)?;
                writeln!(writer, "{}", long_rvds.join("-"))?;
                writeln!(writer, "{}", short_id)?;
                writeln!(writer, "{}", short_rvds.join("-"))?;
                writeln!(writer, "+++++++++++++++++++++++++++++++++++++++++++++++++")?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}
