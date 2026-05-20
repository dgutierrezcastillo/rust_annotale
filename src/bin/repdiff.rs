use anyhow::{Context, Result};
use bio::io::fasta;
use clap::Parser;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "TALE Repeat Differences", long_about = None)]
struct Args {
    /// TALE sequences, complete DNA or AS sequences
    #[arg(short = 't', long = "tale-sequences", required = true)]
    tale_sequences: PathBuf,

    /// The output directory
    #[arg(long = "outdir", default_value = ".")]
    outdir: PathBuf,
}

// Reuse translation for repdiff
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

fn extract_repeats_aa(sequence: &[u8]) -> Vec<Vec<u8>> {
    let mut repeats = Vec::new();
    
    // Check if it's already AA by looking for typical AA characters that don't exist in pure DNA
    let is_aa = sequence.iter().any(|&c| c == b'L' || c == b'P' || c == b'Q' || c == b'R');
    
    let aa_seq = if is_aa { sequence.to_vec() } else { translate(sequence) };
    
    let mut start_idx = None;
    for i in 0..aa_seq.len().saturating_sub(3) {
        if aa_seq[i] == b'L' && aa_seq[i+1] == b'T' && aa_seq[i+2] == b'P' {
            start_idx = Some(i);
            break;
        }
    }

    if let Some(start) = start_idx {
        let mut curr = start;
        while curr + 34 <= aa_seq.len() {
            let repeat_aa = &aa_seq[curr..curr+34];
            repeats.push(repeat_aa.to_vec());
            curr += 34;
        }
    }

    repeats
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !args.outdir.exists() {
        std::fs::create_dir_all(&args.outdir)?;
    }

    let reader = fasta::Reader::from_file(&args.tale_sequences)?;
    let mut tale_data = Vec::new();
    
    for record in reader.records() {
        let rec = record?;
        let repeats = extract_repeats_aa(rec.seq());
        tale_data.push((rec.id().to_string(), repeats));
    }

    let out_path = args.outdir.join("repdiff_matrix.tsv");
    let mut writer = BufWriter::new(File::create(&out_path)?);
    
    // We compute differences between all pairs of TALEs
    // Using Rayon for parallel distance matrix calculation
    let results: Vec<_> = tale_data.par_iter().map(|(id_a, repeats_a)| {
        let mut row_results = Vec::new();
        for (id_b, repeats_b) in &tale_data {
            let mut total_distance = 0;
            // Compare repeat pairs up to the length of the shorter TALE
            let min_len = std::cmp::min(repeats_a.len(), repeats_b.len());
            for i in 0..min_len {
                let dist = strsim::levenshtein(
                    std::str::from_utf8(&repeats_a[i]).unwrap_or(""),
                    std::str::from_utf8(&repeats_b[i]).unwrap_or("")
                );
                total_distance += dist;
            }
            row_results.push((id_b.clone(), total_distance));
        }
        (id_a.clone(), row_results)
    }).collect();

    // Write TSV Matrix
    write!(writer, "TALE")?;
    for (id_b, _) in &tale_data {
        write!(writer, "\t{}", id_b)?;
    }
    writeln!(writer)?;

    for (id_a, row_results) in results {
        write!(writer, "{}", id_a)?;
        for (_, dist) in row_results {
            write!(writer, "\t{}", dist)?;
        }
        writeln!(writer)?;
    }

    println!("Repeat differences calculated. Output written to {:?}", out_path);
    Ok(())
}
