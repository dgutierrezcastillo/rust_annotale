use anyhow::{Context, Result};
use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "Rename TALEs in File", long_about = None)]
struct Args {
    /// A tab-separated table containing the old name in the first column and the new name in the second column.
    #[arg(short = 'r', long = "rename-table", required = true)]
    rename_table: PathBuf,

    /// The input Genbank or GFF3 file that should be renamed.
    #[arg(short = 'i', long = "input", required = true)]
    input: PathBuf,

    /// The output directory, defaults to the current working directory (.)
    #[arg(long = "outdir", default_value = ".")]
    outdir: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Parse the TSV dictionary
    let mut rename_map: HashMap<String, String> = HashMap::new();
    let tsv_file = File::open(&args.rename_table)
        .with_context(|| format!("Failed to open rename table {:?}", args.rename_table))?;
    let tsv_reader = BufReader::new(tsv_file);

    for line in tsv_reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            rename_map.insert(parts[0].to_string(), parts[1].to_string());
        }
    }

    // Prepare input file
    let input_file = File::open(&args.input)
        .with_context(|| format!("Failed to open input file {:?}", args.input))?;
    let input_reader = BufReader::new(input_file);

    // Prepare output file
    if !args.outdir.exists() {
        std::fs::create_dir_all(&args.outdir)
            .with_context(|| format!("Failed to create output directory {:?}", args.outdir))?;
    }
    
    let file_name = args.input.file_name().context("Input file has no name")?;
    let output_path = args.outdir.join(file_name);
    
    let output_file = File::create(&output_path)
        .with_context(|| format!("Failed to create output file {:?}", output_path))?;
    let mut output_writer = BufWriter::new(output_file);

    // Stream and replace
    for line in input_reader.lines() {
        let mut current_line = line?;
        
        // Very basic replace; iterating over all keys is O(K) where K is dictionary size.
        // For a dictionary of TALE names, this is extremely fast.
        for (old_name, new_name) in &rename_map {
            if current_line.contains(old_name) {
                current_line = current_line.replace(old_name, new_name);
            }
        }
        
        writeln!(output_writer, "{}", current_line)?;
    }
    
    output_writer.flush()?;
    println!("Successfully wrote renamed file to {:?}", output_path);

    Ok(())
}
