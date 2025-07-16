use tie_break::analyze::{analyze_alignments, AnalyzeArgs};
use clap::Parser;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = AnalyzeArgs::parse();

    println!("Using bam file: {}", args.input_bam);
    println!("Using contigs file: {}", args.contigs);
    println!("Using output root: {}", args.output_root);

    analyze_alignments(args)
}
