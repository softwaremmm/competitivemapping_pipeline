use crate::{make_bam_index, sort_bam, Result};
use noodles::bam;
use std::fs::read_to_string;
use std::io::Write;
use std::path::Path;
use std::{fs::{File, create_dir_all}, process::Command};

use noodles::bgzf::MultithreadedReader;
use std::num::NonZero;

pub const TEST_DIR: &str = "tests/test_produced_files/";
pub const TEST_DATA: &str = "test_data/";

fn update_expectations() -> bool {
    std::env::var("UPDATE_EXPECTATIONS")
        .map(|val| val == "1" || val.to_lowercase() == "true")
        .unwrap_or(false)
}

pub fn compare_files(expected: &str, result: &str) -> bool {
    let expected_content = match read_to_string(expected) {
        Ok(content) => content,
        Err(_) => {
            println!("Failed to read expected file '{}'", expected);
            "".to_string()
        }
    };
    let result_content = read_to_string(result).unwrap();
    let equal = expected_content == result_content;

    if !equal {
        println!("Files differ:\nExpected: {}\nResult: {}", expected, result);

        if update_expectations() {
            std::fs::write(expected, &result_content).expect("Failed to update expected file");
            println!("Updated expected file: {}", expected);

            return true; // Consider it equal after updating
        }
    }

    return equal;
}

pub fn create_parent_dir<P: AsRef<Path>>(file_path: P) {
    if let Some(parent) = file_path.as_ref().parent() {
        create_dir_all(parent).expect("Failed to create parent directory");
    }
}

pub fn create_bam_from_lines(file_path: &str, lines: Vec<String>) -> Result<()> {
    let sam_path = format!("{}.sam", file_path);

    // create folder if it doesn't exist
    create_parent_dir(file_path);

    // Create a SAM file with the specified reference length and alignments
    let mut sam_file = File::create(&sam_path)?;
    for line in lines {
        writeln!(sam_file, "{}", line)?;
    }

    // Close the SAM file
    drop(sam_file);

    let unsorted_path = format!("{}.unsorted.bam", file_path);

    // Convert SAM to BAM
    Command::new("samtools")
        .args(["view", "-bS", &sam_path, "-o", &unsorted_path])
        .output()
        .expect("Failed to convert SAM to BAM");

    // Sort the BAM file
    sort_bam(&unsorted_path, file_path, None)?;

    make_bam_index(file_path)?;
    // Clean up the intermediate files
    std::fs::remove_file(&sam_path)?;
    std::fs::remove_file(&unsorted_path)?;

    return Ok(());
}

pub fn create_bam(
    file_path: &str,
    ref_length: usize,
    alns: Vec<(usize, usize, usize, String)>,
) -> Result<()> {
    // alignments are tuples of (query, flag, start, cigar)

    let header_lines = vec![
        "@HD\tVN:1.0\tSO:unsorted".to_string(),
        format!("@SQ\tSN:ref\tLN:{}", ref_length),
    ];

    let lines = alns
        .iter()
        .map(|(query_id, flag, start, cigar)| {
            format!(
                "{}\t{}\tref\t{}\t60\t{}\t*\t0\t0\t*\t*",
                query_id, flag, start, cigar
            )
        })
        .collect::<Vec<String>>();

    let lines = [header_lines, lines].concat();

    create_bam_from_lines(file_path, lines)
}

pub fn bam_to_records(file_path: &str) -> Result<Vec<bam::Record>> {
    let decoder = MultithreadedReader::with_worker_count(
        NonZero::new(1).unwrap(),
        File::open(file_path).unwrap(),
    );
    let mut reader = bam::io::Reader::from(decoder);
    let _header = reader.read_header().expect("Failed to read BAM header");

    let results: Vec<_> = reader
        .records()
        .map(|result| {
            result.expect("Failed to read record")
        })
        .collect();

    Ok(results)
}
