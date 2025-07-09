//! This module provides functionality to extract alignment information from BAM records.

use super::Result;
use noodles::bam;
use noodles::sam::alignment::record::data::field::value::Value;
use noodles::sam::alignment::record::{
    cigar::op::Kind, Cigar, QualityScores as QualityScoresTrait,
};

#[derive(Debug, Clone, Default)]
pub struct Alignment {
    pub read_id: String,
    pub is_paired: bool,
    pub is_second_in_pair: bool,
    pub target_id: i32,
    pub is_unmapped: bool,
    pub is_primary: bool,
    pub alignment_score: i32,
    pub ref_start: i32,
    pub query_length: i32,
    pub query_covered: i32,
    pub ref_covered: i32,
    pub gap_compressed_seq_divergence: f32,
    pub expected_error_rate: f32,
    pub ref_id: Option<i32>,    // This requires the references_df to set
    pub ani_group: Option<i32>, // This requires the references_df to set
}

impl PartialEq for Alignment {
    fn eq(&self, other: &Self) -> bool {
        self.read_id == other.read_id
            && self.is_paired == other.is_paired
            && self.is_second_in_pair == other.is_second_in_pair
            && self.target_id == other.target_id
            && self.is_unmapped == other.is_unmapped
            && self.is_primary == other.is_primary
            && self.alignment_score == other.alignment_score
            && self.ref_start == other.ref_start
            && self.query_length == other.query_length
            && self.query_covered == other.query_covered
            && self.gap_compressed_seq_divergence == other.gap_compressed_seq_divergence
            && (self.expected_error_rate - other.expected_error_rate).abs() < f32::EPSILON
    }
}

/// Extra fields for debugging
///
/// Not needed for analysis
#[derive(Debug)]
pub struct AlignmentExtra {
    pub is_properly_segmented: bool, // Means that both reads in pair primary map to the same place
    pub is_supplementary: bool,
    pub mapping_quality: u8,
    pub number_mismatches: i32,
    pub chaining_score: i32,
    pub best_segment_score: i32,
    pub cigar_matches: i32,
    pub cigar_mismatches: i32,
    pub cigar_ins_sum: i32,
    pub cigar_del_sum: i32,
}

fn get_i32_tag(tags: &noodles::bam::record::Data, tag: &[u8; 2]) -> i32 {
    match tags.get(tag) {
        Some(Ok(val)) => {
            if let Some(i) = val.as_int() {
                return i as i32;
            }
            panic!(
                "Tag {} is not an integer {:?}",
                std::str::from_utf8(tag).unwrap(),
                val
            );
        }
        _ => 0,
    }
}

fn get_f32_tag(tags: &noodles::bam::record::Data, tag: &[u8; 2]) -> f32 {
    match tags.get(tag) {
        Some(Ok(val)) => match val {
            Value::Float(inner) => inner,
            _ => panic!(
                "Tag {} is not a float {:?}",
                std::str::from_utf8(tag).unwrap(),
                val
            ),
        },
        _ => -1.0,
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
struct CigarStats {
    matches: i32,
    mismatches: i32,
    insertions: i32,
    ins_sum: i32,
    deletions: i32,
    del_sum: i32,
    skips: i32,
    skip_sum: i32,
    soft_clipped: i32,
    hard_clipped: i32,
    query_aligned_bases: i32,
    ref_aligned_bases: i32,
    query_length: i32,
}

fn process_cigar<T: Cigar>(cigar: T) -> Result<CigarStats> {
    let mut stats = CigarStats {
        matches: 0,
        mismatches: 0,
        insertions: 0,
        ins_sum: 0,
        deletions: 0,
        del_sum: 0,
        skips: 0,
        skip_sum: 0,
        soft_clipped: 0,
        hard_clipped: 0,
        query_aligned_bases: 0,
        ref_aligned_bases: 0,
        query_length: 0,
    };

    if cigar.is_empty() {
        return Ok(stats);
    }

    for op in cigar.iter() {
        let op = op.expect("Every cigar operation should be valid");
        let len: i32 = op.len() as i32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch => {
                stats.matches += len;
            }
            Kind::SequenceMismatch => {
                stats.mismatches += len;
            }
            Kind::Insertion => {
                stats.insertions += 1;
                stats.ins_sum += len;
            }
            Kind::Deletion => {
                stats.deletions += 1;
                stats.del_sum += len;
            }
            Kind::SoftClip => {
                stats.soft_clipped += len;
            }
            Kind::HardClip => {
                stats.hard_clipped += len;
            }
            Kind::Skip => {
                stats.skips += 1;
                stats.skip_sum += len;
            }
            Kind::Pad => {
                panic!("Padding not supported");
            }
        }
    }

    stats.query_aligned_bases = stats.matches + stats.mismatches + stats.ins_sum;
    stats.ref_aligned_bases = stats.matches + stats.mismatches + stats.del_sum + stats.skip_sum;
    stats.query_length = stats.query_aligned_bases + stats.soft_clipped + stats.hard_clipped;

    return Ok(stats);
}

fn calculate_expected_error_rate<T: QualityScoresTrait>(quality_scores: T) -> f32 {
    // Essentially converting quality scores to probabilities and taking mean
    if quality_scores.is_empty() {
        return 0.0;
    }

    let sum: f32 = quality_scores
        .iter()
        .map(|q| {
            let q = q.expect("quality scores should be valid");
            // calculate error probability from quality score

            10.0_f32.powf(-(q as f32) / 10.0)
        })
        .sum();

    return sum / quality_scores.len() as f32;
}

/// Convert a BAM record to an `Alignment` and optionally `AlignmentExtra`
pub fn record_to_alignment_info(
    record: &bam::Record,
    debug: bool,
) -> (Alignment, Option<AlignmentExtra>) {
    let read_id = record
        .name()
        .map(|name| name.to_string())
        .unwrap_or("".to_string());
    let target_id = record
        .reference_sequence_id()
        .map(|id| id.unwrap() as i32)
        .unwrap_or(-1);
    let flags = record.flags();
    let is_unmapped = flags.is_unmapped();

    let is_primary = !(flags.is_secondary() || flags.is_supplementary());
    let is_paired = flags.is_segmented();

    if !is_paired && flags.is_first_segment() && flags.is_last_segment() {
        panic!(
            "Unexpected pairing. {} is segmented but is also first and last segment.",
            read_id
        );
    }
    let is_second_in_pair = flags.is_last_segment();
    if !is_paired & is_second_in_pair {
        panic!(
            "Unexpected pairing. {} is not paired but is also last segment.",
            read_id
        );
    }

    let tags = record.data();
    let alignment_score = get_i32_tag(&tags, b"AS");

    let start = match record.alignment_start() {
        Some(Ok(position)) => position.get() as i32,
        _ => 0,
    };

    let gap_compressed_seq_divergence = get_f32_tag(&tags, b"de");
    let expected_error_rate = calculate_expected_error_rate(record.quality_scores());

    let cigar_stats = process_cigar(record.cigar()).expect("cigar strings should all be valid");

    let aln = Alignment {
        read_id,
        is_paired,
        is_second_in_pair,
        target_id,
        is_unmapped,
        is_primary,
        alignment_score,
        ref_start: start,
        query_length: cigar_stats.query_length,
        query_covered: cigar_stats.query_aligned_bases,
        ref_covered: cigar_stats.ref_aligned_bases,
        gap_compressed_seq_divergence,
        expected_error_rate,
        ..Default::default()
    };

    let aln_extra = match (debug, is_unmapped) {
        (true, false) => Some(AlignmentExtra {
            is_properly_segmented: flags.is_properly_segmented(),
            is_supplementary: flags.is_supplementary(),

            mapping_quality: match record.mapping_quality() {
                Some(quality) => quality.get(),
                None => 0,
            },
            number_mismatches: get_i32_tag(&tags, b"NM"),
            chaining_score: get_i32_tag(&tags, b"s1"),
            best_segment_score: get_i32_tag(&tags, b"ms"),

            cigar_matches: cigar_stats.matches,
            cigar_mismatches: cigar_stats.mismatches,
            cigar_ins_sum: cigar_stats.ins_sum,
            cigar_del_sum: cigar_stats.del_sum,
        }),
        _ => None,
    };

    (aln, aln_extra)
}

#[cfg(test)]
mod tests {
    use std::fs::create_dir_all;

    use super::*;
    use crate::shared_test_functions::{bam_to_records, create_bam_from_lines, TEST_DIR};

    #[test]
    fn test_record_to_alignment_info() {
        let lines = [
            "@HD\tVN:1.0\tSO:unsorted",
            "@SQ\tSN:ref\tLN:1000",
            "1\t0\tref\t101\t60\t20M\t*\t0\t0\tGGTATCCGGTGTCGACCACA\t55555555555555555555\tNM:i:38\tms:i:359\tAS:i:346\tnn:i:0\ttp:A:P\tcm:i:12\ts1:i:113\ts2:i:113\tde:f:0.1\trl:i:30"
            ].iter().map(|s| s.to_string()).collect::<Vec<_>>();

        let expectation = [Alignment {
            read_id: "1".to_string(),
            is_paired: false,
            is_second_in_pair: false,
            target_id: 0,
            is_unmapped: false,
            is_primary: true,
            alignment_score: 346,
            ref_start: 101,
            query_length: 20,
            query_covered: 20,
            ref_covered: 20,
            gap_compressed_seq_divergence: 0.1,
            expected_error_rate: 0.01,
            ..Default::default()
        }];

        create_dir_all(format!("{}/alignments", TEST_DIR))
            .expect("Failed to create test directory");
        let bam_path = format!("{}/alignments/test_record_to_alignment_info.bam", TEST_DIR);
        create_bam_from_lines(&bam_path, lines).expect("Failed to create BAM file");

        let results = bam_to_records(&bam_path)
            .expect("Failed to read BAM records")
            .iter()
            .map(|record| record_to_alignment_info(record, false))
            .collect::<Vec<_>>();

        for (result, expected) in results.iter().zip(expectation.iter()) {
            assert_eq!(result.0, *expected, "Alignment info mismatch");
            assert!(result.1.is_none(), "Extra alignment info should be None");
        }
    }

    #[test]
    fn test_record_to_alignment_info_flags() {
        let lines = [
            "@HD\tVN:1.0\tSO:unsorted",
            "@SQ\tSN:ref\tLN:1000",
            "1\t0\tref\t101\t60\t20M\t*\t0\t0\t*\t*",
            "2\t69\t*\t*\t*\t*\t*\t0\t0\t*\t*",
            "3\t387\tref\t101\t60\t20M\t*\t0\t0\t*\t*",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();

        let expectation = vec![
            Alignment {
                is_paired: false,
                is_second_in_pair: false,
                is_unmapped: false,
                is_primary: true,
                ..Default::default()
            },
            Alignment {
                is_paired: true,
                is_second_in_pair: false,
                is_unmapped: true,
                is_primary: true,
                ..Default::default()
            },
            Alignment {
                is_paired: true,
                is_second_in_pair: true,
                is_unmapped: false,
                is_primary: false,
                ..Default::default()
            },
        ];

        let bam_path = format!(
            "{}/alignments/test_record_to_alignment_info_flags.bam",
            TEST_DIR
        );
        create_bam_from_lines(&bam_path, lines).expect("Failed to create BAM file");

        let results = bam_to_records(&bam_path)
            .expect("Failed to read BAM records")
            .iter()
            .map(|record| record_to_alignment_info(record, false).0)
            .collect::<Vec<_>>();

        for (result, expected) in results.iter().zip(expectation.iter()) {
            assert_eq!(result.is_paired, expected.is_paired, "is_paired mismatch");
            assert_eq!(
                result.is_second_in_pair, expected.is_second_in_pair,
                "is_second_in_pair mismatch"
            );
            assert_eq!(
                result.is_unmapped, expected.is_unmapped,
                "is_unmapped mismatch"
            );
            assert_eq!(
                result.is_primary, expected.is_primary,
                "is_primary mismatch"
            );
        }
    }

    #[test]
    fn test_cigar_stats() {
        let header_lines = vec![
            "@HD\tVN:1.0\tSO:unsorted".to_string(),
            "@SQ\tSN:ref\tLN:1000".to_string(),
        ];
        let lines = [
            "20M", "10H10M", "10S5M5=", "10X10M", "5M10I5M", "5M10D5M", "5M10N5M",
        ]
        .iter()
        .enumerate()
        .map(|(index, cigar)| format!("{}\t0\tref\t10\t60\t{}\t*\t0\t0\t*\t*", index + 1, cigar))
        .collect::<Vec<String>>();
        let lines = [header_lines, lines].concat();

        let expectation = vec![
            CigarStats {
                matches: 20,
                query_aligned_bases: 20,
                ref_aligned_bases: 20,
                query_length: 20,
                ..Default::default()
            },
            CigarStats {
                matches: 10,
                hard_clipped: 10,
                query_aligned_bases: 10,
                ref_aligned_bases: 10,
                query_length: 20,
                ..Default::default()
            },
            CigarStats {
                matches: 10,
                soft_clipped: 10,
                query_aligned_bases: 10,
                ref_aligned_bases: 10,
                query_length: 20,
                ..Default::default()
            },
            CigarStats {
                matches: 10,
                mismatches: 10,
                query_aligned_bases: 20,
                ref_aligned_bases: 20,
                query_length: 20,
                ..Default::default()
            },
            CigarStats {
                matches: 10,
                mismatches: 0,
                query_aligned_bases: 20,
                ref_aligned_bases: 10,
                ins_sum: 10,
                insertions: 1,
                query_length: 20,
                ..Default::default()
            },
            CigarStats {
                matches: 10,
                mismatches: 0,
                query_aligned_bases: 10,
                ref_aligned_bases: 20,
                del_sum: 10,
                deletions: 1,
                query_length: 10,
                ..Default::default()
            },
            CigarStats {
                matches: 10,
                mismatches: 0,
                query_aligned_bases: 10,
                ref_aligned_bases: 20,
                skip_sum: 10,
                skips: 1,
                query_length: 10,
                ..Default::default()
            },
        ];

        let bam_path = format!("{}/alignments/test_cigar_stats.bam", TEST_DIR);
        create_bam_from_lines(&bam_path, lines).expect("Failed to create BAM file");

        let results = bam_to_records(&bam_path)
            .expect("Failed to read BAM records")
            .iter()
            .map(|record| process_cigar(record.cigar()).expect("Failed to process CIGAR"))
            .collect::<Vec<_>>();

        for (result, expected) in results.iter().zip(expectation.iter()) {
            assert_eq!(result, expected);
        }
    }
}
