use std::fs::{create_dir_all, read_to_string};
use tie_break::analyze::{analyze_alignments, AnalyzeArgs};

#[derive(Debug, Clone)]
struct TestSet {
    name: &'static str,
    input_bam: &'static str,
    contigs: &'static str,
    output_root: &'static str,
    parameters: &'static str,
    depth_counts_file: &'static str,
    expected_ref_tie_breaker_order_file: &'static str,
    expected_stats_file: &'static str,
    expected_summary_file: &'static str,
}

const TEST_SETS: &[TestSet] = &[
    TestSet {
        name: "tb_small",
        input_bam: "test_data/tb_small/aln.bam",
        contigs: "test_data/tb_small/contigs.csv",
        output_root: "tests/test_produced_files/tb_small/",
        parameters: "test_data/tb_small/params.yaml",
        depth_counts_file: "test_data/tb_small/depth_counts.csv",
        expected_ref_tie_breaker_order_file: "test_data/tb_small/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/tb_small/stats.yaml",
        expected_summary_file: "test_data/tb_small/alignment_summary.csv",
    },
    TestSet {
        name: "no_mapping",
        input_bam: "test_data/no_mapping/aln.bam",
        contigs: "test_data/no_mapping/contigs.csv",
        output_root: "tests/test_produced_files/no_mapping/",
        parameters: "test_data/no_mapping/params.yaml",
        depth_counts_file: "test_data/no_mapping/depth_counts.csv",
        expected_ref_tie_breaker_order_file: "test_data/no_mapping/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/no_mapping/stats.yaml",
        expected_summary_file: "test_data/no_mapping/alignment_summary.csv",
    },
    TestSet {
        name: "no_good_targets",
        input_bam: "test_data/tb_intracellulare/aln.bam",
        contigs: "test_data/tb_intracellulare/contigs.csv",
        output_root: "tests/test_produced_files/no_good_targets/",
        parameters: "test_data/tb_intracellulare/no_good_targets/params.yaml",
        depth_counts_file: "test_data/tb_intracellulare/no_good_targets/depth_counts.csv",
        expected_ref_tie_breaker_order_file:
            "test_data/tb_intracellulare/no_good_targets/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/tb_intracellulare/no_good_targets/stats.yaml",
        expected_summary_file: "test_data/tb_intracellulare/no_good_targets/alignment_summary.csv",
    },
    TestSet {
        name: "tb_intracellulare",
        input_bam: "test_data/tb_intracellulare/aln.bam",
        contigs: "test_data/tb_intracellulare/contigs.csv",
        output_root: "tests/test_produced_files/tb_intracellulare/",
        parameters: "test_data/tb_intracellulare/params.yaml",
        depth_counts_file: "test_data/tb_intracellulare/depth_counts.csv",
        expected_ref_tie_breaker_order_file:
            "test_data/tb_intracellulare/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/tb_intracellulare/stats.yaml",
        expected_summary_file: "test_data/tb_intracellulare/alignment_summary.csv",
    },
    TestSet {
        name: "with_ani_grouping",
        input_bam: "test_data/tb_intracellulare/aln.bam",
        contigs: "test_data/tb_intracellulare/with_ani_grouping/contigs.csv",
        output_root: "tests/test_produced_files/tb_intracellulare_with_ani_grouping/",
        parameters: "test_data/tb_intracellulare/params.yaml",
        depth_counts_file: "test_data/tb_intracellulare/with_ani_grouping/depth_counts.csv",
        expected_ref_tie_breaker_order_file:
            "test_data/tb_intracellulare/with_ani_grouping/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/tb_intracellulare/with_ani_grouping/stats.yaml",
        expected_summary_file:
            "test_data/tb_intracellulare/with_ani_grouping/alignment_summary.csv",
    },
    TestSet {
        name: "with_partial_ani_grouping",
        input_bam: "test_data/tb_intracellulare/aln.bam",
        contigs: "test_data/tb_intracellulare/with_partial_ani_grouping/contigs.csv",
        output_root: "tests/test_produced_files/tb_intracellulare_with_partial_ani_grouping/",
        parameters: "test_data/tb_intracellulare/params.yaml",
        depth_counts_file: "test_data/tb_intracellulare/with_partial_ani_grouping/depth_counts.csv",
        expected_ref_tie_breaker_order_file:
            "test_data/tb_intracellulare/with_partial_ani_grouping/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/tb_intracellulare/with_partial_ani_grouping/stats.yaml",
        expected_summary_file:
            "test_data/tb_intracellulare/with_partial_ani_grouping/alignment_summary.csv",
    },
    TestSet {
        name: "fortuitum",
        input_bam: "test_data/fortuitum/aln.bam",
        contigs: "test_data/fortuitum/contigs.csv",
        output_root: "tests/test_produced_files/fortuitum/",
        parameters: "test_data/fortuitum/params.yaml",
        depth_counts_file: "test_data/fortuitum/depth_counts.csv",
        expected_ref_tie_breaker_order_file: "test_data/fortuitum/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/fortuitum/stats.yaml",
        expected_summary_file: "test_data/fortuitum/alignment_summary.csv",
    },
];

const THREADS: Option<usize> = Some(20);

fn update_expectations() -> bool {
    std::env::var("UPDATE_EXPECTATIONS")
        .map(|val| val == "1" || val.to_lowercase() == "true")
        .unwrap_or(false)
}

fn compare_files(expected: &str, result: &str) -> bool {
    let expected_content = match read_to_string(expected) {
        Ok(content) => content,
        Err(_) => {
            println!("Failed to read expected file '{expected}'");
            "".to_string()
        }
    };
    let result_content = read_to_string(result).unwrap();
    let equal = expected_content == result_content;

    if !equal {
        println!("Files differ:\nExpected: {expected}\nResult: {result}");

        if update_expectations() {
            std::fs::write(expected, &result_content).expect("Failed to update expected file");
            println!("Updated expected file: {expected}");

            return true; // Consider it equal after updating
        }
    }

    return equal;
}

#[test]
fn test_all() {
    for test_set in TEST_SETS {
        let output_dir = test_set.output_root;
        create_dir_all(output_dir).expect("Failed to create output directory");
        let args = AnalyzeArgs {
            input_bam: test_set.input_bam.to_string(),
            contigs: test_set.contigs.to_string(),
            output_root: output_dir.to_string(),
            parameters: test_set.parameters.to_string(),
            threads: THREADS,
            debug: false,
        };

        analyze_alignments(args).expect("Analysis failed");

        let round1_depth_equal = compare_files(
            test_set.depth_counts_file,
            &format!("{output_dir}depth_counts.csv"),
        );
        let ref_tie_breaker_order_equal = compare_files(
            test_set.expected_ref_tie_breaker_order_file,
            &format!("{output_dir}ref_tie_breaker_order.csv"),
        );
        let stats_equal = compare_files(
            test_set.expected_stats_file,
            &format!("{output_dir}stats.yaml"),
        );
        let depths_equal = compare_files(
            test_set.expected_summary_file,
            &format!("{output_dir}alignment_summary.csv"),
        );

        assert!(
            round1_depth_equal && ref_tie_breaker_order_equal && stats_equal && depths_equal,
            "Test set '{}' failed",
            test_set.name
        );
    }
}

#[test]
fn test_debug() {
    let test_set = TEST_SETS[0].clone();

    let output_dir = &format!("{}debug/", test_set.output_root);
    create_dir_all(output_dir).expect("Failed to create output directory");
    let args = AnalyzeArgs {
        input_bam: test_set.input_bam.to_string(),
        contigs: test_set.contigs.to_string(),
        output_root: output_dir.to_string(),
        parameters: test_set.parameters.to_string(),
        threads: Some(1),
        debug: true,
    };

    analyze_alignments(args).expect("Analysis failed");

    let round1_depth_equal = compare_files(
        test_set.depth_counts_file,
        &format!("{output_dir}depth_counts.csv"),
    );
    let ref_tie_breaker_order_equal = compare_files(
        test_set.expected_ref_tie_breaker_order_file,
        &format!("{output_dir}ref_tie_breaker_order.csv"),
    );
    let stats_equal = compare_files(
        test_set.expected_stats_file,
        &format!("{output_dir}stats.yaml"),
    );
    let depths_equal = compare_files(
        test_set.expected_summary_file,
        &format!("{output_dir}alignment_summary.csv"),
    );
    let debug_equal = compare_files(
        &test_set
            .depth_counts_file
            .replace("depth_counts", "debug_alns_round_2"),
        &format!("{output_dir}debug_alns_round_2.csv"),
    );

    assert!(
        round1_depth_equal
            && ref_tie_breaker_order_equal
            && stats_equal
            && depths_equal
            && debug_equal,
        "Test set '{}' failed",
        test_set.name
    );
}
