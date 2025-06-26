use cm_analyzer::analyze::{analyze_alignments, AnalyzeArgs};
use std::fs::{create_dir_all, read_to_string};

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
        depth_counts_file:
            "test_data/tb_intracellulare/with_ani_grouping/depth_counts.csv",
        expected_ref_tie_breaker_order_file:
            "test_data/tb_intracellulare/with_ani_grouping/ref_tie_breaker_order.csv",
        expected_stats_file: "test_data/tb_intracellulare/with_ani_grouping/stats.yaml",
        expected_summary_file:
            "test_data/tb_intracellulare/with_ani_grouping/alignment_summary.csv",
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

#[test]
fn test_all() {
    if let Some(threads) = THREADS {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }

    for test_set in TEST_SETS {
        let output_dir = test_set.output_root;
        create_dir_all(output_dir).expect("Failed to create output directory");
        let args = AnalyzeArgs {
            input_bam: test_set.input_bam.to_string(),
            contigs: test_set.contigs.to_string(),
            output_root: output_dir.to_string(),
            parameters: test_set.parameters.to_string(),
            threads: None,
            debug: false,
        };

        analyze_alignments(args).expect("Analysis failed");

        let round1_depth_equal = compare_files(
            test_set.depth_counts_file,
            &format!("{}depth_counts.csv", output_dir),
        );
        let ref_tie_breaker_order_equal = compare_files(
            test_set.expected_ref_tie_breaker_order_file,
            &format!("{}ref_tie_breaker_order.csv", output_dir),
        );
        let stats_equal = compare_files(
            test_set.expected_stats_file,
            &format!("{}stats.yaml", output_dir),
        );
        let depths_equal = compare_files(
            test_set.expected_summary_file,
            &format!("{}alignment_summary.csv", output_dir),
        );

        assert!(
            round1_depth_equal && ref_tie_breaker_order_equal && stats_equal && depths_equal,
            "Test set '{}' failed",
            test_set.name
        );
    }
}
