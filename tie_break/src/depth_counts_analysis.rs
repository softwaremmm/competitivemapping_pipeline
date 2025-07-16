//! This module provides functions to analyze depth counts.
//! It summarises this data to estimate of depth and coverage

#![allow(dead_code)]

use super::Result;
use polars::prelude::*;
use rand::{
    distr::{weighted::WeightedIndex, Distribution},
    rngs::StdRng,
    SeedableRng,
};

fn get_median(depth_counts: &[(u32, u32)]) -> f64 {
    let total_count: u32 = depth_counts.iter().map(|(_, v)| v).sum();
    if total_count == 0 {
        return f64::NAN;
    }
    let median = (total_count / 2) + 1; // handles odd correctly
    let mut cumulative_count = 0;
    for (depth, count) in depth_counts {
        cumulative_count += count;
        if cumulative_count >= median {
            return *depth as f64;
        }
    }

    f64::NAN
}

fn robust_mean(depth_counts: &[(u32, u32)]) -> f64 {
    let median = get_median(depth_counts);
    if median.is_nan() {
        return f64::NAN;
    }

    if median >= 3.0 {
        // Want to ignore high outliers, by setting alpha as the upper limit
        let mut alpha = 0;
        let mut current_term = (-median).exp();
        let mut cumulative_probability = current_term;
        while 1.0 - cumulative_probability > 10.0_f64.powf(-5.0) {
            alpha += 1;
            current_term *= median / alpha as f64;
            cumulative_probability += current_term;

            if alpha > 1000 {
                // Just set alpha very high
                alpha = u32::MAX;
                break;
            }
        }

        let mut running_sum = 0;
        let mut n = 0;
        for (depth, count) in depth_counts {
            if *depth > alpha {
                continue;
            }
            running_sum += depth * count;
            n += count;
        }
        return (running_sum as f64) / (n as f64);
    }

    let mut n_1 = 0;
    let mut n_2 = 0;
    let mut n_3 = 0;

    for (depth, count) in depth_counts {
        match depth {
            1 => n_1 = *count,
            2 => n_2 = *count,
            3 => n_3 = *count,
            _ => {}
        }
    }

    // println!("Using lambda adjustment. median: {}, n_1: {}, n_2: {}, n_3: {}", median, n_1, n_2, n_3);

    if n_1 >= 3 && n_2 >= 3 && (n_1 > n_2 || n_3 < 3) {
        return 2.0 * (n_2 as f64) / (n_1 as f64);
    } else if n_2 >= 3 && n_3 >= 3 {
        return 3.0 * (n_3 as f64) / (n_2 as f64);
    }

    return f64::NAN;
}

fn weighted_sample_fixed_size(rng: &mut StdRng, counts: &[u32], sample_size: u32) -> Vec<u32> {
    // Create a weighted distribution based on counts
    let dist = WeightedIndex::new(counts).unwrap();

    // Sample `sample_size` elements without replacement
    let mut sampled_counts = vec![0; counts.len()];
    for _ in 0..sample_size {
        let index = dist.sample(rng);
        sampled_counts[index] += 1;
    }

    sampled_counts
}

fn robust_mean_bootstrap(depth_counts: &[(u32, u32)]) -> [f64; 2] {
    let total_sites: u32 = depth_counts.iter().map(|(_, b)| b).sum();
    if total_sites < 50 {
        return [f64::NAN, f64::NAN];
    }

    let median = get_median(depth_counts);
    if median.is_nan() {
        return [f64::NAN, f64::NAN];
    }

    if median >= 3.0 {
        let robust_mean = robust_mean(depth_counts);
        return [robust_mean, robust_mean];
    }

    let depths: Vec<u32> = depth_counts.iter().map(|(depth, _)| *depth).collect();
    let counts: Vec<u32> = depth_counts.iter().map(|(_, count)| *count).collect();

    let sub_sample_size = total_sites / 50;
    let mut rng = StdRng::seed_from_u64(42);
    let mut sampling_results: Vec<f64> = Vec::new();
    for _i in 0..100 {
        let sampled_counts = weighted_sample_fixed_size(&mut rng, &counts, sub_sample_size);
        let sample_depth_counts: Vec<_> = depths
            .iter()
            .zip(sampled_counts.into_iter())
            .map(|(&depth, count)| (depth, count))
            .collect();
        let robust_mean = robust_mean(&sample_depth_counts);
        if robust_mean.is_finite() {
            sampling_results.push(robust_mean);
        }
    }
    if sampling_results.len() < 50 {
        return [f64::NAN, f64::NAN];
    }
    sampling_results.sort_by(|a, b| a.partial_cmp(b).expect("Have filtered out NANs"));
    let lower = sampling_results[sampling_results.len() / 10];
    let upper = sampling_results[sampling_results.len() * 9 / 10];
    [lower, upper]
}

fn calculate_expected_coverage(mean_depth: f64) -> f64 {
    1.0 - (-mean_depth).exp()
}

fn calculate_expected_coverage_col(
    col: Column,
) -> std::result::Result<Option<Column>, PolarsError> {
    let depths = col.f64()?;
    let expected_coverage = depths.apply(|x| x.map(|f| 100.0 * calculate_expected_coverage(f)));

    Ok(Some(Column::new(
        "expected_coverage".into(),
        expected_coverage,
    )))
}

fn columns_to_depth_counts(
    columns: &[Column],
) -> std::result::Result<Vec<(u32, u32)>, PolarsError> {
    let depths = columns[0]
        .u32()?
        .iter()
        .flatten()
        .collect::<Vec<_>>();
    let counts = columns[1]
        .u32()?
        .iter()
        .flatten()
        .collect::<Vec<_>>();

    Ok(depths.into_iter().zip(counts).collect())
}

fn agg_get_median(
    depths_counts: &mut [Column],
) -> std::result::Result<Option<Column>, PolarsError> {
    let depth_counts = columns_to_depth_counts(depths_counts)?;
    let median = get_median(&depth_counts);
    Ok(Some(Column::new("median".into(), vec![median])))
}

fn agg_get_robust_mean(
    depths_counts: &mut [Column],
) -> std::result::Result<Option<Column>, PolarsError> {
    let depth_counts = columns_to_depth_counts(depths_counts)?;
    let robust_mean = robust_mean(&depth_counts);
    Ok(Some(Column::new("robust_mean".into(), vec![robust_mean])))
}

fn agg_get_robust_mean_range(
    depths_counts: &mut [Column],
) -> std::result::Result<Option<Column>, PolarsError> {
    let mut depth_counts = columns_to_depth_counts(depths_counts)?;
    // sort depth_counts by depth for consistent output
    depth_counts.sort_by(|a, b| a.0.cmp(&b.0));
    let range = robust_mean_bootstrap(&depth_counts);
    Ok(Some(Column::new(
        "robust_mean_range".into(),
        vec![range[0], range[1]],
    )))
}

pub fn summarise_depth(depth_count_df: &DataFrame, reference_df: &DataFrame) -> Result<DataFrame> {
    let agg_df = depth_count_df
        .clone()
        .lazy()
        .sort(["ref_id", "depth_type", "depth"], Default::default())
        .group_by(["ref_id", "depth_type"])
        .agg([
            col("count").sum().alias("bases_covered"),
            (col("depth") * col("count")).sum().alias("total_bases"),
            apply_multiple(
                agg_get_median,
                &[col("depth"), col("count")],
                GetOutput::from_type(DataType::Float64),
                true,
            )
            .alias("median_nonzero_depth"),
            apply_multiple(
                agg_get_robust_mean,
                &[col("depth"), col("count")],
                GetOutput::from_type(DataType::Float64),
                true,
            )
            .alias("robust_depth_estimate"),
            apply_multiple(
                agg_get_robust_mean_range,
                &[col("depth"), col("count")],
                GetOutput::from_type(DataType::Float64),
                true,
            )
            .alias("robust_mean_range"),
        ])
        .with_columns([
            col("robust_mean_range")
                .list()
                .first()
                .alias("robust_mean_lower"),
            col("robust_mean_range")
                .list()
                .last()
                .alias("robust_mean_upper"),
        ])
        .drop([col("robust_mean_range")]);

    let df = agg_df
        .join(
            reference_df
                .clone()
                .lazy()
                .select([cols([
                    "ref_id",
                    "reference",
                    "ani_group",
                    "ref_length",
                    "species",
                ])])
                .unique(None, Default::default()),
            [col("ref_id")],
            [col("ref_id")],
            JoinArgs::new(JoinType::Left).with_coalesce(JoinCoalesce::CoalesceColumns),
        )
        .with_columns([
            (col("bases_covered").cast(DataType::Float32) * lit(100.0) / col("ref_length"))
                .round(3)
                .alias("coverage"),
            (col("total_bases").cast(DataType::Float32) / col("ref_length"))
                .round(3)
                .alias("mean_depth"),])
        .with_columns([
            col("mean_depth")
                .map(
                    calculate_expected_coverage_col,
                    GetOutput::from_type(DataType::Float64),
                )
                .round(2)
                .alias("simple_expected_coverage"),
            col("robust_depth_estimate")
                .map(
                    calculate_expected_coverage_col,
                    GetOutput::from_type(DataType::Float64),
                )
                .round(2)
                .alias("robust_expected_coverage"),
            col("robust_mean_lower")
                .map(
                    calculate_expected_coverage_col,
                    GetOutput::from_type(DataType::Float64),
                )
                .round(2)
                .alias("robust_expected_coverage_lower_bound"),
            col("robust_mean_upper")
                .map(
                    calculate_expected_coverage_col,
                    GetOutput::from_type(DataType::Float64),
                )
                .round(2)
                .alias("robust_expected_coverage_upper_bound"),])
        .with_columns([
            // combine the robust columns into single column MEAN (LOWER - UPPER)
            concat_str(
                [
                    col("robust_depth_estimate").round(2),
                    lit(" ("),
                    col("robust_mean_lower").round(2),
                    lit("-"),
                    col("robust_mean_upper").round(2),
                    lit(")"),
                ],
                "",
                true,
            )
            .alias("robust_depth_estimate"),
            concat_str(
                [
                    col("robust_expected_coverage"),
                    lit(" ("),
                    col("robust_expected_coverage_lower_bound"),
                    lit("-"),
                    col("robust_expected_coverage_upper_bound"),
                    lit(")"),
                ],
                "",
                true,
            )
            .alias("robust_expected_coverage"),
        ])
        .drop([cols([
            "bases_covered",
            "total_bases",
            "robust_mean_lower",
            "robust_expected_coverage_lower_bound",
            "robust_mean_upper",
            "robust_expected_coverage_upper_bound",
        ])]);

    let front_cols = [
        "ref_id",
        "reference",
        "ani_group",
        "species",
        "ref_length",
        "depth_type",
        "coverage",
        "mean_depth",
        "median_nonzero_depth",
        "simple_expected_coverage",
        "robust_depth_estimate",
        "robust_expected_coverage",
    ];

    let final_df = df
        .select([cols(front_cols), all().exclude(front_cols)])
        .sort(["species", "depth_type"], SortMultipleOptions::default())
        .collect()?;

    Ok(final_df)
}

pub fn order_best_refs(summarised_df: &DataFrame, min_coverage_pc: f64) -> Result<DataFrame> {
    // Order the df by best_covered_pc and convert references column to vector of strings
    let ordered_refs = summarised_df
        .clone()
        .lazy()
        .filter(col("depth_type").eq(lit("unique")))
        .sort(
            ["coverage"],
            SortMultipleOptions::default().with_order_descending(true),
        )
        .filter(col("coverage").gt(min_coverage_pc))
        .select([cols(["ref_id", "reference", "ani_group"])])
        .collect()?;
    Ok(ordered_refs)
}

pub fn count_alns(alns_csv_file: &str, reference_df: &DataFrame) -> Result<DataFrame> {

    // schema used to force some types
    let schema = Schema::from_iter(vec![
        Field::new("target_id".into(), DataType::UInt32),
    ]);

    let alns_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_infer_schema_length(None)
        .with_schema_overwrite(Some(std::sync::Arc::new(schema)))
        .try_into_reader_with_file_path(Some(alns_csv_file.into()))?
        .finish()?;

    let aln_counts_df = alns_df
        .lazy()
        .group_by(["target_id"])
        .agg([len().alias("count")])
        .join(
            reference_df
                .clone()
                .lazy()
                .select([cols([
                    "target_id",
                    "ref_id",
                ])])
                .unique(None, Default::default()),
            [col("target_id")],
            [col("target_id")],
            JoinArgs::new(JoinType::Left).with_coalesce(JoinCoalesce::CoalesceColumns),
        )
        .group_by(["ref_id"])
        .agg([
            col("count").sum().alias("reads"),
        ])
        .sort(["ref_id"], Default::default())
        .collect()?;

    return Ok(aln_counts_df);

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        save_csv,
        shared_test_functions::{compare_files, create_parent_dir, TEST_DATA, TEST_DIR},
    };

    struct TestSet {
        name: &'static str,
        depth_counts_file: &'static str,
        refs_file: &'static str,
        expectation_file: &'static str,
    }
    const TEST_SETS: &[TestSet] = &[
        TestSet {
            name: "simple",
            depth_counts_file: "simple/depth_counts.csv",
            refs_file: "simple/references.csv",
            expectation_file: "simple/summarised_depths.csv",
        },
        TestSet {
            name: "full",
            depth_counts_file: "full/depth_counts.csv",
            refs_file: "full/references.csv",
            expectation_file: "full/summarised_depths.csv",
        },
        TestSet {
            name: "full",
            depth_counts_file: "robust/depth_counts.csv",
            refs_file: "robust/references.csv",
            expectation_file: "robust/summarised_depths.csv",
        },
    ];

    fn get_abs_input_path(path: &str) -> String {
        format!("{}/{}/{}", TEST_DATA, "summarise_depth", path)
    }
    fn get_abs_output_path(path: &str) -> String {
        format!("{}/{}/{}", TEST_DIR, "summarise_depth", path)
    }

    fn read_csv(file_path: &str) -> DataFrame {
        CsvReadOptions::default()
            .with_infer_schema_length(None)
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(file_path.into()))
            .unwrap()
            .finish()
            .unwrap()
    }

    #[test]
    fn test_summarise_depth_all() {
        for test_set in TEST_SETS {
            let depth_counts_file = get_abs_input_path(test_set.depth_counts_file);
            let reference_file = get_abs_input_path(test_set.refs_file);
            let expected_summary_file = get_abs_input_path(test_set.expectation_file);

            let depth_counts_df: DataFrame = read_csv(&depth_counts_file)
                .lazy()
                .with_columns([
                    col("depth").cast(DataType::UInt32),
                    col("count").cast(DataType::UInt32),
                ])
                .collect()
                .unwrap();
            let reference_df: DataFrame = read_csv(&reference_file);

            let summarised_df = summarise_depth(&depth_counts_df, &reference_df).unwrap();
            // save
            let output_path = get_abs_output_path(test_set.expectation_file);
            create_parent_dir(&output_path);
            save_csv(&summarised_df, output_path.clone()).unwrap();

            assert!(compare_files(&expected_summary_file, &output_path));
        }
    }

    #[test]
    fn test_count_alns() {
        let alns_file = get_abs_input_path("count_alns/alns.csv");
        let reference_file = get_abs_input_path("count_alns/references.csv");
        let expected_counts_file = get_abs_input_path("count_alns/aln_counts.csv");

        let reference_df: DataFrame = read_csv(&reference_file);
        let aln_counts_df = count_alns(&alns_file, &reference_df).unwrap();

        // save
        let output_path = get_abs_output_path("count_alns/aln_counts.csv");
        create_parent_dir(&output_path);
        save_csv(&aln_counts_df, output_path.clone()).unwrap();

        assert!(compare_files(&expected_counts_file, &output_path));
    }
}
