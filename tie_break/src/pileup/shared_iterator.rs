//! Module for counting read depths from alignment files
//!
//! Currently set to use alignments in CSV format
//! But can be adapted to use BAM files fairly easily
//!
//! Used to be multithreaded with BAM files,
//! but that is much slower with csvs due to repeated reading and filtering csv file
//! And often majority of alignments are to the top target_id anyway

use crate::filter_reads::filter_counts::Signals;
use crate::Result;
use polars::prelude::*;
use std::collections::HashMap;
use std::time::SystemTime;

use super::csv_reader::CsvPileupIterator;

#[macro_export]
macro_rules! hashmap_to_dataframe {
    ($hashmap:expr, $key_name:expr, $value_name:expr) => {{
        use polars::prelude::*;
        let (keys, values): (Vec<_>, Vec<_>) = $hashmap.into_iter().unzip();
        DataFrame::new(vec![
            Series::new($key_name, keys).into(),
            Series::new($value_name, values).into(),
        ])
    }};
}

fn depth_count_hashmap_to_dataframe(depth_counts: HashMap<(u32, u32), u32>) -> Result<DataFrame> {
    let (tid_depths, counts): (Vec<_>, Vec<_>) = depth_counts.into_iter().unzip();
    let (tids, depths): (Vec<_>, Vec<_>) = tid_depths.into_iter().unzip();
    Ok(DataFrame::new(vec![
        Series::new("target_id".into(), tids).into(),
        Series::new("depth".into(), depths).into(),
        Series::new("count".into(), counts).into(),
    ])?)
}

pub fn get_depth_counts_round1(
    reference_df: &DataFrame,
    alns_path: &str,
) -> Result<DataFrame> {
    let now = SystemTime::now();

    let unique_df = get_depth_counts(reference_df, alns_path, Some(vec![Signals::Unique]))?;
    let winner_df = get_depth_counts(reference_df, alns_path, Some(vec![Signals::Unique, Signals::Winner]))?;

    let combined = concat(
        [
            unique_df
                .lazy()
                .with_column(lit("unique").alias("depth_type")),
            winner_df
                .lazy()
                .with_column(lit("winner").alias("depth_type")),
        ],
        UnionArgs::default(),
    )?
    .filter(col("depth").gt(0))
    .sort(
            ["ref_id", "depth", "depth_type"],
            SortMultipleOptions::default(),
        )
    .collect()?;

    println!("Reading took {:?} overall", now.elapsed().unwrap());
    Ok(combined)
}


pub fn get_depth_counts(reference_df: &DataFrame, alns_paths: &str, signals: Option<Vec<Signals>>) -> Result<DataFrame> {
    let now = SystemTime::now();

    let depth_iter = CsvPileupIterator::from_path(alns_paths, None, signals)
        .expect("Failed to create iterator from path");

    let mut depth_counts: HashMap<(u32, u32), u32> = HashMap::new();

    for (tid, _pos, depth) in depth_iter {
        let count = depth_counts.entry((tid, depth)).or_insert(0);
        *count += 1;
    }

    let df = depth_count_hashmap_to_dataframe(depth_counts)?;

    let final_df = df
        .lazy()
        .join(
            reference_df
                .clone()
                .lazy()
                .select([cols(["target_id", "ref_id"])]),
            [col("target_id")],
            [col("target_id")],
            JoinArgs::new(JoinType::Left).with_coalesce(JoinCoalesce::CoalesceColumns),
        )
        .group_by(["ref_id", "depth"])
        .agg([col("count").sum().alias("count")])
        .with_column(lit("final").alias("depth_type"))
        .sort(
            ["ref_id", "depth", "depth_type"],
            SortMultipleOptions::default(),
        )
        .select([cols(["ref_id", "depth_type", "depth", "count"])])
        .collect()?;

    println!("Reading took {:?} overall", now.elapsed().unwrap());
    println!("Depth: {final_df:?}");
    Ok(final_df)
}

#[cfg(test)]
mod tests {
    use std::fs::create_dir_all;

    use super::*;
    use crate::save_csv;
    use crate::shared_test_functions::{compare_files, TEST_DATA, TEST_DIR};

    fn reference_df() -> DataFrame {
        let references_path = format!("{}/{}", TEST_DATA, "pileup/references.csv");

        // schema used to force some types
        let schema = Schema::from_iter(vec![
            Field::new("target_id".into(), DataType::UInt32),
            Field::new("ref_id".into(), DataType::UInt32),
        ]);

        let df: DataFrame = CsvReadOptions::default()
            .with_infer_schema_length(None)
            .with_has_header(true)
            .with_schema_overwrite(Some(std::sync::Arc::new(schema)))
            .try_into_reader_with_file_path(Some(references_path.into()))
            .unwrap()
            .finish()
            .unwrap();

        return df;
    }

    #[test]
    fn test_hashmap_to_dataframe() {
        let mut hashmap: HashMap<u32, u32> = HashMap::new();
        hashmap.insert(1, 10);
        hashmap.insert(2, 20);
        hashmap.insert(3, 30);

        let df = hashmap_to_dataframe!(hashmap, "key".into(), "value".into()).unwrap();

        println!("{df:?}");

        assert_eq!(df.shape(), (3, 2));
        let key: u32 = df.column("key").unwrap().u32().unwrap().get(0).unwrap();
        let value: u32 = df.column("value").unwrap().u32().unwrap().get(0).unwrap();
        assert!(key * 10 == value);
    }

    #[test]
    fn test_get_depth_counts() {
        let reference_df = reference_df();
        let alns_path = format!("{}/{}", TEST_DATA, "pileup/alns_2.csv");

        let expected_depth_file = format!(
            "{}/{}",
            TEST_DATA, "pileup/expected_unique_depth_counts.csv"
        );

        create_dir_all(format!("{}/{}", TEST_DIR, "pileup"))
            .expect("Failed to create output directory");
        let output_path = format!("{}/{}", TEST_DIR, "pileup/unique_depth_counts.csv");

        let df = get_depth_counts(&reference_df, &alns_path, Some(vec![Signals::Unique])).unwrap();
        save_csv(&df, &output_path).unwrap();

        assert!(compare_files(&expected_depth_file, &output_path));
    }

    #[test]
    fn test_get_depth_counts_round1() {
        let reference_df = reference_df();
        let alns_path = format!("{}/{}", TEST_DATA, "pileup/alns_2.csv");

        let expected_depth_file = format!(
            "{}/{}",
            TEST_DATA, "pileup/expected_combined_depth_counts.csv"
        );

        create_dir_all(format!("{}/{}", TEST_DIR, "pileup"))
            .expect("Failed to create output directory");
        let output_path = format!("{}/{}", TEST_DIR, "pileup/combined_depth_counts.csv");

        let df = get_depth_counts_round1(&reference_df, &alns_path).unwrap();
        save_csv(&df, &output_path).unwrap();

        assert!(compare_files(&expected_depth_file, &output_path));
    }
}
