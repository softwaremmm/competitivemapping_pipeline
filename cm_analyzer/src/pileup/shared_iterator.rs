//! Module for counting read depths from alignment files
//!
//! Currently set to use alignments in CSV format
//! But can be adapted to use BAM files fairly easily
//!
//! Used to be multithreaded with BAM files,
//! but that is much slower with csvs due to repeated reading and filtering csv file
//! And often majority of alignments are to the top target_id anyway

use crate::Result;
use itertools::merge_join_by;
use itertools::EitherOrBoth::{Both, Left, Right};
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
    unique_alns_path: &str,
    winner_alns_path: &str,
) -> Result<DataFrame> {
    let now = SystemTime::now();

    let unique_depth_iter = CsvPileupIterator::from_path(unique_alns_path, None)
        .expect("Failed to create iterator from path");
    let winner_depth_iter = CsvPileupIterator::from_path(winner_alns_path, None)
        .expect("Failed to create iterator from path");

    let merged_iter = merge_join_by(
        unique_depth_iter,
        winner_depth_iter,
        |a, b| (a.0, a.1).cmp(&(b.0, b.1)), // Compare by target_id and pos
    );

    let mut unique_depth_counts: HashMap<(u32, u32), u32> = HashMap::new();
    let mut winner_depth_counts: HashMap<(u32, u32), u32> = HashMap::new();

    merged_iter
        .map(|item| match item {
            Both((tid, _, unique_depth), (_, _, winner_depth)) => {
                (tid, Some(unique_depth), Some(winner_depth + unique_depth))
            }
            Left((tid, _, unique_depth)) => (tid, Some(unique_depth), Some(unique_depth)),
            Right((tid, _, winner_depth)) => (tid, None, Some(winner_depth)),
        })
        .for_each(|(tid, unique_depth, winner_depth)| {
            if let Some(unique_depth) = unique_depth {
                let count = unique_depth_counts.entry((tid, unique_depth)).or_insert(0);
                *count += 1;
            }
            if let Some(winner_depth) = winner_depth {
                let count = winner_depth_counts.entry((tid, winner_depth)).or_insert(0);
                *count += 1;
            }
        });

    let unique_df = depth_count_hashmap_to_dataframe(unique_depth_counts)?;
    let winner_df = depth_count_hashmap_to_dataframe(winner_depth_counts)?;
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
    .collect()?;

    let final_df = combined
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
        .group_by(["ref_id", "depth_type", "depth"])
        .agg([col("count").sum().alias("count")])
        .sort(
            ["ref_id", "depth", "depth_type"],
            SortMultipleOptions::default(),
        )
        .collect()?;

    println!("Reading took {:?} overall", now.elapsed().unwrap());
    Ok(final_df)
}

/// Currently set up to use unique and winner
fn read_depth_count_for_region<I>(
    winner_iter: I,
    unique_iter: I,
    target_id: u32, // used as a check
) -> Result<DataFrame>
where
    I: Iterator<Item = (u32, u32, u32)>,
{
    let mut unique_depth_counts: HashMap<u32, u32> = HashMap::new();
    let mut winner_depth_counts: HashMap<u32, u32> = HashMap::new();

    let winner_iter = winner_iter
        .filter(|(tid, _, _)| *tid == target_id)
        .map(|(_, pos, depth)| (pos, depth));
    let unique_iter = unique_iter
        .filter(|(tid, _, _)| *tid == target_id)
        .map(|(_, pos, depth)| (pos, depth));

    let merged = merge_join_by(unique_iter, winner_iter, |a, b| a.0.cmp(&b.0));

    merged
        .map(|item| match item {
            Both((_, unique_depth), (_, winner_depth)) => {
                (Some(unique_depth), Some(winner_depth + unique_depth))
            }
            Left((_, unique_depth)) => (Some(unique_depth), Some(unique_depth)),
            Right((_, winner_depth)) => (None, Some(winner_depth)),
        })
        .for_each(|(unique_depth, winner_depth)| {
            if let Some(unique_depth) = unique_depth {
                let count = unique_depth_counts.entry(unique_depth).or_insert(0);
                *count += 1;
            }
            if let Some(winner_depth) = winner_depth {
                let count = winner_depth_counts.entry(winner_depth).or_insert(0);
                *count += 1;
            }
        });

    let unique_df = hashmap_to_dataframe!(unique_depth_counts, "depth".into(), "count".into())?;
    let winner_df = hashmap_to_dataframe!(winner_depth_counts, "depth".into(), "count".into())?;
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
    .with_column(lit(target_id).alias("target_id"))
    .filter(col("depth").gt(0))
    .collect()?;

    Ok(combined)
}

/// Currently set up to use unique and winner
/// Not in use yet
/// This function is a more generic version of the above
#[allow(dead_code)]
fn read_depth_count_for_region_custom<I>(
    depth_iterators: Vec<I>,
    target_id: u32,
    combining_logic: fn(Vec<u32>) -> Vec<u32>,
    combined_labels: Vec<String>,
) -> Result<DataFrame>
where
    I: Iterator<Item = (u32, u32)>,
{
    let mut depth_counts: Vec<HashMap<u32, u32>> = Vec::new();

    let iterator = merge_pos_iterators(depth_iterators);

    for (_pos, depths) in iterator {
        let combined_depths = combining_logic(depths);
        for (i, depth) in combined_depths.iter().enumerate() {
            let count = depth_counts[i].entry(depth.clone()).or_insert(0);
            *count += 1;
        }
    }

    // check that combined_labels and depth_counts have the same length
    if depth_counts.len() != combined_labels.len() {
        return Err("Length of combined_labels and depth_counts do not match".into());
    }

    let dataframes = depth_counts
        .into_iter()
        .zip(combined_labels)
        .map(|(depth_count, label)| {
            let df = hashmap_to_dataframe!(depth_count, "depth".into(), "count".into()).unwrap();
            df.lazy().with_column(lit(label).alias("depth_type"))
        })
        .collect::<Vec<LazyFrame>>();

    let combined = concat(dataframes, UnionArgs::default())?
        .with_column(lit(target_id).alias("target_id"))
        .filter(col("depth").gt(0))
        .collect()?;

    Ok(combined)
}

// Generic function to merge multiple iterators over positions and depths
// Unlikely to be as efficient as just mergining two with merge_join_by
// but this is a more generic solution
fn merge_pos_iterators<I>(iterators: Vec<I>) -> impl Iterator<Item = (u32, Vec<u32>)>
where
    I: Iterator<Item = (u32, u32)>,
{
    let mut peekable_iters: Vec<_> = iterators.into_iter().map(|iter| iter.peekable()).collect();

    std::iter::from_fn(move || {
        let mut min_pos: Option<u32> = None;

        for peekable_iter in &mut peekable_iters {
            if let Some((pos, _)) = peekable_iter.peek() {
                if min_pos.is_none() || *pos < min_pos.unwrap() {
                    min_pos = Some(*pos);
                }
            }
        }

        // If all iterators are exhausted, return None
        if min_pos.is_none() {
            return None;
        }
        let min_pos = min_pos.unwrap();

        let mut depths = vec![];
        for peekable_iter in &mut peekable_iters {
            if let Some((pos, depth)) = peekable_iter.peek() {
                if *pos == min_pos {
                    depths.push(*depth);
                    peekable_iter.next();
                    continue;
                }
            }
            depths.push(0);
        }
        Some((min_pos, depths))
    })
}

pub fn get_depth_counts_single(reference_df: &DataFrame, alns_paths: &str) -> Result<DataFrame> {
    let now = SystemTime::now();

    let depth_iter = CsvPileupIterator::from_path(alns_paths, None)
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
    println!("Depth: {:?}", final_df);
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

        println!("{:?}", df);

        assert_eq!(df.shape(), (3, 2));
        let key: u32 = df.column("key").unwrap().u32().unwrap().get(0).unwrap();
        let value: u32 = df.column("value").unwrap().u32().unwrap().get(0).unwrap();
        assert!(key * 10 == value);
    }

    #[test]
    fn test_get_depth_counts_simple() {
        let reference_df = reference_df();
        let unique_alns_path = format!("{}/{}", TEST_DATA, "pileup/unique_alns.csv");
        let expected_depth_file = format!(
            "{}/{}",
            TEST_DATA, "pileup/expected_unique_depth_counts.csv"
        );

        create_dir_all(format!("{}/{}", TEST_DIR, "pileup"))
            .expect("Failed to create output directory");
        let output_path = format!("{}/{}", TEST_DIR, "pileup/unique_depth_counts.csv");

        let df = get_depth_counts_single(&reference_df, &unique_alns_path).unwrap();
        save_csv(&df, &output_path).unwrap();

        assert!(compare_files(&expected_depth_file, &output_path));
    }

    #[test]
    fn test_get_depth_counts() {
        let reference_df = reference_df();
        let unique_alns_path = format!("{}/{}", TEST_DATA, "pileup/unique_alns.csv");
        let winner_alns_path = format!("{}/{}", TEST_DATA, "pileup/winner_alns.csv");
        let expected_depth_file = format!(
            "{}/{}",
            TEST_DATA, "pileup/expected_combined_depth_counts.csv"
        );

        create_dir_all(format!("{}/{}", TEST_DIR, "pileup"))
            .expect("Failed to create output directory");
        let output_path = format!("{}/{}", TEST_DIR, "pileup/combined_depth_counts.csv");

        let df =
            get_depth_counts_round1(&reference_df, &unique_alns_path, &winner_alns_path).unwrap();
        save_csv(&df, &output_path).unwrap();

        assert!(compare_files(&expected_depth_file, &output_path));
    }
}
