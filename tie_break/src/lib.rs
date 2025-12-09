use noodles::bam;
use noodles::bam::bai;
use polars::prelude::*;
use std::time::SystemTime;
use std::{collections::HashMap, fs::File};

use std::process::Command;

pub mod alignments;
pub mod analyze; // top level module
pub mod depth_counts_analysis;
pub mod filter_reads;
pub mod parameters;
pub mod pileup;

#[cfg(test)]
pub mod shared_test_functions;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[macro_export]
macro_rules! struct_to_dataframe {
    ($input:expr, [$($field:ident),+]) => {
        {
            let len = $input.len().to_owned();

            // Extract the field values into separate vectors
            $(let mut $field = Vec::with_capacity(len);)*

            for e in $input.into_iter() {
                $($field.push(e.$field);)*
            }
            df! {
                $(stringify!($field) => $field,)*
            }
        }
    };
}

pub fn save_csv<S>(df: &DataFrame, path: S) -> Result<()>
where
    S: AsRef<std::path::Path>,
{
    CsvWriter::new(File::create(path)?).finish(&mut df.clone())?;
    Ok(())
}

pub fn save_parquet<S>(df: &mut DataFrame, path: S) -> Result<()>
where
    S: AsRef<std::path::Path>,
{
    ParquetWriter::new(File::create(path)?)
        .with_compression(ParquetCompression::Zstd(None))
        .finish(df)?;
    Ok(())
}

pub fn make_bam_index(bam_path: &str) -> Result<()> {
    let index = bam::fs::index(bam_path).unwrap();
    bai::fs::write(format!("{bam_path}.bai"), &index).unwrap();
    Ok(())
}

/// Sort using samtools as it already has code to keep memory low
pub fn sort_bam(bam_path: &str, output_path: &str, threads: Option<usize>) -> Result<()> {
    let now = SystemTime::now();
    let threads = threads.unwrap_or(1);
    let status = Command::new("samtools")
        .args([
            "sort",
            "-@",
            &threads.to_string(),
            "-o",
            output_path,
            bam_path,
        ])
        .status()?;
    make_bam_index(output_path)?;
    if status.success() {
        println!(
            "BAM {} sorted successfully in {:?}.",
            output_path,
            now.elapsed()?
        );
    } else {
        eprintln!("samtools sort failed with exit code: {status}");
    }
    Ok(())
}

pub fn make_reference_df(bam_path: &str, contigs_csv_path: &str) -> Result<DataFrame> {
    let mut contigs_df: DataFrame = CsvReadOptions::default()
        .with_infer_schema_length(None)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(contigs_csv_path.into()))?
        .finish()?;

    // Add ref_id
    let ref_index = contigs_df
        .clone()
        .lazy()
        .select([col("reference")])
        .unique_stable(None, UniqueKeepStrategy::First)
        .with_row_index("ref_id", Some(1));

    contigs_df = contigs_df
        .lazy()
        .join(
            ref_index,
            [col("reference")],
            [col("reference")],
            JoinArgs::default(),
        )
        .collect()?;

    // check if contigs_df has species column.
    if contigs_df.get_column_index("species").is_none() {
        contigs_df = contigs_df
            .lazy()
            .with_columns([col("reference").alias("species")])
            .collect()?;
    }
    // check if contigs_df has ani_group column.
    if contigs_df.get_column_index("ani_group").is_none() {
        contigs_df = contigs_df
            .lazy()
            .with_columns([(col("ref_id").cast(DataType::Int64) * lit(-1)).alias("ani_group")])
            .collect()?;
    }
    else {
        // need to replace null ani_group with -1 * ref_id
        contigs_df = contigs_df
            .lazy()
            .with_columns([when(col("ani_group").is_null())
                .then(col("ref_id").cast(DataType::Int64) * lit(-1))
                .otherwise(col("ani_group"))
                .alias("ani_group")])
            .collect()?;
    }

    let contigs_df = contigs_df
        .lazy()
        .select([
            cols(["rname", "reference", "ref_id", "species", "ani_group"]),
            col("totallength").alias("ref_length"),
        ])
        .collect()?;

    let header = bam::io::reader::Builder
        .build_from_path(bam_path)?
        .read_header()?;

    let ref_sequences = header
        .reference_sequences()
        .keys()
        .map(|k| k.to_string())
        .collect::<Vec<String>>();

    let ref_df: DataFrame = df! {
        "rname" => ref_sequences
    }?
    .lazy()
    .with_row_index("target_id", Some(0))
    .join(
        contigs_df.lazy(),
        [col("rname")],
        [col("rname")],
        JoinArgs::default(),
    )
    .collect()?;

    println!("References loaded");

    Ok(ref_df)
}

pub fn get_target_id_to_ref_and_ani_group(ref_df: &DataFrame) -> Result<HashMap<i32, (i32, i32)>> {
    let target_ids = ref_df.column("target_id")?.u32()?.into_no_null_iter();
    let ref_ids = ref_df.column("ref_id")?.u32()?.into_no_null_iter();
    let ani_groups = ref_df.column("ani_group")?.i64()?.into_no_null_iter();
    Ok(target_ids
        .zip(ref_ids)
        .zip(ani_groups)
        .map(|((target_id, ref_id), ani_group)| {
            (target_id as i32, (ref_id as i32, ani_group as i32))
        })
        .collect())
}

pub fn take_random_subset(
    df: &DataFrame,
    n: i32,
    target_col: &str,
    ref_df: &DataFrame,
) -> Result<DataFrame> {
    let n = Series::new("".into(), &[n]);
    let subset = df
        .clone()
        .lazy()
        .select([col(target_col).unique_stable()])
        .collect()?
        .sample_n(&n, false, true, None)?;
    let subset_df = df
        .clone()
        .lazy()
        .join(
            subset.lazy(),
            [col(target_col)],
            [col(target_col)],
            JoinArgs::new(JoinType::Right),
        )
        .join(
            ref_df
                .clone()
                .lazy()
                .select([cols(["ref_id", "reference", "species"])])
                .unique(None, UniqueKeepStrategy::First),
            [col("ref_id")],
            [col("ref_id")],
            JoinArgs::default(),
        )
        .collect()?;
    Ok(subset_df)
}
