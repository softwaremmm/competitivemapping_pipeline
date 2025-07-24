use crate::depth_counts_analysis::{count_alns, order_best_refs, summarise_depth};
use crate::filter_reads::filter_bam;
use crate::filter_reads::filter_counts::{OverallStats, Signals};
use crate::parameters::Params;
use crate::pileup::shared_iterator::{get_depth_counts, get_depth_counts_round1};
use crate::{make_reference_df, save_csv};
use clap::Parser;
use polars::prelude::*;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::fs::File;
use std::time::SystemTime;

/// Analyzes competitive mapping data
#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
pub struct AnalyzeArgs {
    /// Sets the input alignment file
    #[clap(short, long, value_name = "BAM")]
    pub input_bam: String,

    /// Sets the input contigs CSV file
    #[clap(short, long, value_name = "CSV")]
    pub contigs: String,

    #[clap(short, long, value_name = "OUTPUT_ROOT")]
    pub output_root: String,

    #[clap(short, long, value_name = "PARAMETERS_YAML")]
    pub parameters: String,

    #[clap(short, long)]
    pub threads: Option<usize>,

    #[clap(short, long)]
    pub debug: bool,
}

pub fn analyze_alignments(args: AnalyzeArgs) -> Result<(), Box<dyn std::error::Error>> {
    let program_start = SystemTime::now();
    let mut stats = OverallStats::default();

    let params: Params = serde_yaml::from_reader(
        File::open(args.parameters)
            .map_err(|e| format!("Failed to read params file. Error: {e}"))?,
    )?;

    let reference_df = make_reference_df(&args.input_bam, &args.contigs)?;
    save_csv(
        &reference_df,
        format!("{}{}", args.output_root, "references.csv"),
    )?;
    if args.debug {
        println!("Reference DataFrame: {reference_df:?}");
    }

    let (alns_path, input_stats, round1_stats) = filter_bam(
        &args.input_bam,
        &reference_df,
        params.round_1,
        &args.output_root,
        args.threads,
        None,
        args.debug,
    )?;
    stats.input_stats = input_stats;
    stats.filter_round1 = round1_stats;

    // Depth counts are ref_id, depth_type, depth, count
    let depth_counts = get_depth_counts_round1(&reference_df, &alns_path, args.threads)?;
    let round1_summarised_depth = summarise_depth(&depth_counts, &reference_df)?;

    let ref_tie_breaker_order =
        order_best_refs(&round1_summarised_depth, params.min_unique_coverage)?;
    if args.debug {
        println!("Ref tie breaker order: {ref_tie_breaker_order:?}");
    }
    save_csv(
        &ref_tie_breaker_order,
        format!("{}{}", args.output_root, "ref_tie_breaker_order.csv"),
    )?;

    let (best_alns_path, _, round2_stats) = filter_bam(
        &args.input_bam,
        &reference_df,
        params.round_2,
        &args.output_root,
        args.threads,
        Some(ref_tie_breaker_order),
        args.debug,
    )?;

    stats.filter_round2 = round2_stats;
    let stats_file = File::create(format!("{}{}", args.output_root, "stats.yaml"))?;
    serde_yaml::to_writer(stats_file, &stats.sorted())?;

    let final_depth_counts = get_depth_counts(&reference_df, &best_alns_path, None)?;
    let round2_summarised_depth = summarise_depth(&final_depth_counts, &reference_df)?;

    save_csv(
        &concat(
            [final_depth_counts.lazy(), depth_counts.lazy()],
            Default::default(),
        )?
        .collect()?,
        format!("{}{}", args.output_root, "depth_counts.csv"),
    )?;

    // Now just need to combine all into one table

    let depth_type_order = df!(
        "depth_type" => ["final", "unique", "winner", "good"],
        "depth_type_order" => [3, 2, 1, 0]
    )?;
    let ref_order = round2_summarised_depth
        .clone()
        .lazy()
        .sort(["mean_depth"], SortMultipleOptions::default())
        .select([col("ref_id")])
        .with_row_index("ref_order", Some(1));


    let pool = ThreadPoolBuilder::new()
        .num_threads(args.threads.unwrap_or(1))
        .build()
        .unwrap();

    let read_count_dfs: Vec<_> = pool.install(|| {
        let jobs = vec![
            (&best_alns_path, None, "final"),
            (&alns_path, Some(vec![Signals::Unique]), "unique"),
            (
                &alns_path,
                Some(vec![Signals::Unique, Signals::Winner]),
                "winner",
            ),
            (
                &alns_path,
                Some(vec![Signals::Unique, Signals::Winner, Signals::Shared]),
                "good",
            ),
        ];

        jobs.into_par_iter()
            .map(|(path, signals, read_type)| {
                let count_df = count_alns(path, &reference_df, signals.clone())
                    .unwrap_or_else(|_| panic!("Failed to count alignments for {signals:?}"));
                count_df.lazy()
                    .with_column(lit(read_type).alias("depth_type"))
            })
            .collect()
    });

    let read_counts = concat(
        read_count_dfs,
        UnionArgs::default(),
    )?;

    let combined_summary = concat(
        [
            round2_summarised_depth.clone().lazy(),
            round1_summarised_depth.clone().lazy(),
        ],
        UnionArgs::default(),
    )?
    .join(
        read_counts,
        [col("ref_id"), col("depth_type")],
        [col("ref_id"), col("depth_type")],
        JoinArgs::new(JoinType::Left),
    )
    .join(
        depth_type_order.lazy(),
        [col("depth_type")],
        [col("depth_type")],
        JoinArgs::new(JoinType::Left),
    )
    .join(
        ref_order,
        [col("ref_id")],
        [col("ref_id")],
        JoinArgs::new(JoinType::Left),
    )
    .sort(
        ["depth_type_order", "ref_order", "mean_depth"],
        SortMultipleOptions::default()
            .with_order_descending(true)
            .with_nulls_last(true),
    )
    .drop(["depth_type_order", "ref_order", "ref_id"])
    .select([cols([
        "depth_type",
        "reference",
        "species",
        "ani_group",
        "ref_length",
        "reads",
        "coverage",
        "mean_depth",
        "median_nonzero_depth",
        "simple_expected_coverage",
        "robust_depth_estimate",
        "robust_expected_coverage",
    ])])
    .collect()?;

    save_csv(
        &combined_summary,
        format!("{}{}", args.output_root, "alignment_summary.csv"),
    )?;

    println!("Total time: {:?}", program_start.elapsed().unwrap());

    return Ok(());
}
