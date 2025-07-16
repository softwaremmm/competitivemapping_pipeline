//! This module reads alignments from a BAM file, filters them and then writes to output BAM files.
//!
//! It has some complexity due to the multithreaded nature.
//! The noodles is multithreaded already.
//! Rayon threadpool is used to manage the main filtering step
//! Other threads which mainly wait on crossbeam channels are spawned directly
//! as they are not CPU bound.
use super::Result;
use noodles::bam;
use noodles::bgzf::MultithreadedReader;
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::BufWriter;
use std::num::NonZero;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::SystemTime;

use crate::alignments::{record_to_alignment_info, Alignment};

use crate::parameters::FilterParams;
use crossbeam::channel::{bounded, Receiver, Sender};
use rayon::ThreadPoolBuilder;

pub mod filter_counts;
use crate::get_target_id_to_ref_and_ani_group;
use filter_counts::{
    FilterCounts, FilterResult, FilterRoundStats, InputStats, ReadStats, ReadType,
};

type SignalFiles = Vec<(String, String)>;

/// Take all the alignments for a single read (or read pair) and filter them based on the provided parameters.
///
/// This function will return a tuple containing:
/// - A vector of filtered BAM records with alignment info.
/// - A `FilterCounts` object containing the counts of filtered alignments.
/// - A signal string indicating the filtering result (e.g., "all_fail", "unique", "winner", "shared", "best").
/// - A `ReadType` indicating the type of read (Mapped, Unmapped, HalfMapped).
fn filter_read_alns(
    aln_group: Vec<bam::Record>,
    params: &FilterParams,
    tid_to_ref_id_and_ani_group: &HashMap<i32, (i32, i32)>,
    ani_group_ordering: &Option<HashMap<i32, u32>>, // Where a lower score is better
    debug: bool,
) -> (Vec<Alignment>, ReadStats) {
    // For paired reads, each part of the pair will have the filter stats calculated separately.
    let mut max_query_coverage = [0.0, 0.0]; // (first in pair, second in pair)
    let mut min_divergence = [1.0, 1.0];
    let mut primary_score = [0, 0];
    let mut expected_divergence = [0.0, 0.0];

    let mut some_unmapped = false;
    let mut alns: Vec<Alignment> = aln_group
        .into_iter()
        .map(|record| {
            let (mut aln, _aln_extra) = record_to_alignment_info(&record, false);
            if aln.is_unmapped {
                some_unmapped = true;
            } else {
                let (ref_id, ani_group) = tid_to_ref_id_and_ani_group
                    .get(&aln.target_id)
                    .unwrap_or_else(|| {
                        panic!("ANI group not found for target_id {}", aln.target_id)
                    });
                aln.ref_id = Some(*ref_id);
                aln.ani_group = Some(*ani_group);
            }
            aln
        })
        .filter(|aln| !aln.is_unmapped)
        .collect();
    let read_type = match (some_unmapped, !alns.is_empty()) {
        (false, _) => ReadType::Mapped,
        (true, false) => ReadType::Unmapped,
        (true, true) => ReadType::HalfMapped,
    };

    for aln in alns.iter() {
        let query_coverage = (aln.query_covered as f32) / (aln.query_length as f32);

        let pair_index = aln.pair_index();
        max_query_coverage[pair_index] = f32::max(max_query_coverage[pair_index], query_coverage);
        min_divergence[pair_index] = f32::min(
            min_divergence[pair_index],
            aln.gap_compressed_seq_divergence,
        );

        if aln.is_primary {
            primary_score[pair_index] = aln.alignment_score;
            expected_divergence[pair_index] = aln.expected_error_rate;
        }
    }

    // Query coverage filters
    let min_query_coverage = params.min_query_coverage.unwrap_or(0.0);
    let min_query_coverage_pc_of_max =
        max_query_coverage.map(|x| x * params.min_query_coverage_pc_of_max.unwrap_or(0.0));

    // Divergence filters
    let max_divergence = params.max_divergence.unwrap_or(1.0);
    let max_divergence_from_expected =
        expected_divergence.map(|x| x + params.max_divergence_from_expected.unwrap_or(1.0));
    let max_divergence_from_best =
        min_divergence.map(|x| x + params.max_divergence_from_best.unwrap_or(1.0));

    // Score filters
    let min_score_fraction = params.min_score_fraction.unwrap_or(0.0);

    // Also need to consider sharing
    let mut filter_counts: FilterCounts = FilterCounts::default();
    let share_threshold = params.share_threshold.unwrap_or(1.0);
    let mut ani_groups: HashSet<i32> = HashSet::new();
    let mut strong_ani_groups: HashSet<i32> = HashSet::new();

    for aln in alns.iter_mut() {
        let query_coverage = (aln.query_covered as f32) / (aln.query_length as f32);
        let pair_index = aln.pair_index();

        if query_coverage < min_query_coverage {
            filter_counts.low_query_cov += 1;
            aln.filter_result = FilterResult::LowQueryCoverage;
            continue;
        }
        if query_coverage < min_query_coverage_pc_of_max[pair_index] {
            filter_counts.low_query_cov_pc_of_max += 1;
            aln.filter_result = FilterResult::LowQueryCoveragePcOfMax;
            continue;
        }
        if aln.gap_compressed_seq_divergence > max_divergence {
            filter_counts.high_divergence += 1;
            aln.filter_result = FilterResult::HighDivergence;
            continue;
        }
        if aln.gap_compressed_seq_divergence > max_divergence_from_expected[pair_index] {
            filter_counts.high_divergence_from_expected += 1;
            aln.filter_result = FilterResult::HighDivergenceFromExpected;
            continue;
        }
        if aln.gap_compressed_seq_divergence > max_divergence_from_best[pair_index] {
            filter_counts.high_divergence_from_best += 1;
            aln.filter_result = FilterResult::HighDivergenceFromBest;
            continue;
        }

        let score_fraction = aln.alignment_score as f32 / primary_score[pair_index] as f32;
        if score_fraction < min_score_fraction {
            filter_counts.low_score += 1;
            aln.filter_result = FilterResult::LowScore;
            continue;
        }

        ani_groups.insert(aln.ani_group.expect("ANI group should be set"));
        if score_fraction < share_threshold {
            filter_counts.weak += 1;
            aln.filter_result = FilterResult::WeakScore;
            continue;
        }

        strong_ani_groups.insert(aln.ani_group.expect("ANI group should be set"));
        aln.filter_result = FilterResult::Passed;
    }

    if !debug {
        alns.retain(|aln| aln.filter_result == FilterResult::Passed);
    }

    if !params.allow_reads_multiple_alns_per_ref {
        fn calculate_hash<T: Hash>(t: &T) -> u64 {
            let mut s = DefaultHasher::new();
            t.hash(&mut s);
            s.finish()
        }
        // Keep the best by alignment score
        alns.sort_by_key(|aln| {
            (
                aln.ref_id.expect("Ref ID should be set"),
                -aln.alignment_score,
                calculate_hash(&aln.ref_start),
                calculate_hash(&aln.target_id),
            )
        });
        // sorting by target_id and ref_start to ensure stable sorting
        // but using hash to avoid bias

        let mut seen_refs = [HashSet::new(), HashSet::new()];

        for aln in alns.iter_mut() {
            if aln.filter_result != FilterResult::Passed {
                continue;
            }
            let ref_id = aln.ref_id.expect("Ref ID should be set");
            let seen_refs = &mut seen_refs[aln.pair_index()];
            if seen_refs.contains(&ref_id) {
                filter_counts.read_has_better_aln_for_ref += 1;
                aln.filter_result = FilterResult::ReadHasBetterAlnForRef;
                continue;
            }
            seen_refs.insert(ref_id);
        }

        if !debug {
            alns.retain(|aln| aln.filter_result == FilterResult::Passed);
        }
    }

    if !strong_ani_groups.is_empty() {
        // If ordering provided then need to select best ANI group
        if let Some(ani_group_ordering) = ani_group_ordering {
            strong_ani_groups.retain(|ani_group| ani_group_ordering.contains_key(ani_group));

            let best_group = strong_ani_groups
                .clone()
                .into_iter()
                .min_by_key(|ani_group| {
                    ani_group_ordering
                        .get(ani_group)
                        .expect("Just checked that this ani_group is in the ordering")
                });

            for aln in alns.iter_mut() {
                if aln.filter_result != FilterResult::Passed {
                    continue;
                }
                let ani_group = aln.ani_group.expect("ANI group should be set");
                if Some(ani_group) == best_group {
                    continue;
                }

                if strong_ani_groups.contains(&ani_group) {
                    filter_counts.ref_lost_draw += 1;
                    aln.filter_result = FilterResult::RefLostDraw;
                    continue;
                }
                filter_counts.target_excluded += 1;
                aln.filter_result = FilterResult::TargetExcluded;
            }

            if !debug {
                alns.retain(|aln| aln.filter_result == FilterResult::Passed);
            }
        }
    }

    filter_counts.passed = alns.iter().filter(|aln| aln.filter_result == FilterResult::Passed).count();

    let signal = match (
        ani_groups.len(),
        strong_ani_groups.len(),
        ani_group_ordering.is_some(),
    ) {
        (_, 0, _) => "all_fail".to_string(),
        (_, _, true) => "best".to_string(),
        (1, _, _) => "unique".to_string(),
        (_, 1, _) => "winner".to_string(),
        (_, _, _) => "shared".to_string(),
    };

    assert!(signal == "all_fail" || filter_counts.passed > 0);

    return (
        alns,
        ReadStats {
            read_type,
            signal,
            filter_counts,
        },
    );
}

#[allow(clippy::too_many_arguments)]
fn process_read_aln_groups(
    aln_group_rx: Receiver<Vec<bam::Record>>,
    thread_pool: &rayon::ThreadPool,
    max_jobs: usize,
    params: FilterParams,
    tid_to_ref_id_and_ani_group: Arc<HashMap<i32, (i32, i32)>>,
    tie_break_order: Arc<Option<HashMap<i32, u32>>>,
    stats_tx: Sender<ReadStats>,
    write_tx: Sender<(Alignment, String)>,
    debug_tx: Option<Sender<Vec<Alignment>>>,
) {
    let job_counter = Arc::new(AtomicUsize::new(0));
    let debug = debug_tx.is_some();
    for group in aln_group_rx {
        while job_counter.load(Ordering::Relaxed) >= max_jobs {
            thread::sleep(std::time::Duration::from_millis(10));
        }
        job_counter.fetch_add(1, Ordering::Relaxed);

        let inner_counter = job_counter.clone();
        let stats_tx = stats_tx.clone();
        let write_tx = write_tx.clone();
        let debug_tx = debug_tx.clone();
        let tid_to_ref_id_and_ani_group = tid_to_ref_id_and_ani_group.clone();
        let tie_break_order = tie_break_order.clone();
        thread_pool.spawn(move || {
            let (mut alns, read_stats) = filter_read_alns(
                group,
                &params,
                &tid_to_ref_id_and_ani_group,
                &tie_break_order,
                debug,
            );

            let signal = read_stats.signal.clone();
            stats_tx.send(read_stats).expect("Stats thread crashed");

            if debug {
                if let Some(debug_tx) = &debug_tx {
                    debug_tx
                        .send(alns.clone())
                        .expect("Debug thread crashed");
                }
                alns.retain(|aln| aln.filter_result == FilterResult::Passed);
            }

            for aln in alns {
                write_tx
                    .send((aln, signal.clone()))
                    .expect("Processor thread crashed");
            }
            inner_counter.fetch_sub(1, Ordering::Relaxed);
        });
    }
}

fn write_alns_to_csv(file_paths: &[(String, String)], alns_rx: Receiver<(Alignment, String)>) {
    let mut writers: HashMap<String, csv::Writer<BufWriter<File>>> = HashMap::new();
    for (signal, path) in file_paths {
        let mut writer = csv::Writer::from_writer(BufWriter::new(File::create(path).unwrap()));

        writer
            .write_record([
                "query_name",
                "is_second_in_pair",
                "target_id",
                "query_length",
                "ref_start",
                "ref_end",
            ])
            .expect("Failed to write header to CSV");
        writers.insert(signal.clone(), writer);
    }

    for (aln, signal) in alns_rx {
        if let Some(writer) = writers.get_mut(&signal) {
            writer
                .write_record(&[
                    aln.read_id,
                    aln.is_second_in_pair.to_string(),
                    aln.target_id.to_string(),
                    aln.query_length.to_string(),
                    aln.ref_start.to_string(),
                    (aln.ref_start + aln.ref_covered - 1).to_string(),
                ])
                .expect("Failed to write alignment to CSV");
        } else {
            panic!("No writer found for signal {signal}");
        }
    }

    // Flush all writers
    for writer in writers.values_mut() {
        writer.flush().expect("Failed to flush CSV writer");
    }
}

fn write_debug_to_csv(file_path: &str, alns_rx: Receiver<Vec<Alignment>>, active: bool) {
    if ! active {
        return;
    }
    let mut writer = csv::Writer::from_writer(BufWriter::new(File::create(file_path).unwrap()));

    for aln_group in alns_rx {
        for aln in aln_group {
            writer.serialize(aln).expect("Failed to serialize alignment");
        }
    }

    writer.flush().expect("Failed to flush debug CSV writer");
}

fn process_stats(
    stats_rx: Receiver<ReadStats>,
    shared_stats_object: Arc<Mutex<(InputStats, FilterRoundStats)>>,
) {
    let mut input_stats = InputStats::default();
    let mut round_stats = FilterRoundStats::default();
    let mut signal_counts = HashMap::new();
    for read_stats in stats_rx {
        round_stats.add_filter_counts(&read_stats.filter_counts);
        if read_stats.filter_counts.passed > 0 {
            round_stats.passed_reads += 1;

            // Count the signal for this read
            let signal_count = signal_counts.entry(read_stats.signal.clone()).or_insert(0);
            *signal_count += read_stats.filter_counts.passed;
        }
        match read_stats.read_type {
            ReadType::Mapped => input_stats.mapped_reads += 1,
            ReadType::Unmapped => input_stats.unmapped_reads += 1,
            ReadType::HalfMapped => input_stats.half_mapped_reads += 1,
        }
    }

    // Now update the shared stats object
    let mut shared_stats = shared_stats_object.lock().unwrap();
    shared_stats.0.mapped_reads = input_stats.mapped_reads;
    shared_stats.0.unmapped_reads = input_stats.unmapped_reads;
    shared_stats.0.half_mapped_reads = input_stats.half_mapped_reads;
    shared_stats.1.passed_reads = round_stats.passed_reads;
    shared_stats.1.filter_counts = round_stats.filter_counts;
    shared_stats.1.signal_counts = signal_counts.into_iter().collect();
}

fn get_ani_group_order(tie_break_order: Option<DataFrame>) -> Result<Option<HashMap<i32, u32>>> {
    if let Some(tie_break_order) = tie_break_order {
        let df = tie_break_order
            .lazy()
            .select([col("ani_group")])
            .unique_stable(None, Default::default())
            .with_row_index("ranking", Some(0))
            .collect()?;
        let ani_groups = df
            .column("ani_group")?
            .i64()?
            .into_no_null_iter()
            .map(|ani_group| ani_group as i32)
            .collect::<Vec<i32>>();
        let rankings = df
            .column("ranking")?
            .u32()?
            .into_no_null_iter()
            .collect::<Vec<u32>>();
        let mapping = ani_groups
            .into_iter()
            .zip(rankings)
            .collect::<HashMap<i32, u32>>();
        return Ok(Some(mapping));
    }

    return Ok(None);
}

/// Determine the number of threads to use for reading and processing
fn determine_threads(requested_threads: Option<usize>) -> (usize, usize) {
    let available_threads = match requested_threads {
        Some(threads) => threads,
        None => thread::available_parallelism()
            .unwrap_or(std::num::NonZeroUsize::MIN)
            .get(),
    };

    if available_threads <= 2 {
        return (1, 1);
    }

    let process_threads = available_threads / 2;
    let read_threads = if available_threads % 2 == 0 {
        process_threads
    } else {
        process_threads + 1
    };
    return (read_threads, process_threads);
}

/// Function filters alignments in bam file based on the provided parameters.
pub fn filter_bam(
    bam_path: &str,
    ref_df: &DataFrame,
    params: FilterParams,
    output_root: &str,
    threads: Option<usize>,
    tie_break_order: Option<DataFrame>,
    debug: bool,
) -> Result<(SignalFiles, InputStats, FilterRoundStats)> {
    let now = SystemTime::now();

    let is_round_2 = tie_break_order.is_some();

    let (read_threads, process_threads) = determine_threads(threads);
    if debug {
        println!(
            "Using {read_threads} read threads, {process_threads} process threads"
        );
    }

    // Reading is using a multithreaded reader
    let decoder = MultithreadedReader::with_worker_count(
        NonZero::new(read_threads).unwrap(),
        File::open(bam_path)?,
    );
    let mut reader = bam::io::Reader::from(decoder);
    let _header = reader.read_header()?; // Need to read header otherwise error when reading records

    // Set up rayon pool to be used for the main filtering step
    let pool = ThreadPoolBuilder::new()
        .num_threads(process_threads)
        .build()
        .unwrap();

    // Set up channels for communication
    let max_jobs = 500;
    let (aln_group_tx, aln_group_rx) = bounded::<Vec<bam::Record>>(max_jobs); // Groups
    let (stats_tx, stats_rx) = bounded::<ReadStats>(max_jobs); // Results
    let (write_tx, write_rx) = bounded::<(Alignment, String)>(max_jobs); // Results for CSV
    let (debug_tx, debug_rx) = bounded::<Vec<Alignment>>(max_jobs);

    // Set up shared stats objects
    let shared_stats_object = Arc::new(Mutex::new((
        InputStats::default(),
        FilterRoundStats::default(),
    )));
    let tid_to_ref_id_and_ani_group = Arc::new(get_target_id_to_ref_and_ani_group(ref_df)?);
    let tie_break_order = Arc::new(get_ani_group_order(tie_break_order)?);

    let signals = if is_round_2 {
        vec!["best"]
    } else {
        vec!["unique", "winner", "shared"]
    };
    let signal_paths: SignalFiles = signals
        .iter()
        .map(|signal| {
            let path = format!("{output_root}{signal}_alns.csv");
            (signal.to_string(), path)
        })
        .collect();

    let process_handle = thread::spawn({
        let stats_tx = stats_tx.clone();
        let write_tx = write_tx.clone();
        let debug_tx = if debug {
            Some(debug_tx.clone())
        } else { None };
        move || {
            process_read_aln_groups(
                aln_group_rx,
                &pool,
                max_jobs,
                params,
                tid_to_ref_id_and_ani_group,
                tie_break_order,
                stats_tx,
                write_tx,
                debug_tx,
            );
        }
    });

    // Writing here is just writing to csv so only need one thread
    let write_handle = thread::spawn({
        let file_paths = signal_paths.clone();
        let alns_rx = write_rx.clone();
        move || {
            write_alns_to_csv(&file_paths, alns_rx);
        }
    });

    let debug_handle = thread::spawn({
        let round_index = if is_round_2 {"round_2"} else {"round_1"};
        let path = format!("{output_root}debug_alns_{round_index}.csv");
        let debug_rx = debug_rx.clone();
        move || {
            write_debug_to_csv(&path, debug_rx, debug);
        }
    });

    // Tracking stats is not CPU bound, so we can use a separate thread
    let stats_handle = thread::spawn({
        let shared_stats_object = shared_stats_object.clone();
        move || {
            process_stats(stats_rx, shared_stats_object);
        }
    });

    // Main thread reads bam alignments and groups by read
    let mut buffer = Vec::new();
    let mut current_query: Option<String> = None;
    let mut top_level_counts = InputStats::default();
    let update_frequency = if debug { 100000 } else { 1000000 };
    for (j, record) in reader.records().enumerate() {
        if j % update_frequency == 0 && j > 0 {
            println!(
                "Read {} records. Time elapsed {}",
                j,
                now.elapsed().unwrap().as_secs()
            );
            if debug {
                println!(
                    "process: {}, writing: {}, stats: {}",
                    aln_group_tx.len(),
                    write_tx.len(),
                    stats_tx.len()
                );
            }
        }
        let record = record.unwrap();
        top_level_counts.total_alns += 1;

        let read_id = record
            .name()
            .map(|name| name.to_string())
            .unwrap_or("".to_string());

        match &current_query {
            Some(cur) if *cur == read_id => buffer.push(record),
            Some(_) => {
                aln_group_tx.send(buffer).expect("Worker thread crashed");
                top_level_counts.total_reads += 1;
                buffer = vec![record];
                current_query = Some(read_id);
            }
            None => {
                current_query = Some(read_id);
                buffer.push(record);
            }
        }
    }
    if !buffer.is_empty() {
        aln_group_tx.send(buffer).expect("Worker thread crashed");
        top_level_counts.total_reads += 1;
    }

    drop(aln_group_tx);
    drop(write_tx);
    drop(stats_tx);
    drop(debug_tx);

    process_handle
        .join()
        .expect("Failed to join process handle");
    write_handle.join().expect("Failed to join write handle");
    stats_handle.join().expect("Failed to join stats handle");
    debug_handle.join().expect("Failed to join debug handle");

    println!("Reading bam took {:?}", now.elapsed().unwrap());

    // Copy the final stats
    let shared_stats = shared_stats_object.lock().unwrap();
    let mut input_stats = shared_stats.0.clone();
    input_stats.total_reads = top_level_counts.total_reads;
    input_stats.total_alns = top_level_counts.total_alns;
    let filter_stats = shared_stats.1.clone();

    Ok((signal_paths, input_stats, filter_stats))
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::shared_test_functions::{bam_to_records, create_bam_from_lines, TEST_DIR};

    #[test]
    fn test_filter_reads_one_per_ref() -> Result<()> {
        let params = FilterParams {
            allow_reads_multiple_alns_per_ref: false,
            share_threshold: Some(0.9),
            ..Default::default()
        };
        let tid_to_ref_id_and_ani_group = HashMap::from([(0, (0, 0)), (1, (1, 1)), (2, (2, 2))]);

        let sam_lines = [
            "@HD\tVN:1.6\tSO:queryname",
            "@SQ\tSN:ref1\tLN:1000",
            "@SQ\tSN:ref2\tLN:1000",
            "1\t0\tref1\t1\t60\t100M\t*\t0\t0\t*\t*\tNM:i:21\tms:i:2244\tAS:i:100\tde:f:0.0",
            "1\t256\tref1\t1\t60\t100M\t*\t0\t0\t*\t*\tNM:i:21\tms:i:2244\tAS:i:99\tde:f:0.00",
        ]
        .iter()
        .map(|line| line.to_string())
        .collect::<Vec<String>>();

        let file_path = format!("{TEST_DIR}/test_filter_reads/one_per_ref1.bam");
        create_bam_from_lines(&file_path, sam_lines)?;
        let records = bam_to_records(&file_path)?;

        let (filtered_alns, read_stats) =
            filter_read_alns(records, &params, &tid_to_ref_id_and_ani_group, &None, false);

        assert_eq!(read_stats.read_type, ReadType::Mapped);
        assert_eq!(read_stats.signal, "unique");
        assert_eq!(read_stats.filter_counts.passed, 1);
        assert_eq!(filtered_alns.len(), 1);

        // Repeat but with two references
        let sam_lines = [
            "@HD\tVN:1.6\tSO:queryname",
            "@SQ\tSN:ref1\tLN:1000",
            "@SQ\tSN:ref2\tLN:1000",
            "1\t0\tref1\t1\t60\t100M\t*\t0\t0\t*\t*\tNM:i:21\tms:i:2244\tAS:i:100\tde:f:0.0",
            "1\t256\tref1\t1\t60\t100M\t*\t0\t0\t*\t*\tNM:i:21\tms:i:2244\tAS:i:99\tde:f:0.00",
            "1\t256\tref2\t1\t60\t100M\t*\t0\t0\t*\t*\tNM:i:21\tms:i:2244\tAS:i:99\tde:f:0.00",
        ]
        .iter()
        .map(|line| line.to_string())
        .collect::<Vec<String>>();

        let file_path = format!("{TEST_DIR}/test_filter_reads/one_per_ref2.bam");
        create_bam_from_lines(&file_path, sam_lines)?;
        let records = bam_to_records(&file_path)?;

        let (filtered_alns, read_stats) =
            filter_read_alns(records, &params, &tid_to_ref_id_and_ani_group, &None, false);

        assert_eq!(read_stats.read_type, ReadType::Mapped);
        assert_eq!(read_stats.signal, "shared");
        assert_eq!(read_stats.filter_counts.passed, 2);
        assert_eq!(filtered_alns.len(), 2);

        Ok(())
    }
}
