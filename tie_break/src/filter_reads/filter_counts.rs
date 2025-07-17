use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq)]
pub enum ReadType {
    Mapped,
    Unmapped,
    HalfMapped,
}

// Signal converts to u8 for csv writing/reading
#[repr(u8)]
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq, Copy, Hash)]
#[serde(into = "u8", try_from = "u8")]
pub enum Signals {
    Unique = 0,
    Winner = 1,
    Shared = 2,
    Best = 3,
    Failed = 4,
}

impl From<Signals> for u8 {
    fn from(s: Signals) -> u8 {
        s as u8
    }
}

impl TryFrom<u8> for Signals {
    type Error = String;

    fn try_from(v: u8) -> Result<Self, Self::Error> {
        match v {
            0 => Ok(Signals::Unique),
            1 => Ok(Signals::Winner),
            2 => Ok(Signals::Shared),
            3 => Ok(Signals::Best),
            4 => Ok(Signals::Failed),
            _ => Err(format!("Invalid value for Signals enum: {v}")),
        }
    }
}

// implement to_string for Signals
impl std::fmt::Display for Signals {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Signals::Unique => write!(f, "unique"),
            Signals::Winner => write!(f, "winner"),
            Signals::Shared => write!(f, "shared"),
            Signals::Best => write!(f, "best"),
            Signals::Failed => write!(f, "failed"),
        }
    }
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct OverallStats {
    pub input_stats: InputStats,
    pub filter_round1: FilterRoundStats,
    pub filter_round2: FilterRoundStats,
}

impl OverallStats {
    pub fn sorted(self) -> Self {
        Self {
            input_stats: self.input_stats,
            filter_round1: self.filter_round1.sorted(),
            filter_round2: self.filter_round2.sorted(),
        }
    }
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct InputStats {
    pub total_reads: usize,
    pub total_alns: usize,
    pub mapped_reads: usize,
    pub unmapped_reads: usize,
    pub half_mapped_reads: usize, // only one read in pair is mapped
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct FilterRoundStats {
    pub passed_reads: usize, // reads where at least one alignment passed
    pub filter_counts: FilterCounts,
    pub signal_counts: IndexMap<String, usize>,
}

impl FilterRoundStats {
    /// Returns a clone with `signal_counts` sorted by key
    pub fn sorted(self) -> Self {
        let mut sorted_counts: Vec<_> = self.signal_counts.into_iter().collect();
        sorted_counts.sort_by(|a, b| a.1.cmp(&b.1).reverse());

        let sorted_signal_counts: IndexMap<_, _> = sorted_counts.into_iter().collect();

        Self {
            signal_counts: sorted_signal_counts,
            ..self
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
pub enum FilterResult {
    #[default]
    Passed,
    Unmapped,
    LowQueryCoverage,
    LowQueryCoveragePcOfMax,
    HighDivergence,
    HighDivergenceFromExpected,
    HighDivergenceFromBest,
    LowScore,
    WeakScore,
    ReadHasBetterAlnForRef,
    RefLostDraw,
    TargetExcluded,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct FilterCounts {
    pub unmapped: usize,
    pub low_query_cov: usize,
    pub low_query_cov_pc_of_max: usize,
    pub high_divergence: usize,
    pub high_divergence_from_expected: usize,
    pub high_divergence_from_best: usize,
    pub low_score: usize,
    pub weak: usize,
    pub read_has_better_aln_for_ref: usize,
    pub ref_lost_draw: usize,
    pub target_excluded: usize,
    pub passed: usize,
}

impl FilterRoundStats {
    pub fn add_filter_counts(&mut self, other: &FilterCounts) {
        self.filter_counts.low_query_cov += other.low_query_cov;
        self.filter_counts.low_query_cov_pc_of_max += other.low_query_cov_pc_of_max;
        self.filter_counts.high_divergence += other.high_divergence;
        self.filter_counts.high_divergence_from_expected += other.high_divergence_from_expected;
        self.filter_counts.high_divergence_from_best += other.high_divergence_from_best;
        self.filter_counts.low_score += other.low_score;
        self.filter_counts.weak += other.weak;
        self.filter_counts.read_has_better_aln_for_ref += other.read_has_better_aln_for_ref;
        self.filter_counts.ref_lost_draw += other.ref_lost_draw;
        self.filter_counts.target_excluded += other.target_excluded;
        self.filter_counts.passed += other.passed;
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ReadStats {
    pub read_type: ReadType,
    pub signal: Signals,
    pub filter_counts: FilterCounts,
}
