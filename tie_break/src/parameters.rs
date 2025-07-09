//! This is the data structure required for the parameters yaml file

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Params {
    pub round_1: FilterParams,
    pub round_2: FilterParams,
    pub min_unique_coverage: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, Default)]
pub struct FilterParams {
    pub allow_reads_multiple_alns_per_ref: bool,
    pub min_query_coverage: Option<f32>,
    pub min_query_coverage_pc_of_max: Option<f32>,
    pub max_divergence: Option<f32>,
    pub max_divergence_from_expected: Option<f32>,
    pub max_divergence_from_best: Option<f32>,
    pub min_score_fraction: Option<f32>,
    pub share_threshold: Option<f32>,
}
