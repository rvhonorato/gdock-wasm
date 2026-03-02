// GA Parameters
pub const POPULATION_SIZE: u64 = 150;
pub const MUTATION_RATE: f64 = 0.1;
pub const CROSSOVER_RATE: f64 = 0.6;
pub const TOURNAMENT_SIZE: u64 = 3;
pub const ELITISM_COUNT: usize = 5; // Preserve top 5 individuals
pub const MAX_GENERATIONS: u64 = 250;
// NOTE: Changing the MAX_DISPLACEMENT from 15 to 50 greatly increased the score in debug from 0.7
// to 0.9!!
pub const MAX_DISPLACEMENT: f64 = 50.0;
pub const RANDOM_SEED: u64 = 42;

// Early stopping parameters
pub const ENABLE_EARLY_STOPPING: bool = true;
pub const CONVERGENCE_THRESHOLD: f64 = 0.001; // 0.1% improvement threshold
pub const CONVERGENCE_WINDOW: u64 = 10; // Stop if no improvement for this many generations

// Default weights
//  These defaults are based on a score calibration,
//   check the `calibration/` directory of the base
//   repository for more information
pub const DEFAULT_W_VDW: f64 = 0.4;
pub const DEFAULT_W_ELEC: f64 = 0.05;
pub const DEFAULT_W_DESOLV: f64 = 3.4;
pub const DEFAULT_W_AIR: f64 = 1.0;

/// Energy function weights for scoring
#[derive(Debug, Clone, Copy)]
pub struct EnergyWeights {
    pub vdw: f64,
    pub elec: f64,
    pub desolv: f64,
    pub air: f64,
}

impl Default for EnergyWeights {
    fn default() -> Self {
        Self {
            vdw: DEFAULT_W_VDW,
            elec: DEFAULT_W_ELEC,
            desolv: DEFAULT_W_DESOLV,
            air: DEFAULT_W_AIR,
        }
    }
}

impl EnergyWeights {
    pub fn new(vdw: f64, elec: f64, desolv: f64, air: f64) -> Self {
        Self {
            vdw,
            elec,
            desolv,
            air,
        }
    }
}

// Energy function cutoffs and parameters
pub const VDW_CUTOFF: f64 = 12.0; // Angstroms
pub const SOFTCORE_ALPHA: f64 = 0.5; // Softness parameter for soft-core LJ
pub const DESOLV_CUTOFF: f64 = 8.0; // Angstroms
pub const ELEC_CUTOFF: f64 = 15.0; // Angstroms
pub const ELEC_MIN_DISTANCE: f64 = 1.0; // Angstroms - avoid singularities
pub const AIR_FORCE_CONSTANT: f64 = 10.0; // kcal/mol per Å²

// VdW radii
// `https://en.wikipedia.org/wiki/Van_der_Waals_radius`
pub const HYDROGEN_RADIUS: f64 = 1.2;
pub const CARBON_RADIUS: f64 = 1.7;
pub const NITROGEN_RADIUS: f64 = 1.55;
pub const OXYGEN_RADIUS: f64 = 1.52;
pub const SULFUR_RADIUS: f64 = 1.8;
