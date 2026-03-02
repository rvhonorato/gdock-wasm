use crate::chromosome::Chromosome;
use std::f64::consts::PI;

const HALL_OF_FAME_MAX_SIZE: usize = 500;
const HALL_OF_FAME_TOP_K: usize = 10;
const UNIQUENESS_ROTATION_THRESHOLD: f64 = 0.2; // ~11 degrees
const UNIQUENESS_TRANSLATION_THRESHOLD: f64 = 2.0; // 2 Å

#[derive(Debug, Clone)]
pub struct HallOfFameEntry {
    pub genes: [f64; 6],
    pub fitness: f64,
    pub vdw: f64,
    pub elec: f64,
    pub desolv: f64,
    pub air: f64,
}

#[derive(Debug)]
pub struct HallOfFame {
    entries: Vec<HallOfFameEntry>,
    max_size: usize,
}

impl HallOfFame {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            max_size: HALL_OF_FAME_MAX_SIZE,
        }
    }

    pub fn try_add(&mut self, genes: &[f64], fitness: f64, vdw: f64, elec: f64, desolv: f64, air: f64) -> bool {
        if genes.len() != 6 {
            return false;
        }
        let new_genes: [f64; 6] = [genes[0], genes[1], genes[2], genes[3], genes[4], genes[5]];
        if !self.is_unique(&new_genes) {
            return false;
        }
        self.entries.push(HallOfFameEntry { genes: new_genes, fitness, vdw, elec, desolv, air });
        if self.entries.len() > self.max_size {
            self.prune();
        }
        true
    }

    fn is_unique(&self, new_genes: &[f64; 6]) -> bool {
        for entry in &self.entries {
            if Self::genes_are_similar(new_genes, &entry.genes) {
                return false;
            }
        }
        true
    }

    fn genes_are_similar(a: &[f64; 6], b: &[f64; 6]) -> bool {
        for i in 0..3 {
            if Self::angular_difference(a[i], b[i]) > UNIQUENESS_ROTATION_THRESHOLD {
                return false;
            }
        }
        for i in 3..6 {
            if (a[i] - b[i]).abs() > UNIQUENESS_TRANSLATION_THRESHOLD {
                return false;
            }
        }
        true
    }

    fn angular_difference(a: f64, b: f64) -> f64 {
        let diff = (a - b).abs();
        if diff > PI { 2.0 * PI - diff } else { diff }
    }

    fn prune(&mut self) {
        self.entries
            .sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap());
        self.entries.truncate(self.max_size);
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn entries(&self) -> &[HallOfFameEntry] {
        &self.entries
    }

    pub fn add_from_population(&mut self, chromosomes: &[Chromosome]) {
        let mut indices: Vec<usize> = (0..chromosomes.len()).collect();
        indices.sort_by(|&a, &b| {
            chromosomes[a]
                .fitness
                .partial_cmp(&chromosomes[b].fitness)
                .unwrap()
        });

        let mut added = 0;
        for &idx in &indices {
            if added >= HALL_OF_FAME_TOP_K {
                break;
            }
            let chr = &chromosomes[idx];
            if self.try_add(&chr.genes, chr.fitness, chr.vdw, chr.elec, chr.desolv, chr.air) {
                added += 1;
            }
        }
    }
}
