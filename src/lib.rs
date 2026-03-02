use serde::Serialize;
use wasm_bindgen::prelude::*;

mod chromosome;
mod clustering;
mod constants;
mod fitness;
mod ga;
mod hall_of_fame;
mod population;
mod restraints;
mod structure;
mod toppar;

const NUM_OUTPUT_MODELS: usize = 5;

#[wasm_bindgen]
pub struct ScoreResult {
    pub vdw: f64,
    pub elec: f64,
    pub desolv: f64,
    pub total: f64,
}

#[wasm_bindgen]
pub fn score_structures(
    receptor_pdb: &str,
    ligand_pdb: &str,
    w_vdw: f64,
    w_elec: f64,
    w_desolv: f64,
) -> Result<ScoreResult, JsValue> {
    let receptor_model = structure::read_pdb_from_str(receptor_pdb);
    let ligand_model = structure::read_pdb_from_str(ligand_pdb);

    let receptor = receptor_model
        .0
        .into_iter()
        .next()
        .ok_or_else(|| JsValue::from_str("Receptor PDB contains no atoms"))?;

    let ligand = ligand_model
        .0
        .into_iter()
        .next()
        .ok_or_else(|| JsValue::from_str("Ligand PDB contains no atoms"))?;

    let vdw = fitness::vdw_energy(&receptor, &ligand);
    let elec = fitness::elec_energy(&receptor, &ligand);
    let desolv = fitness::desolv_energy(&receptor, &ligand);
    let total = w_vdw * vdw + w_elec * elec + w_desolv * desolv;

    Ok(ScoreResult {
        vdw,
        elec,
        desolv,
        total,
    })
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PoseResult {
    rank: usize,
    fitness: f64,
    vdw: f64,
    elec: f64,
    desolv: f64,
    air: f64,
    cluster_size: usize,
    ligand_pdb: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct DockingOutput {
    generations_run: u32,
    /// Cluster centres sorted by fitness (diverse poses)
    clustered_poses: Vec<PoseResult>,
    /// Top entries from Hall of Fame sorted purely by fitness
    ranked_poses: Vec<PoseResult>,
}

#[wasm_bindgen]
pub fn run_docking(
    receptor_pdb: &str,
    ligand_pdb: &str,
    restraint_pairs_json: &str,
    max_generations: Option<u64>,
    seed: Option<u64>,
) -> Result<JsValue, JsValue> {
    use rand::SeedableRng;

    let receptor = structure::read_pdb_from_str(receptor_pdb)
        .0
        .into_iter()
        .next()
        .ok_or_else(|| JsValue::from_str("Receptor PDB contains no atoms"))?;

    let ligand = structure::read_pdb_from_str(ligand_pdb)
        .0
        .into_iter()
        .next()
        .ok_or_else(|| JsValue::from_str("Ligand PDB contains no atoms"))?;

    let pairs: Vec<(i32, i32)> = if restraint_pairs_json.trim().is_empty()
        || restraint_pairs_json.trim() == "[]"
    {
        Vec::new()
    } else {
        let parsed: Vec<Vec<i32>> = serde_json::from_str(restraint_pairs_json)
            .map_err(|e| JsValue::from_str(&format!("Invalid restraints JSON: {}", e)))?;
        parsed
            .into_iter()
            .filter(|p| p.len() == 2)
            .map(|p| (p[0], p[1]))
            .collect()
    };

    let restraint_list = restraints::create_restraints_from_pairs(&receptor, &ligand, &pairs);
    let weights = constants::EnergyWeights::default();

    let actual_seed = seed.unwrap_or(constants::RANDOM_SEED);
    let actual_max_generations = max_generations.unwrap_or(constants::MAX_GENERATIONS);

    let mut rng = rand::rngs::StdRng::seed_from_u64(actual_seed);

    let mut pop = population::Population::new(
        Vec::with_capacity(constants::POPULATION_SIZE as usize),
        receptor.clone(),
        ligand.clone(),
        restraint_list,
        weights,
    );

    for _ in 0..constants::POPULATION_SIZE {
        let c = chromosome::Chromosome::new(&mut rng);
        pop.chromosomes.push(c);
    }

    let mut hall_of_fame = hall_of_fame::HallOfFame::new();
    let mut generation_count = 0u64;
    let mut generations_without_improvement = 0u64;
    let mut last_best_score = f64::MAX;

    while generation_count < actual_max_generations {
        pop.eval_fitness();

        hall_of_fame.add_from_population(&pop.chromosomes);

        let best_fitness = pop.get_min_fittest().fitness;

        if generation_count > 0 {
            let improvement = (last_best_score - best_fitness) / last_best_score.abs();
            if improvement < constants::CONVERGENCE_THRESHOLD {
                generations_without_improvement += 1;
            } else {
                generations_without_improvement = 0;
            }
        }
        last_best_score = best_fitness;

        if constants::ENABLE_EARLY_STOPPING
            && generations_without_improvement >= constants::CONVERGENCE_WINDOW
        {
            break;
        }

        generation_count += 1;
        pop = pop.evolve(&mut rng);
    }

    // Final eval + final HoF collection
    pop.eval_fitness();
    hall_of_fame.add_from_population(&pop.chromosomes);

    let hof_entries = hall_of_fame.entries().to_vec();

    // Build receptor+ligand complexes from each HoF entry for FCC clustering
    let complexes: Vec<structure::Molecule> = hof_entries
        .iter()
        .map(|entry| {
            let docked_ligand = ligand
                .clone()
                .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
                .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
            structure::combine_molecules(&receptor, &docked_ligand)
        })
        .collect();

    // FCC clustering — cluster centres sorted by fitness (best first), like the CLI
    let cluster_config = clustering::ClusteringConfig::default();
    let clusters = clustering::cluster_structures(&complexes, &cluster_config);

    let mut cluster_centers: Vec<(usize, f64, usize)> = clusters
        .iter()
        .map(|c| (c.center_idx, hof_entries[c.center_idx].fitness, c.size))
        .collect();
    cluster_centers.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut clustered_sel: Vec<(usize, usize)> = cluster_centers
        .iter()
        .take(NUM_OUTPUT_MODELS)
        .map(|(idx, _, size)| (*idx, *size))
        .collect();

    let mut used: std::collections::HashSet<usize> = clustered_sel.iter().map(|(i, _)| *i).collect();

    // Fill to NUM_OUTPUT_MODELS with best HoF entries not yet selected
    if clustered_sel.len() < NUM_OUTPUT_MODELS {
        let mut remaining: Vec<(usize, f64)> = hof_entries
            .iter()
            .enumerate()
            .filter(|(i, _)| !used.contains(i))
            .map(|(i, e)| (i, e.fitness))
            .collect();
        remaining.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (idx, _) in remaining.iter().take(NUM_OUTPUT_MODELS - clustered_sel.len()) {
            clustered_sel.push((*idx, 1));
            used.insert(*idx);
        }
    }

    // Ranked by score: top NUM_OUTPUT_MODELS from HoF sorted purely by fitness
    let mut ranked_by_fitness: Vec<(usize, f64)> = hof_entries
        .iter()
        .enumerate()
        .map(|(i, e)| (i, e.fitness))
        .collect();
    ranked_by_fitness.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    fn build_pose(rank: usize, hof_idx: usize, cluster_size: usize, hof_entries: &[hall_of_fame::HallOfFameEntry], ligand: &structure::Molecule) -> PoseResult {
        let entry = &hof_entries[hof_idx];
        let docked = ligand
            .clone()
            .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
            .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
        PoseResult { rank, fitness: entry.fitness, vdw: entry.vdw, elec: entry.elec, desolv: entry.desolv, air: entry.air, cluster_size, ligand_pdb: docked.to_pdb_string() }
    }

    let clustered_poses: Vec<PoseResult> = clustered_sel
        .iter()
        .enumerate()
        .map(|(i, (hof_idx, cluster_size))| build_pose(i + 1, *hof_idx, *cluster_size, &hof_entries, &ligand))
        .collect();

    let ranked_poses: Vec<PoseResult> = ranked_by_fitness
        .iter()
        .take(NUM_OUTPUT_MODELS)
        .enumerate()
        .map(|(i, (hof_idx, _))| build_pose(i + 1, *hof_idx, 0, &hof_entries, &ligand))
        .collect();

    let output = DockingOutput {
        generations_run: generation_count as u32,
        clustered_poses,
        ranked_poses,
    };

    serde_wasm_bindgen::to_value(&output).map_err(|e| JsValue::from_str(&e.to_string()))
}
