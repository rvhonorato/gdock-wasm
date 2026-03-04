use gdock::chromosome;
use gdock::clustering;
use gdock::constants;
use gdock::fitness;
use gdock::population;
use gdock::restraints;
use gdock::runner::{run_ga, select_models};
use gdock::structure;
use serde::Serialize;
use wasm_bindgen::prelude::*;

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

    let pairs: Vec<(i32, i32)> =
        if restraint_pairs_json.trim().is_empty() || restraint_pairs_json.trim() == "[]" {
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
        structure::Molecule::new(), // no reference structure in wasm
        restraint_list,
        weights,
        None, // no debug evaluator in wasm
    );

    for _ in 0..constants::POPULATION_SIZE {
        let c = chromosome::Chromosome::new(&mut rng);
        pop.chromosomes.push(c);
    }

    let ga_result = run_ga(pop, &mut rng, actual_max_generations, |_, _| {});
    let hof_entries = ga_result.hall_of_fame.entries().to_vec();
    let generation_count = ga_result.generations_run;

    let selected = select_models(&hof_entries, &receptor, &ligand);

    fn build_pose(
        rank: usize,
        hof_idx: usize,
        cluster_size: usize,
        hof_entries: &[gdock::hall_of_fame::HallOfFameEntry],
        ligand: &structure::Molecule,
    ) -> PoseResult {
        let entry = &hof_entries[hof_idx];
        let docked = ligand
            .clone()
            .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
            .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
        PoseResult {
            rank,
            fitness: entry.fitness,
            vdw: entry.vdw,
            elec: entry.elec,
            desolv: entry.desolv,
            air: entry.air,
            cluster_size,
            ligand_pdb: docked.to_pdb_string(),
        }
    }

    let clustered_poses: Vec<PoseResult> = selected
        .clustered
        .iter()
        .enumerate()
        .map(|(i, (hof_idx, cluster_size))| {
            build_pose(i + 1, *hof_idx, *cluster_size, &hof_entries, &ligand)
        })
        .collect();

    let ranked_poses: Vec<PoseResult> = selected
        .ranked
        .iter()
        .enumerate()
        .map(|(i, hof_idx)| build_pose(i + 1, *hof_idx, 0, &hof_entries, &ligand))
        .collect();

    let output = DockingOutput {
        generations_run: generation_count as u32,
        clustered_poses,
        ranked_poses,
    };

    serde_wasm_bindgen::to_value(&output).map_err(|e| JsValue::from_str(&e.to_string()))
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ClusterOutput {
    center_idx: usize,
    members: Vec<usize>,
    size: usize,
}

/// Compute inter-chain FCC contacts for a receptor + ligand complex.
///
/// Returns a JS array of contact strings ("chainA resA chainB resB").
/// Call this once per model while iterating over the ensemble.
#[wasm_bindgen]
pub fn compute_contacts(receptor_pdb: &str, ligand_pdb: &str) -> Result<JsValue, JsValue> {
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

    let complex = structure::combine_molecules(&receptor, &ligand);
    let mut contacts: Vec<String> = clustering::calculate_contacts(&complex, 5.0)
        .into_iter()
        .collect();
    contacts.sort();

    serde_wasm_bindgen::to_value(&contacts).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Cluster structures from pre-computed contact sets using FCC.
///
/// `contacts_json` — JSON array of arrays of contact strings, one inner array
/// per model (as returned by repeated calls to `compute_contacts`).
///
/// Returns an array of cluster objects `{ centerIdx, members, size }` sorted
/// by cluster size descending. Only clusters with `size >= min_cluster_size`
/// are included.
#[wasm_bindgen]
pub fn cluster_from_contacts(
    contacts_json: &str,
    fcc_cutoff: Option<f64>,
    min_cluster_size: Option<usize>,
) -> Result<JsValue, JsValue> {
    use std::collections::HashSet;

    let raw: Vec<Vec<String>> = serde_json::from_str(contacts_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid JSON: {}", e)))?;

    let contact_sets: Vec<HashSet<String>> =
        raw.into_iter().map(|v| v.into_iter().collect()).collect();

    let n = contact_sets.len();
    let cutoff = fcc_cutoff.unwrap_or(0.60);
    let strictness = 0.75_f64;
    let min_size = min_cluster_size.unwrap_or(4);

    // Build FCC neighbor graph (O(n²))
    let mut neighbors: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let xi = &contact_sets[i];
            let xj = &contact_sets[j];
            if xi.is_empty() || xj.is_empty() {
                continue;
            }
            let common = xi.intersection(xj).count() as f64;
            let fcc = common / xi.len() as f64;
            let fcc_v = common / xj.len() as f64;
            if fcc >= cutoff && fcc_v >= cutoff * strictness {
                neighbors[i].insert(j);
            }
            if fcc_v >= cutoff && fcc >= cutoff * strictness {
                neighbors[j].insert(i);
            }
        }
    }

    // Greedy clustering (mirrors cluster_elements in gdock::clustering)
    let mut used: HashSet<usize> = HashSet::new();
    let mut clusters: Vec<ClusterOutput> = Vec::new();

    loop {
        let center = (0..n).filter(|k| !used.contains(k)).max_by(|&a, &b| {
            let ca = neighbors[a].iter().filter(|x| !used.contains(x)).count();
            let cb = neighbors[b].iter().filter(|x| !used.contains(x)).count();
            ca.cmp(&cb).then_with(|| b.cmp(&a)) // tie-break: larger index wins
        });

        let Some(center) = center else { break };

        let mut members: Vec<usize> = vec![center];
        for &nb in &neighbors[center] {
            if !used.contains(&nb) && nb != center {
                members.push(nb);
                used.insert(nb);
            }
        }
        used.insert(center);
        members.sort();

        clusters.push(ClusterOutput {
            center_idx: center,
            size: members.len(),
            members,
        });
    }

    clusters.sort_by(|a, b| b.size.cmp(&a.size));
    let output: Vec<ClusterOutput> = clusters
        .into_iter()
        .filter(|c| c.size >= min_size)
        .collect();

    serde_wasm_bindgen::to_value(&output).map_err(|e| JsValue::from_str(&e.to_string()))
}
