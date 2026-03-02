//! FCC (Fraction of Common Contacts) clustering — sequential version for WASM.

use std::collections::{HashMap, HashSet};

use crate::structure::Molecule;

#[derive(Debug, Clone)]
pub struct ClusteringConfig {
    pub contact_distance: f64,
    pub fcc_cutoff: f64,
    pub strictness: f64,
    pub min_cluster_size: usize,
}

impl Default for ClusteringConfig {
    fn default() -> Self {
        Self {
            contact_distance: 5.0,
            fcc_cutoff: 0.60,
            strictness: 0.75,
            min_cluster_size: 4,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ClusterResult {
    pub center_idx: usize,
    pub members: Vec<usize>,
    pub size: usize,
}

struct Element {
    neighbors: HashSet<usize>,
}

type ResidueCoordMap = HashMap<(char, i16), Vec<(f64, f64, f64)>>;

pub fn calculate_contacts(molecule: &Molecule, contact_distance: f64) -> HashSet<String> {
    let mut residues: ResidueCoordMap = HashMap::new();

    for atom in &molecule.0 {
        if atom.element.trim() == "H" {
            continue;
        }
        residues
            .entry((atom.chainid, atom.resseq))
            .or_default()
            .push((atom.x, atom.y, atom.z));
    }

    let mut contacts = HashSet::new();
    let residue_keys: Vec<_> = residues.keys().collect();

    for i in 0..residue_keys.len() {
        for j in (i + 1)..residue_keys.len() {
            let (chain_a, res_a) = residue_keys[i];
            let (chain_b, res_b) = residue_keys[j];

            if chain_a == chain_b {
                continue;
            }

            let atoms_a = &residues[residue_keys[i]];
            let atoms_b = &residues[residue_keys[j]];

            'outer: for (xa, ya, za) in atoms_a {
                for (xb, yb, zb) in atoms_b {
                    let dist =
                        ((xa - xb).powi(2) + (ya - yb).powi(2) + (za - zb).powi(2)).sqrt();
                    if dist <= contact_distance {
                        let contact = if chain_a < chain_b {
                            format!("{} {} {} {}", chain_a, res_a, chain_b, res_b)
                        } else {
                            format!("{} {} {} {}", chain_b, res_b, chain_a, res_a)
                        };
                        contacts.insert(contact);
                        break 'outer;
                    }
                }
            }
        }
    }

    contacts
}

fn calculate_fcc(x: &HashSet<String>, y: &HashSet<String>) -> (f64, f64) {
    if x.is_empty() || y.is_empty() {
        return (0.0, 0.0);
    }
    let common = x.intersection(y).count() as f64;
    (common / x.len() as f64, common / y.len() as f64)
}

fn calculate_pairwise_fcc(contact_sets: &[HashSet<String>]) -> Vec<(usize, usize, f64, f64)> {
    let n = contact_sets.len();
    let mut results = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let (fcc, fcc_v) = calculate_fcc(&contact_sets[i], &contact_sets[j]);
            results.push((i, j, fcc, fcc_v));
        }
    }
    results
}

fn create_elements(
    pairwise_fcc: Vec<(usize, usize, f64, f64)>,
    n_structures: usize,
    config: &ClusteringConfig,
) -> HashMap<usize, Element> {
    let mut elements: HashMap<usize, Element> = (0..n_structures)
        .map(|i| (i, Element { neighbors: HashSet::new() }))
        .collect();

    for (i, j, fcc, fcc_v) in pairwise_fcc {
        if fcc >= config.fcc_cutoff && fcc_v >= config.fcc_cutoff * config.strictness {
            elements.get_mut(&i).unwrap().neighbors.insert(j);
        }
        if fcc_v >= config.fcc_cutoff && fcc >= config.fcc_cutoff * config.strictness {
            elements.get_mut(&j).unwrap().neighbors.insert(i);
        }
    }

    elements
}

fn cluster_elements(mut elements: HashMap<usize, Element>) -> Vec<ClusterResult> {
    let mut used: HashSet<usize> = HashSet::new();
    let mut clusters = Vec::new();

    loop {
        let clusterable: Vec<usize> = elements
            .keys()
            .filter(|k| !used.contains(k))
            .copied()
            .collect();

        if clusterable.is_empty() {
            break;
        }

        let center = clusterable
            .iter()
            .max_by(|&&a, &&b| {
                let ca = elements[&a].neighbors.iter().filter(|n| !used.contains(n)).count();
                let cb = elements[&b].neighbors.iter().filter(|n| !used.contains(n)).count();
                ca.cmp(&cb).then_with(|| b.cmp(&a))
            })
            .copied()
            .unwrap();

        let neighbors: Vec<usize> = elements[&center]
            .neighbors
            .iter()
            .filter(|n| !used.contains(n))
            .copied()
            .collect();

        let mut members: Vec<usize> = vec![center];
        for neighbor in &neighbors {
            if *neighbor != center {
                members.push(*neighbor);
                used.insert(*neighbor);
            }
        }
        used.insert(center);
        members.sort();

        clusters.push(ClusterResult {
            center_idx: center,
            members: members.clone(),
            size: members.len(),
        });

        elements.remove(&center);
    }

    clusters.sort_by(|a, b| b.size.cmp(&a.size));
    clusters
}

pub fn cluster_structures(structures: &[Molecule], config: &ClusteringConfig) -> Vec<ClusterResult> {
    if structures.is_empty() {
        return Vec::new();
    }

    let contact_sets: Vec<HashSet<String>> = structures
        .iter()
        .map(|mol| calculate_contacts(mol, config.contact_distance))
        .collect();

    let pairwise_fcc = calculate_pairwise_fcc(&contact_sets);
    let elements = create_elements(pairwise_fcc, structures.len(), config);
    let clusters = cluster_elements(elements);

    clusters
        .into_iter()
        .filter(|c| c.size >= config.min_cluster_size)
        .collect()
}
