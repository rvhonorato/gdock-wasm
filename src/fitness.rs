use crate::constants::{
    AIR_FORCE_CONSTANT, DESOLV_CUTOFF, ELEC_CUTOFF, ELEC_MIN_DISTANCE, SOFTCORE_ALPHA, VDW_CUTOFF,
};
use crate::restraints;
use crate::structure;
use crate::structure::Atom;

fn softcore_lj_potential(atom1: &Atom, atom2: &Atom, distance: f64, alpha: f64) -> f64 {
    let epsilon_ij = (atom1.epsilon * atom2.epsilon).sqrt();
    let rmin_ij = atom1.rmin2 + atom2.rmin2;

    let rmin6 = rmin_ij.powi(6);
    let r6 = distance.powi(6);
    let r_eff6 = r6 + alpha * rmin6;

    let ratio = rmin6 / r_eff6;
    epsilon_ij * (ratio.powi(2) - 2.0 * ratio)
}

pub fn vdw_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let mut energy = 0.0;

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            if atom1.serial != atom2.serial {
                let dist = structure::distance(atom1, atom2);
                if dist > 0.0 && dist < VDW_CUTOFF {
                    let vdw = softcore_lj_potential(atom1, atom2, dist, SOFTCORE_ALPHA);
                    let capped_vdw = vdw.clamp(-50.0, 500.0);
                    energy += capped_vdw;
                }
            }
        }
    }
    energy
}

pub fn desolv_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let mut energy = 0.0;

    let get_asp = |element: &str| -> f64 {
        match element.trim() {
            "C" => 0.012,
            "N" => -0.160,
            "O" => -0.160,
            "S" => 0.012,
            "H" => 0.0,
            _ => 0.0,
        }
    };

    for atom1 in &receptor.0 {
        let asp1 = get_asp(&atom1.element);
        if asp1.abs() < 0.001 {
            continue;
        }

        let mut burial_count = 0;
        for atom2 in &ligand.0 {
            let dist = structure::distance(atom1, atom2);
            if dist < DESOLV_CUTOFF {
                burial_count += 1;
            }
        }

        let burial = (burial_count as f64 / 10.0).min(1.0);
        energy -= asp1 * burial * atom1.vdw_radius * atom1.vdw_radius * 4.0 * std::f64::consts::PI;
    }

    for atom2 in &ligand.0 {
        let asp2 = get_asp(&atom2.element);
        if asp2.abs() < 0.001 {
            continue;
        }

        let mut burial_count = 0;
        for atom1 in &receptor.0 {
            let dist = structure::distance(atom1, atom2);
            if dist < DESOLV_CUTOFF {
                burial_count += 1;
            }
        }

        let burial = (burial_count as f64 / 10.0).min(1.0);
        energy -= asp2 * burial * atom2.vdw_radius * atom2.vdw_radius * 4.0 * std::f64::consts::PI;
    }

    energy
}

pub fn air_energy(
    restraints: &[crate::restraints::Restraint],
    receptor: &structure::Molecule,
    ligand: &structure::Molecule,
) -> f64 {
    let mut energy = 0.0;

    let lower_bound = 0.0;
    let upper_bound = 7.0;

    for restraint in restraints {
        let ca_receptor = receptor
            .0
            .iter()
            .find(|a| a.resseq == restraint.0.resseq && a.name.trim() == "CA");

        let ca_ligand = ligand
            .0
            .iter()
            .find(|a| a.resseq == restraint.1.resseq && a.name.trim() == "CA");

        if let (Some(ca1), Some(ca2)) = (ca_receptor, ca_ligand) {
            let dist = structure::distance(ca1, ca2);

            let violation = if dist < lower_bound {
                lower_bound - dist
            } else if dist > upper_bound {
                dist - upper_bound
            } else {
                0.0
            };

            energy += AIR_FORCE_CONSTANT * violation * violation;
        }
    }

    energy
}

pub fn elec_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let k = 332.0636;
    let mut energy = 0.0;

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            let dist = structure::distance(atom1, atom2);
            if dist > ELEC_MIN_DISTANCE && dist < ELEC_CUTOFF {
                let charge1 = atom1.charge;
                let charge2 = atom2.charge;
                energy += k * (charge1 * charge2) / (dist * dist);
            }
        }
    }

    energy
}

pub fn satisfaction_ratio(
    restraints: &[restraints::Restraint],
    receptor: &structure::Molecule,
    ligand: &structure::Molecule,
) -> f64 {
    let satisfied = restraints
        .iter()
        .filter(|&x| x.is_satisfied(receptor, ligand))
        .count();

    let total = restraints.len();
    if total > 0 {
        satisfied as f64 / total as f64
    } else {
        0.0
    }
}
