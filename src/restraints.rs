use crate::structure;

#[derive(Debug, Clone)]
pub struct Restraint(pub structure::Atom, pub structure::Atom);

impl Restraint {
    fn new(atom1: structure::Atom, atom2: structure::Atom) -> Self {
        Self(atom1, atom2)
    }

    pub fn is_satisfied(
        &self,
        receptor: &structure::Molecule,
        ligand: &structure::Molecule,
    ) -> bool {
        let ca_receptor = receptor
            .0
            .iter()
            .find(|x| x.resseq == self.0.resseq && x.name.trim() == "CA");

        let ca_ligand = ligand
            .0
            .iter()
            .find(|x| x.resseq == self.1.resseq && x.name.trim() == "CA");

        if let (Some(ca1), Some(ca2)) = (ca_receptor, ca_ligand) {
            let dist = structure::distance(ca1, ca2);
            return dist <= 7.0;
        }

        false
    }
}

/// Create restraints from user-specified residue pairs
pub fn create_restraints_from_pairs(
    mol1: &structure::Molecule,
    mol2: &structure::Molecule,
    pairs: &[(i32, i32)],
) -> Vec<Restraint> {
    let mut restraints = Vec::new();

    for (res1, res2) in pairs {
        let ca_atom_1 = mol1
            .0
            .iter()
            .find(|atom| atom.name.trim() == "CA" && atom.resseq as i32 == *res1);

        let ca_atom_2 = mol2
            .0
            .iter()
            .find(|atom| atom.name.trim() == "CA" && atom.resseq as i32 == *res2);

        if let (Some(atom1), Some(atom2)) = (ca_atom_1, ca_atom_2) {
            restraints.push(Restraint::new(atom1.clone(), atom2.clone()));
        } else {
            eprintln!(
                "Warning: Restraint pair {}:{} not found in structures",
                res1, res2
            );
        }
    }

    restraints
}
