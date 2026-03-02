use crate::constants;
use crate::toppar;

use std::collections::HashSet;

#[cfg(not(target_arch = "wasm32"))]
use std::fs::File;
#[cfg(not(target_arch = "wasm32"))]
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Model(pub Vec<Molecule>);

impl Default for Model {
    fn default() -> Self {
        Self::new()
    }
}

impl Model {
    pub fn new() -> Model {
        Model(Vec::new())
    }
}

#[derive(Debug, Clone)]
pub struct Molecule(pub Vec<Atom>);

#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: i32,
    pub name: String,
    pub altloc: char,
    pub resname: String,
    pub chainid: char,
    pub resseq: i16,
    pub icode: char,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f32,
    pub tempfactor: f32,
    pub element: String,
    pub charge: f64,
    pub vdw_radius: f64,
    pub epsilon: f64,
    pub rmin2: f64,
    pub eps_1_4: f64,
    pub rmin2_1_4: f64,
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

impl Molecule {
    pub fn new() -> Molecule {
        Molecule(Vec::new())
    }

    pub fn center_of_mass(&self) -> (f64, f64, f64) {
        let num_atoms = self.0.len() as f64;
        let sum_x = self.0.iter().map(|atom| atom.x).sum::<f64>();
        let sum_y = self.0.iter().map(|atom| atom.y).sum::<f64>();
        let sum_z = self.0.iter().map(|atom| atom.z).sum::<f64>();
        (sum_x / num_atoms, sum_y / num_atoms, sum_z / num_atoms)
    }

    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        for atom in &mut self.0 {
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
        }
    }
}

impl Atom {
    fn new() -> Atom {
        Atom {
            serial: 0,
            name: String::new(),
            altloc: ' ',
            resname: String::new(),
            chainid: ' ',
            resseq: 0,
            icode: ' ',
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 0.0,
            tempfactor: 0.0,
            element: String::new(),
            charge: 0.0,
            vdw_radius: 0.0,
            epsilon: 0.0,
            rmin2: 0.0,
            eps_1_4: 0.0,
            rmin2_1_4: 0.0,
        }
    }
}

/// Parse an ATOM line from a PDB file using fixed-column format
fn process_atom_line(line: &str) -> Option<Atom> {
    if !line.starts_with("ATOM") || line.len() < 54 {
        return None;
    }

    let get_str = |start: usize, end: usize| -> &str {
        if line.len() >= end {
            &line[start..end]
        } else {
            ""
        }
    };

    let mut atom = Atom::new();

    atom.serial = get_str(6, 11).trim().parse().ok()?;
    atom.name = get_str(12, 16).trim().to_string();
    atom.altloc = get_str(16, 17).chars().next().unwrap_or(' ');
    atom.resname = get_str(17, 20).trim().to_string();
    atom.chainid = get_str(21, 22).chars().next().unwrap_or(' ');
    atom.resseq = get_str(22, 26).trim().parse().ok()?;
    atom.icode = get_str(26, 27).chars().next().unwrap_or(' ');
    atom.x = get_str(30, 38).trim().parse().ok()?;
    atom.y = get_str(38, 46).trim().parse().ok()?;
    atom.z = get_str(46, 54).trim().parse().ok()?;
    atom.occupancy = get_str(54, 60).trim().parse().unwrap_or(1.0);
    atom.tempfactor = get_str(60, 66).trim().parse().unwrap_or(0.0);
    atom.element = get_str(76, 78).trim().to_string();

    if atom.element.is_empty() {
        atom.element = atom.name.chars().next().unwrap_or(' ').to_string();
    }

    atom.vdw_radius = match atom.element.trim() {
        "H" => constants::HYDROGEN_RADIUS,
        "C" => constants::CARBON_RADIUS,
        "N" => constants::NITROGEN_RADIUS,
        "O" => constants::OXYGEN_RADIUS,
        _ => 1.0,
    };

    let atom_type = toppar::get_atom(atom.resname.as_str(), atom.name.as_str());

    if let Some(v) = atom_type {
        atom.epsilon = toppar::get_epsilon(v).unwrap_or(0.0);
        atom.rmin2 = toppar::get_rmin2(v).unwrap_or(0.0);
        atom.eps_1_4 = toppar::get_eps_1_4(v).unwrap_or(0.0);
        atom.rmin2_1_4 = toppar::get_rmin2_1_4(v).unwrap_or(0.0);
        atom.charge = toppar::get_charge(v).unwrap_or(0.0);
    }

    Some(atom)
}

fn parse_lines<I>(lines: I) -> Model
where
    I: Iterator<Item = String>,
{
    let mut model = Model::new();
    let mut molecule = Molecule::new();
    let mut has_model_records = false;

    for line in lines {
        if line.starts_with("MODEL") {
            has_model_records = true;
            molecule = Molecule::new();
        } else if line.starts_with("ENDMDL") {
            model.0.push(molecule.clone());
            molecule = Molecule::new();
        } else if let Some(atom) = process_atom_line(&line) {
            molecule.0.push(atom);
        }
    }

    if !has_model_records {
        model.0.push(molecule);
    }

    model
}

/// Reads a PDB from an in-memory string. Works in both native and WASM targets.
pub fn read_pdb_from_str(content: &str) -> Model {
    parse_lines(content.lines().map(|l| l.to_string()))
}

/// Reads a PDB file from disk. Not available in WASM.
#[cfg(not(target_arch = "wasm32"))]
pub fn read_pdb(pdb_file: &str) -> Model {
    let file = File::open(pdb_file).expect("Cannot open file");
    parse_lines(BufReader::new(file).lines().map_while(Result::ok))
}

pub fn distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let dx = atom1.x - atom2.x;
    let dy = atom1.y - atom2.y;
    let dz = atom1.z - atom2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

pub fn filter_by_resseq_vec(molecule: &Molecule, resseq_vec: &HashSet<i16>) -> Molecule {
    let mut filtered_molecule = Molecule::new();
    for atom in &molecule.0 {
        if resseq_vec.contains(&atom.resseq) {
            filtered_molecule.0.push(atom.clone());
        }
    }
    filtered_molecule
}
