use wasm_bindgen::prelude::*;

mod constants;
mod fitness;
mod restraints;
mod structure;
mod toppar;

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
