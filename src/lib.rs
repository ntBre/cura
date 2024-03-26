use openff_toolkit::ForceField;
use rdkit_rs::ROMol;

pub const PROGRESS_INTERVAL: usize = 1000;

pub mod parse;
pub mod query;
pub mod store;
pub mod table;

/// Load an OpenFF [ForceField] from `forcefield` and return a sequence of
/// parameter_id, SMIRKS pattern pairs corresponding to its `parameter_type`
/// [ParameterHandler].
fn load_forcefield(
    forcefield: String,
    parameter_type: String,
) -> Vec<(String, ROMol)> {
    ForceField::load(&forcefield)
        .unwrap()
        .get_parameter_handler(&parameter_type)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), ROMol::from_smarts(&p.smirks())))
        .collect()
}
