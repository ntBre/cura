use std::collections::{HashMap, HashSet};

use openff_toolkit::ForceField;
use rdkit_rs::{find_smarts_matches_mol, ROMol};

pub const PROGRESS_INTERVAL: usize = 1000;

pub mod board;
pub mod parse;
pub mod query;
pub mod store;
pub mod table;

/// Parameter ID in a force field
type Pid = String;

/// SMIRKS string
type Smirks = String;

/// Database record
pub struct Molecule {
    /// None when inserting into the database but set when retrieving
    id: Option<usize>,
    smiles: String,
    natoms: usize,
    elements: i128,
}

impl Molecule {
    pub fn new(smiles: String, natoms: usize, elements: i128) -> Self {
        Self {
            id: None,
            smiles,
            natoms,
            elements,
        }
    }
}

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

/// returns a map of chemical environments in `mol` to their matching parameter
/// ids. matching starts with the first parameter and proceeds through the
/// whole sequence of parameters, so this should follow the SMIRNOFF typing
/// rules
pub fn find_matches_full(
    params: &[(String, ROMol)],
    mol: &ROMol,
) -> HashMap<Vec<usize>, String> {
    let mut matches = HashMap::new();
    for (id, smirks) in params {
        let env_matches = find_smarts_matches_mol(mol, smirks);
        for mut mat in env_matches {
            if mat.first().unwrap() > mat.last().unwrap() {
                mat.reverse();
            }
            matches.insert(mat, id.clone());
        }
    }
    matches
}

/// returns the set of parameter ids matching `mol`. matching starts with the
/// first parameter and proceeds through the whole sequence of parameters, so
/// this should follow the SMIRNOFF typing rules
pub fn find_matches(
    params: &[(String, ROMol)],
    mol: &ROMol,
) -> HashSet<String> {
    find_matches_full(params, mol).into_values().collect()
}

/// Encode the elements in `mol` as a 128-bit set
pub fn get_elements(mol: &ROMol) -> i128 {
    let mut ret = 0;
    for i in mol.elements() {
        ret |= 1 << i;
    }
    ret
}

pub fn bits_to_elements(bits: i128) -> Vec<usize> {
    let mut ret = Vec::new();
    for i in 0..128 {
        if (bits & (1 << i)) != 0 {
            ret.push(i);
        }
    }
    ret
}

#[cfg(test)]
mod tests {
    use self::{store::store, table::Table};

    use super::*;

    #[test]
    fn test_get_elements() {
        let mol = ROMol::from_smiles("CCO");
        let got = bits_to_elements(get_elements(&mol));
        let want = vec![6, 8];
        assert_eq!(got, want);
    }

    #[test]
    fn test_store() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut table = Table::create(tmp.path()).unwrap();
        // this file has multiple entries, but only one SMILES
        store(&mut table, "testfiles/small.sdf");
        let got = table.get_smiles().unwrap().len();
        let want = 1;
        assert_eq!(got, want);

        let res = table.with_molecules(|mol| mol.natoms).unwrap();
        assert_eq!(res.len(), 1);
    }
}
