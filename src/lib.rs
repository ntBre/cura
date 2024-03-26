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
        let got = table.get_moldata().unwrap().len();
        let want = 1;
        assert_eq!(got, want);
    }
}
