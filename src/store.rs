use std::sync::atomic::AtomicUsize;

use log::{trace, warn};
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::{fragment::recap_decompose, RDError, ROMol, SDMolSupplier};

use crate::{get_elements, table::Table, Molecule, PROGRESS_INTERVAL};

pub fn store(table: &mut Table, molecule_file: &str, tag: String) {
    trace!("initializing mol supplier from {}", molecule_file);
    let m = SDMolSupplier::new(molecule_file).unwrap();

    let progress = AtomicUsize::new(0);

    // these are from Lily's recap example script
    let so_inp = ROMol::from_smiles("S(=O)(=O)*");
    let so_out = ROMol::from_smiles("S(=O)([O-])");

    let map_op = |mol: Result<ROMol, RDError>| -> Vec<Molecule> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let Ok(mut mol) = mol else {
            warn!("error loading molecule, skipping");
            return Vec::new();
        };
        mol.openff_clean();
        mol.to_smiles()
            .split('.')
            .flat_map(|s| {
                let mut mol = ROMol::from_smiles(s);
                mol.openff_clean();
                recap_decompose(&mol, None, Some(4), None)
                    .get_leaves()
                    .into_values()
                    .map(|mut m| {
                        m = m
                            .replace_substructs(&so_inp, &so_out, true)
                            .remove(0);
                        Molecule::new(
                            m.to_smiles(),
                            None,
                            m.num_atoms(),
                            get_elements(&m),
                            tag.clone(),
                        )
                    })
                    .chain(std::iter::once(Molecule::new(
                        mol.to_smiles(),
                        Some(mol.to_inchi_key()),
                        mol.num_atoms(),
                        get_elements(&mol),
                        tag.clone(),
                    )))
            })
            .collect()
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    table.insert_molecules(results).unwrap();
}
