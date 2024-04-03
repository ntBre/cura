use std::sync::atomic::AtomicUsize;

use log::{trace, warn};
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::{RDError, ROMol, SDMolSupplier};

use crate::{get_elements, table::Table, Molecule, PROGRESS_INTERVAL};

pub fn store(table: &mut Table, molecule_file: &str) {
    trace!("initializing mol supplier from {}", molecule_file);
    let m = SDMolSupplier::new(molecule_file).unwrap();

    let progress = AtomicUsize::new(0);

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
            .map(|s| {
                let mut mol = ROMol::from_smiles(s);
                mol.openff_clean();
                Molecule::new(
                    mol.to_smiles(),
                    mol.num_atoms(),
                    get_elements(&mol),
                )
            })
            .collect()
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    table.insert_molecules(results).unwrap();
}
