use std::sync::atomic::AtomicUsize;

use log::{trace, warn};
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::{RDError, ROMol, SDMolSupplier};

use crate::{get_elements, table::Table, Molecule, PROGRESS_INTERVAL};

pub fn store(table: &mut Table, molecule_file: &str) {
    trace!("initializing mol supplier from {}", molecule_file);
    let m = SDMolSupplier::new(molecule_file).unwrap();

    let progress = AtomicUsize::new(0);

    let map_op = |mol: Result<ROMol, RDError>| -> Option<Molecule> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let Ok(mut mol) = mol else {
            warn!("error loading molecule, skipping");
            return None;
        };
        mol.openff_clean();

        Some(Molecule {
            smiles: mol.to_smiles(),
            natoms: mol.num_atoms(),
            elements: get_elements(&mol),
            moldata: mol.to_json(),
        })
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    table.insert_molecules(results).unwrap();
}
