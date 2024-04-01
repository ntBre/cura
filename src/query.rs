use std::{
    collections::{HashMap, HashSet},
    path::Path,
    sync::{atomic::AtomicUsize, mpsc},
    thread,
};

use log::{info, trace};
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::{find_smarts_matches_mol, ROMol};
use rsearch::{load_want, print_output, write_output};

use crate::{load_forcefield, table::Table, PROGRESS_INTERVAL};

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
            trace!("{mat:?} => {id}");
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

pub fn query(
    table: Table,
    forcefield: String,
    parameter_type: String,
    search_params: String,
    output_dir: Option<String>,
) {
    info!("loading moldata from database");
    const CHAN_SIZE: usize = 1024;
    let (sender, receiver) = mpsc::sync_channel::<String>(CHAN_SIZE);
    let th = thread::spawn(move || table.send_smiles(sender));

    info!("loading force field and parameters");
    let params = load_forcefield(forcefield, parameter_type);

    let want = load_want(&search_params);

    info!("processing data");
    let progress = AtomicUsize::new(0);
    let map_op = |smiles: String| -> Vec<(String, String)> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let mut mol = ROMol::from_smiles(&smiles);

        trace!("calling clean");
        mol.openff_clean(); // avoids pre-condition violation on match

        trace!("calling find_matches");
        let matches = find_matches(&params, &mol);

        let mut res: Vec<(String, String)> = Vec::new();
        for pid in matches.intersection(&want) {
            res.push((pid.to_string(), smiles.clone()));
        }
        res
    };
    let results: Vec<_> =
        receiver.into_iter().par_bridge().flat_map(map_op).collect();

    th.join().unwrap().unwrap();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for (pid, mol) in results {
        res.entry(pid.to_string()).or_default().push(mol);
    }

    if let Some(dir) = output_dir {
        let dir = Path::new(&dir);
        if !dir.exists() && std::fs::create_dir_all(dir).is_err() {
            eprintln!("failed to create output dir {dir:?}");
            eprintln!("falling back to stdout");
            print_output(res);
            return;
        }
        write_output(dir, res);
    } else {
        print_output(res);
    }
}
