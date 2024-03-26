use std::{collections::HashMap, path::Path, sync::atomic::AtomicUsize};

use log::{info, trace};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rdkit_rs::ROMol;
use rsearch::{find_matches, load_want, print_output, write_output};

use crate::{load_forcefield, table::Table, PROGRESS_INTERVAL};

pub fn query(
    table: &mut Table,
    forcefield: String,
    parameter_type: String,
    search_params: String,
    output_dir: Option<String>,
) {
    info!("loading moldata from database");
    let data = table.get_moldata().unwrap();

    info!("loading force field and parameters");
    let params = load_forcefield(forcefield, parameter_type);

    let want = load_want(&search_params);

    info!("processing data");
    let progress = AtomicUsize::new(0);
    let map_op = |mol: &String| -> Vec<(String, String)> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let mut mol = ROMol::from_json(mol);

        trace!("calling clean");
        mol.openff_clean(); // avoids pre-condition violation on match

        trace!("calling find_matches");
        let matches = find_matches(&params, &mol);

        let mut res: Vec<(String, String)> = Vec::new();
        let mut smiles = None;
        for pid in matches.intersection(&want) {
            if smiles.is_none() {
                smiles = Some(mol.to_smiles());
            }
            res.push((pid.to_string(), smiles.clone().unwrap()));
        }
        res
    };
    let results: Vec<_> = data.par_iter().flat_map(map_op).collect();

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
