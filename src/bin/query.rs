use std::{collections::HashMap, sync::atomic::AtomicUsize};

use cura::{Table, PROGRESS_INTERVAL};
use log::{info, trace};
use openff_toolkit::ForceField;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rdkit_rs::ROMol;
use rsearch::{find_matches, load_want, write_output};

fn main() {
    env_logger::init();
    let table = Table::open("try.sqlite").unwrap();
    info!("loading moldata from database");
    let data = table.get_moldata().unwrap();

    info!("loading force field and parameters");
    let ff = ForceField::load("openff-2.1.0.offxml").unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), ROMol::from_smarts(&p.smirks())));
    }

    let want = load_want("/home/brent/omsf/chembl/input/want.params");

    info!("processing data");
    let progress = AtomicUsize::new(0);
    let map_op = |mol: &String| -> Vec<(String, String)> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let mut mol = ROMol::from_json(&mol);
        trace!("calling clean");

        // necessary to avoid pre-condition violation on match
        mol.openff_clean();

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

    write_output("output", res);
}
