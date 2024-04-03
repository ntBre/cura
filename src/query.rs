use std::{
    collections::{HashMap, HashSet},
    path::Path,
    sync::{atomic::AtomicUsize, mpsc},
    thread,
};

use log::{info, trace};
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::ROMol;

use crate::{
    find_matches, load_forcefield, table::Table, Molecule, PROGRESS_INTERVAL,
};

/// load a sequence of newline-separated entries from `path` and collect them
/// into a HashSet
pub fn load_want(path: &str) -> HashSet<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to open `{path}` with `{e:?}`"))
        .lines()
        .map(|s| s.trim().to_owned())
        .collect()
}

pub fn print_output(res: HashMap<String, Vec<String>>) {
    for (pid, moles) in res {
        for mol in moles {
            println!("{pid}\t{mol}");
        }
    }
}

pub fn write_output(dir: impl AsRef<Path>, res: HashMap<String, Vec<String>>) {
    use std::io::Write;
    for (pid, moles) in res {
        let path = dir.as_ref().join(pid).with_extension("smiles");
        let mut f = std::fs::File::create(path).unwrap();
        for mol in moles {
            writeln!(f, "{mol}").unwrap();
        }
    }
}

pub fn query(
    table: &mut Table,
    forcefield: String,
    parameter_type: String,
    search_params: String,
    output_dir: Option<String>,
) {
    info!("loading moldata from database");
    const CHAN_SIZE: usize = 1024;
    let (sender, receiver) = mpsc::sync_channel(CHAN_SIZE);
    let results: Vec<_> = thread::scope(|s| {
        s.spawn(|| table.send_molecules(sender));

        info!("loading force field and parameters");
        let params = load_forcefield(forcefield, parameter_type);

        let want = load_want(&search_params);

        info!("processing data");
        let progress = AtomicUsize::new(0);
        let map_op = |mol: Molecule| -> Vec<(String, String)> {
            let cur =
                progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if cur % PROGRESS_INTERVAL == 0 {
                eprint!("{cur} complete\r");
            }
            let mut romol = ROMol::from_smiles(&mol.smiles);

            trace!("calling clean");
            romol.openff_clean(); // avoids pre-condition violation on match

            trace!("calling find_matches");
            let matches = find_matches(&params, &romol);

            let mut res: Vec<(String, String)> = Vec::new();
            for pid in matches.intersection(&want) {
                res.push((pid.to_string(), mol.smiles.clone()));
            }
            res
        };

        receiver.into_iter().par_bridge().flat_map(map_op).collect()
    });

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
