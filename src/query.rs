use std::{
    collections::{HashMap, HashSet},
    path::Path,
    sync::{atomic::AtomicUsize, mpsc},
    thread,
};

use log::{info, trace};
use openff_toolkit::ForceField as OFF;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::ROMol;

use crate::{
    find_matches, load_forcefield, table::Table, ForceField, Match, Molecule,
    Pid, Smirks, PROGRESS_INTERVAL,
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

pub enum Filter {
    Elements(i128),
    Inchi(HashSet<String>),
    Natoms(usize),
}

impl Filter {
    /// Returns true if `mol` satisfies the condition of the filter and should
    /// be retained.
    fn apply(&self, mol: &Molecule) -> bool {
        match self {
            Filter::Inchi(set) => {
                !mol.inchikey.as_ref().is_some_and(|key| set.contains(key))
            }
            Filter::Natoms(n) => mol.natoms <= *n,
            Filter::Elements(mask) => (mol.elements | mask) == *mask,
        }
    }
}

pub fn query(
    table: &mut Table,
    forcefield: String,
    parameter_type: String,
    search_params: String,
    filters: &[Filter],
) {
    info!("loading force field and parameters");

    // TODO sad to load this twice but I want to store the smirks in the db
    let pid_to_smirks: HashMap<Pid, Smirks> = OFF::load(&forcefield)
        .unwrap()
        .get_parameter_handler(&parameter_type)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), p.smirks()))
        .collect();

    let params = load_forcefield(&forcefield, &parameter_type);

    let want = load_want(&search_params);
    info!("loading moldata from database");
    const CHAN_SIZE: usize = 1024;
    let (sender, receiver) = mpsc::sync_channel(CHAN_SIZE);
    let progress = AtomicUsize::new(0);
    let skipped = AtomicUsize::new(0);
    let results: Vec<_> = thread::scope(|s| {
        s.spawn(|| table.send_molecules(sender));

        info!("processing data");
        let map_op = |mol: Molecule| -> Vec<(Pid, usize)> {
            let cur =
                progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if cur % PROGRESS_INTERVAL == 0 {
                eprint!("{cur} complete\r");
            }

            if !filters.iter().all(|f| f.apply(&mol)) {
                skipped.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                return Vec::new();
            }

            let mut romol = ROMol::from_smiles(&mol.smiles);

            trace!("calling clean");
            romol.openff_clean(); // avoids pre-condition violation on match

            trace!("calling find_matches");
            let matches = find_matches(&params, &romol);

            matches
                .intersection(&want)
                .map(|pid| (pid.to_owned(), mol.id.unwrap()))
                .collect()
        };

        receiver.into_iter().par_bridge().flat_map(map_op).collect()
    });
    let progress = progress.into_inner();
    let skipped = skipped.into_inner();

    eprintln!("filtered {skipped}/{progress} records");

    let mut res: HashMap<Pid, Match> = HashMap::new();
    for (pid, mol_id) in results {
        res.entry(pid.to_string())
            .or_insert(Match {
                smirks: pid_to_smirks[&pid].to_owned(),
                pid,
                molecules: Vec::new(),
            })
            .molecules
            .push(mol_id);
    }

    table
        .insert_forcefield(ForceField {
            id: None,
            name: forcefield,
            matches: res.into_values().collect(),
        })
        .unwrap();
}
