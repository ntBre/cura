use openff_toolkit::ForceField;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rdkit_rs::{fingerprint::tanimoto, ROMol};
use rsearch::{
    cluster::{dbscan, Label},
    find_matches_full,
    utils::{fragment, make_fps, Report},
};
use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
    sync::atomic::AtomicUsize,
    sync::atomic::Ordering,
};

use log::{info, trace};

use crate::table::Table;

type Pid = String;
type Smirks = String;

pub fn load_mols(
    smiles: Vec<&str>,
    max_atoms: usize,
    do_fragment: bool,
    pid: &str,
    mol_map: &[(String, ROMol)],
) -> Vec<ROMol> {
    // this gives each of the "fragments" from the original smiles
    let mut mols: Vec<_> = smiles.into_iter().map(ROMol::from_smiles).collect();

    info!("collected {} initial molecules", mols.len());

    if do_fragment {
        mols = fragment(mols);
    }

    let too_big = AtomicUsize::new(0);
    let no_match = AtomicUsize::new(0);
    let overlap = AtomicUsize::new(0);
    let fragments = AtomicUsize::new(0);

    let mut ret: Vec<_> = mols
        .into_par_iter()
        .flat_map(|mut mol| {
            mol.openff_clean();
            fragments.fetch_add(1, Ordering::Relaxed);
            if mol.num_atoms() > max_atoms {
                too_big.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            let matches = find_matches_full(mol_map, &mol);
            if !matches.iter().any(|(_, p)| p.as_str() == pid) {
                trace!(
                    "pid {pid} not found in {:?} for {}",
                    matches.into_values().collect::<HashSet<String>>(),
                    mol.to_smiles(),
                );
                no_match.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            Some((mol.to_smiles(), mol))
        })
        .collect();

    info!(
        "expanded initial smiles to {} fragments",
        fragments.into_inner()
    );

    let presort = ret.len();

    ret.sort_by_key(|(smiles, _mol)| smiles.clone());
    ret.dedup_by_key(|(smiles, _mol)| smiles.clone());

    info!(
        "filtered {} for size, {} for smirks, {} for inchi, {} for duplicates",
        too_big.into_inner(),
        no_match.into_inner(),
        overlap.into_inner(),
        presort - ret.len(),
    );

    let (_, ret): (Vec<String>, _) = ret.into_iter().unzip();

    ret
}

pub fn parse(
    table: &mut Table,
    input: PathBuf,
    forcefield: String,
    parameter_type: String,
    target: String,
) {
    let s = std::fs::read_to_string(&input)
        .unwrap_or_else(|e| panic!("failed to read {:?} for {e}", input));
    let smiles: Vec<_> = s.lines().collect();

    info!("{} initial smiles", smiles.len());

    let parameters = ForceField::load(&forcefield)
        .unwrap()
        .get_parameter_handler(&parameter_type)
        .unwrap()
        .parameters();

    let map: HashMap<Pid, Smirks> =
        parameters.iter().map(|p| (p.id(), p.smirks())).collect();

    let mol_map: Vec<(Pid, ROMol)> = parameters
        .into_iter()
        .map(|p| (p.id(), ROMol::from_smarts(&p.smirks())))
        .collect();

    const MAX_ATOMS: usize = 80;
    const DO_FRAGMENT: bool = false;
    const MORGAN_RADIUS: u32 = 4;

    let mols = load_mols(smiles, MAX_ATOMS, DO_FRAGMENT, &target, &mol_map);

    info!("making fingerprints");

    let fps: Vec<_> = make_fps(&mols, MORGAN_RADIUS);
    let nfps = fps.len();
    let distance_fn = |i, j| {
        if i == j {
            0.0
        } else {
            1.0 - tanimoto(&fps[i], &fps[j])
        }
    };

    info!("running DBSCAN");

    const DBSCAN_EPS: f64 = 0.5;
    const DBSCAN_MIN_PTS: usize = 1;

    let labels = dbscan(nfps, nfps, distance_fn, DBSCAN_EPS, DBSCAN_MIN_PTS);

    let max = match labels
        .iter()
        .filter_map(|l| match l {
            Label::Cluster(n) => Some(n),
            _ => None,
        })
        .max()
    {
        Some(n) => *n,
        None => {
            dbg!(labels);
            eprintln!("error: all noise points, exiting");
            std::process::exit(1);
        }
    };

    let mut clusters: Vec<Vec<usize>> = vec![vec![]; max + 1];
    let mut noise = 0;
    for (i, l) in labels.iter().enumerate() {
        match l {
            Label::Cluster(n) => clusters[*n].push(i),
            _ => noise += 1,
        }
    }

    info!("generating report");

    let output = input.with_extension("html");
    Report {
        args: std::env::args().collect::<Vec<_>>(),
        max,
        nfps,
        noise,
        clusters,
        mols,
        parameter: &target,
        map,
        mol_map,
    }
    .generate(output)
    .unwrap();
}
