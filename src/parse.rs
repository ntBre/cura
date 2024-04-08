use std::{collections::HashMap, path::PathBuf};

use clustrs::{dbscan, Label};
use log::info;
use openff_toolkit::ForceField;
use rdkit_rs::{fingerprint::tanimoto, ROMol};

use crate::{load_mols, make_fps, table::Table, Pid, Report, Smirks};

pub fn parse(
    _table: &mut Table, // TODO use this
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
