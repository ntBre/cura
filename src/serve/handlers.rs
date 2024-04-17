use std::{
    collections::HashMap,
    fs::File,
    ops::AddAssign,
    sync::{Arc, Mutex},
    time::Instant,
};

use askama::Template;
use axum::{
    extract::{Path, Query, State},
    http::StatusCode,
    response::{Html, Redirect},
    Json,
};
use clustrs::{dbscan, Label};
use log::debug;
use openff_toolkit::ForceField;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rdkit_rs::{fingerprint::tanimoto, ROMol};

use crate::{
    find_matches_full, load_forcefield, make_fps,
    serve::{
        templates::{Cluster, DrawMol, ErrorPage, Index, Param, Preview},
        AppState,
    },
    Pid,
};

pub(crate) async fn index(
    State(state): State<Arc<Mutex<AppState>>>,
) -> Html<String> {
    let state = state.lock().unwrap();
    let mut parameter_ids: Vec<_> = state
        .table
        .get_forcefield(&state.forcefield)
        .unwrap()
        .matches
        .into_iter()
        .map(|m| (m.pid, m.molecules.len()))
        .collect();
    parameter_ids.sort_by_key(|(pid, _)| {
        let mut chars = pid.chars();
        let prefix = chars.next().unwrap();
        let number: String =
            chars.by_ref().take_while(|c| c.is_numeric()).collect();
        let suffix: Vec<_> = chars.collect();
        (prefix, number.parse::<usize>().unwrap(), suffix)
    });
    let (parameter_ids, cluster_counts): (Vec<_>, Vec<_>) =
        parameter_ids.into_iter().unzip();
    let ds_size = state.table.get_dataset_size().unwrap();

    // get all dataset entries and map back to pids
    let entries: Vec<(String, Pid)> =
        state.table.get_dataset_entries().unwrap();
    let mut map = HashMap::new();
    for (_smiles, pid) in entries {
        map.entry(pid).or_insert(0).add_assign(1);
    }

    let mut pid_counts = Vec::with_capacity(parameter_ids.len());
    for pid in &parameter_ids {
        pid_counts.push(map.get(pid).copied().unwrap_or(0));
    }

    Index {
        parameter_ids,
        molecule_counts: cluster_counts,
        pid_counts,
        ds_size,
    }
    .render()
    .unwrap()
    .into()
}

/// returns the generated clustering report from [Report::write] as a String,
/// along with the number of clusters
fn make_cluster_report(
    ff: &str,
    mols: Vec<ROMol>,
    eps: f64,
    min_pts: usize,
) -> Result<Report, std::io::Error> {
    // TODO pass these in
    const MORGAN_RADIUS: u32 = 4;

    debug!("building map");

    // pid to smirks
    let map: HashMap<String, String> = ForceField::load(ff)
        .unwrap()
        .get_parameter_handler("ProperTorsions") // TODO
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), p.smirks()))
        .collect();

    // pid to mol
    let mol_map = load_forcefield(ff, "ProperTorsions"); // TODO

    debug!("loading molecules");

    let fps: Vec<_> = make_fps(&mols, MORGAN_RADIUS);
    let nfps = fps.len();
    let distance_fn = |i, j| {
        if i == j {
            0.0
        } else {
            1.0 - tanimoto(&fps[i], &fps[j])
        }
    };

    let labels = dbscan(nfps, nfps, distance_fn, eps, min_pts);

    let max = *labels
        .iter()
        .filter_map(|l| match l {
            Label::Cluster(n) => Some(n),
            _ => None,
        })
        .max()
        .unwrap_or(&0);

    let mut clusters: Vec<Vec<usize>> = vec![vec![]; max + 1];
    let mut noise = Vec::new();
    for (i, l) in labels.iter().enumerate() {
        match l {
            Label::Cluster(n) => clusters[*n].push(i),
            _ => noise.push(i),
        }
    }

    let noise_pts = noise.len();
    if clusters[0].is_empty() {
        eprintln!("warning: all noise points: {labels:?}");
        clusters[0] = noise;
    }

    Ok(Report {
        max,
        nfps,
        noise: noise_pts,
        clusters,
        mols,
        map,
        mol_map,
    })
}

/// Get the list of SMILES matching `pid` in `ffname`. Uses the cached value in
/// `state` if available, or computes and caches the result if not.
fn get_smiles_list(
    state: &mut std::sync::MutexGuard<'_, AppState>,
    pid: &String,
    ffname: &str,
) -> Vec<String> {
    let ps = state.param_states.get(pid);
    if ps.is_none() {
        let collect = state.table.get_smiles_matching(ffname, pid).unwrap();
        state
            .param_states
            .entry(pid.clone())
            .or_default()
            .smiles_list = Some(collect);
    }
    let ps = state.param_states.get_mut(pid).unwrap();
    ps.smiles_list.clone().unwrap()
}

fn get_param<T: std::str::FromStr>(
    params: &HashMap<String, String>,
    key: &str,
    def: T,
) -> T {
    match params.get(key).map(|s| s.parse()) {
        Some(Ok(t)) => t,
        _ => def,
    }
}

struct Report {
    max: usize,
    nfps: usize,
    noise: usize,
    clusters: Vec<Vec<usize>>,
    mols: Vec<ROMol>,
    map: HashMap<String, String>,
    mol_map: Vec<(String, ROMol)>,
}

pub(crate) async fn cluster(
    State(state): State<Arc<Mutex<AppState>>>,
    Path(pid): Path<String>,
    Query(params): Query<HashMap<String, String>>,
) -> Html<String> {
    let mut state = state.lock().unwrap();
    let ffname = state.forcefield.clone();
    let Some(smarts) = state.pid_to_smarts.get(&pid).cloned() else {
        return ErrorPage { pid }.render().unwrap().into();
    };

    let mols: Vec<_> = get_smiles_list(&mut state, &pid, &ffname)
        .into_par_iter()
        .map(|s| {
            let mut mol = ROMol::from_smiles(&s);
            mol.openff_clean();
            mol
        })
        .collect();

    const DBSCAN_EPS: f64 = 0.5;
    const DBSCAN_MIN_PTS: usize = 1;

    let eps = get_param(&params, "eps", DBSCAN_EPS);
    let min_pts = get_param(&params, "min_pts", DBSCAN_MIN_PTS);

    let start = Instant::now();

    let Report {
        max,
        nfps,
        noise,
        mut clusters,
        mols,
        map,
        mol_map,
    } = make_cluster_report(&ffname, mols, eps, min_pts).unwrap();

    clusters.sort_by_key(|c| mols[c[0]].num_atoms());

    Cluster {
        pid,
        smarts,
        eps,
        min_pts,
        max,
        nfps,
        noise,
        clusters,
        mols,
        map,
        mol_map,
        time: start.elapsed().as_millis() as f64 / 1000.0,
    }
    .render()
    .unwrap()
    .into()
}

/// return the index of the smallest molecule in cluster
pub(crate) fn find_smallest(mols: &[ROMol], cluster: &[usize]) -> usize {
    cluster
        .iter()
        .enumerate()
        .map(|(i, c)| (i, &mols[*c]))
        .min_by_key(|(_, mol)| mol.num_atoms())
        .unwrap()
        .0
}

pub(crate) async fn param(
    State(state): State<Arc<Mutex<AppState>>>,
    Path(pid): Path<String>,
    Query(params): Query<HashMap<String, String>>,
) -> Html<String> {
    let mut state = state.lock().unwrap();
    let ffname = state.forcefield.clone();
    let Some(smarts) = state.pid_to_smarts.get(&pid).cloned() else {
        return ErrorPage { pid }.render().unwrap().into();
    };

    const MAX_DRAW: usize = 50; /* the maximum number of mols to draw */
    let max_draw = match params.get("max").map(|s| s.parse()) {
        Some(Ok(n)) => n,
        _ => MAX_DRAW,
    };

    let mut mols: Vec<_> = get_smiles_list(&mut state, &pid, &ffname)
        .into_par_iter()
        .map(|s| {
            let mut mol = ROMol::from_smiles(&s);
            mol.openff_clean();
            let natoms = mol.num_atoms();
            (mol, s, natoms)
        })
        .collect();
    mols.sort_by_key(|(_mol, _smiles, natoms)| *natoms);

    let mol_map = load_forcefield(&state.forcefield, "ProperTorsions");
    let total_mols = mols.len();
    let mols = mols
        .into_iter()
        .take(max_draw)
        .map(|(mut mol, smiles, _natoms)| {
            mol.openff_clean();
            let (hl_atoms, _pid) = find_matches_full(&mol_map, &mol)
                .into_iter()
                .find(|(_atoms, param_id)| param_id == &pid)
                .unwrap_or_default();
            let svg = mol.draw_svg(300, 300, "", &hl_atoms);
            let natoms = mol.num_atoms();
            DrawMol {
                smiles,
                natoms,
                svg,
            }
        })
        .collect();
    let tmpl = Param {
        smarts,
        pid: pid.clone(),
        mols,
        total_mols,
    };
    let slot = state.param_states.get_mut(&pid).unwrap();
    slot.param_page = Some(tmpl.clone());
    tmpl.render().unwrap().into()
}

#[derive(serde::Deserialize)]
pub(crate) struct AddMolecule {
    smiles: String,
    pid: Pid,
}

pub(crate) async fn add_molecule(
    State(state): State<Arc<Mutex<AppState>>>,
    Json(body): Json<AddMolecule>,
) -> StatusCode {
    let state = state.lock().unwrap();
    let AddMolecule { smiles, pid } = body;
    debug!("adding smiles: `{smiles}` to database for pid: `{pid}`");
    match state.table.add_to_dataset(smiles, pid) {
        Ok(_) => StatusCode::OK,
        Err(e) => {
            debug!("error adding to dataset: {e:?}");
            return StatusCode::BAD_REQUEST;
        }
    }
}

pub(crate) async fn reset_dataset(
    State(state): State<Arc<Mutex<AppState>>>,
) -> Redirect {
    let state = state.lock().unwrap();
    state.table.reset_dataset().unwrap();
    Redirect::to("/")
}

pub(crate) async fn preview_dataset(
    State(state): State<Arc<Mutex<AppState>>>,
) -> Html<String> {
    let state = state.lock().unwrap();
    let mols = state
        .table
        .get_dataset_entries()
        .unwrap()
        .into_iter()
        .map(|(s, _pid)| {
            let mut mol = ROMol::from_smiles(&s);
            mol.openff_clean();
            DrawMol {
                smiles: s,
                natoms: mol.num_atoms(),
                svg: mol.draw_svg(300, 300, "", &[]),
            }
        })
        .collect();
    Preview { mols }.render().unwrap().into()
}

pub(crate) async fn export_dataset(
    State(state): State<Arc<Mutex<AppState>>>,
    body: String,
) -> Redirect {
    // expecting a request body like `filename="dataset.smi"`
    let &[label, filename] = &body.split('=').collect::<Vec<_>>()[..] else {
        panic!("filename not in body: {body:?}");
    };
    assert_eq!(label, "filename");
    let state = state.lock().unwrap();
    let mut out = File::create(filename).unwrap();
    use std::io::Write;
    for (smiles, pid) in state.table.get_dataset_entries().unwrap() {
        writeln!(out, "{pid} {smiles}").unwrap();
    }
    Redirect::to("/")
}
