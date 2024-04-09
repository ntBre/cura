use std::collections::HashMap;
use std::ops::Index as _;

use askama::Template;
use rdkit_rs::ROMol;

#[derive(Template)]
#[template(path = "index.html")]
pub(crate) struct Index {
    pub(crate) parameter_ids: Vec<String>,
    pub(crate) molecule_counts: Vec<usize>,
}

#[derive(Clone)]
pub struct DrawMol {
    pub smiles: String,
    pub natoms: usize,
    pub svg: String,
}

#[derive(Clone, Template)]
#[template(path = "param.html")]
pub(crate) struct Param {
    pub smarts: String,
    pub pid: String,
    pub total_mols: usize,
    pub mols: Vec<DrawMol>,
}

#[derive(Clone, Template)]
#[template(path = "cluster.html")]
pub(crate) struct Cluster {
    pub smarts: String,
    pub pid: String,
    pub eps: f64,
    pub min_pts: usize,
    pub max: usize,
    pub nfps: usize,
    pub noise: usize,
    pub clusters: Vec<Vec<usize>>,
    pub mols: Vec<ROMol>,
    pub map: HashMap<String, String>,
    pub mol_map: Vec<(String, ROMol)>,

    /// The duration of the clustering process
    pub time: f64,
}

impl Cluster {
    pub fn make_svg(&self, mol: &ROMol) -> String {
        let mut hl_atoms = Vec::new();
        let pid = &self.pid;
        if self.map.contains_key(pid) {
            let tmp = crate::find_matches_full(&self.mol_map, mol);
            let got = tmp.iter().find(|(_atoms, param_id)| param_id == &pid);
            if let Some((atoms, _pid)) = got {
                hl_atoms.clone_from(atoms);
            } else {
                panic!("smirks doesn't match any more");
            }
        }
        mol.draw_svg(400, 300, "", &hl_atoms)
    }
}

#[derive(Template)]
#[template(path = "error.html")]
pub(crate) struct ErrorPage {
    pub(crate) pid: String,
}
