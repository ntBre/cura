use askama::Template;

#[derive(Template)]
#[template(path = "index.html")]
pub(crate) struct Index {
    pub(crate) parameter_ids: Vec<String>,
    pub(crate) cluster_counts: Vec<usize>,
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
    pub body: String,
}

#[derive(Template)]
#[template(path = "error.html")]
pub(crate) struct ErrorPage {
    pub(crate) pid: String,
}
