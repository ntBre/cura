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

#[derive(Clone)]
pub enum Body {
    SmilesList {
        total_mols: usize,
        mols: Vec<DrawMol>,
    },
    #[allow(unused)]
    Report(String),
}

#[derive(Clone, Template)]
#[template(path = "param.html")]
pub(crate) struct Param {
    pub smarts: String,
    pub pid: String,
    pub body: Body,
}

#[derive(Template)]
#[template(path = "error.html")]
pub(crate) struct ErrorPage {
    pub(crate) pid: String,
}
