use std::fs::read_to_string;
use std::path::Path;

use serde::Deserialize;

#[derive(Clone, Deserialize)]
pub(crate) struct Dbscan {
    /// The maximum acceptable distance between a core point of a cluster and
    /// one of its neighbors.
    #[serde(default = "Dbscan::default_eps")]
    pub(crate) epsilon: f64,

    /// The minimum number of points required to form a dense region
    #[serde(default = "Dbscan::default_min_pts")]
    pub(crate) min_pts: usize,
}

impl Dbscan {
    fn default_eps() -> f64 {
        0.5
    }

    fn default_min_pts() -> usize {
        1
    }
}

#[derive(Deserialize)]
pub(crate) struct Config {
    /// The maximum number of atoms to consider
    pub(crate) max_atoms: usize,

    /// The force field to use for parameter labeling.
    pub(crate) forcefield: String,

    /// Morgan fingerprinting radius
    pub(crate) radius: u32,

    #[serde(rename = "parameter")]
    pub(crate) parameters: Vec<Parameter>,
}

#[derive(Clone, Deserialize)]
pub(crate) struct Parameter {
    /// The file of SMILES strings to read as input, one SMILES per line.
    pub(crate) smiles: String,

    /// The parameter identifier to use when highlighting atoms in the
    /// molecules.
    pub(crate) id: String,

    /// The `Parameter` type for which to extract parameters. Allowed options
    /// are valid arguments to `ForceField.get_parameter_handler`, such as
    /// Bonds, Angles, or ProperTorsions.
    #[serde(rename = "type", default = "default_ptype")]
    pub(crate) typ: String,

    /// [DBSCAN] parameters
    #[serde(default = "default_dbscan")]
    pub(crate) dbscan: Dbscan,

    /// Whether or not to fragment the molecules before the fingerprinting
    /// analysis.
    #[serde(default = "default_fragment")]
    pub(crate) fragment: bool,
}

fn default_fragment() -> bool {
    true
}

fn default_ptype() -> String {
    "ProperTorsions".to_owned()
}

fn default_dbscan() -> Dbscan {
    Dbscan {
        epsilon: 0.5,
        min_pts: 1,
    }
}

impl Config {
    pub(crate) fn load(path: impl AsRef<Path>) -> Self {
        toml::from_str(&read_to_string(path).unwrap()).unwrap()
    }
}
