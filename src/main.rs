//! this module is based on chembl-search but with a focus on storing molecules
//! from multiple sources (not just ChEMBL) in a SQLite database for later
//! queries instead of reading and processing the source SDF each time.
//! currently, the stored representation is very simple: just the SMILES and
//! JSON representation of the ROMol.

use std::path::Path;
use std::sync::atomic::AtomicUsize;

use clap::Parser;
use log::trace;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::{RDError, ROMol, SDMolSupplier};
use rusqlite::Connection;

#[derive(Parser)]
struct Cli {
    /// The OpenFF force field to load from the toolkit.
    #[arg(short, long, default_value = "openff-2.1.0.offxml")]
    forcefield: String,

    /// The `Parameter` type for which to extract parameters. Allowed options
    /// are valid arguments to `ForceField.get_parameter_handler`, such as
    /// Bonds, Angles, or ProperTorsions.
    #[arg(short, long, default_value = "ProperTorsions")]
    parameter_type: String,

    /// The path to the SDF file from which to read Molecules.
    #[arg(short, long, default_value = "chembl_33.sdf")]
    molecule_file: String,

    /// The path to the file listing the parameter identifiers to match against,
    /// one per line. These must correspond to parameters in the provided
    /// forcefield.
    #[arg(short, long, default_value = "want.params")]
    search_params: String,

    /// The number of threads to use. Defaults to the number of logical CPUs as
    /// detected by rayon.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Where to write output SMILES files, one for each input parameter. If
    /// false, print the output to stdout.
    #[arg(short, long)]
    output_dir: Option<String>,
}

struct Table {
    conn: Connection,
}

impl Table {
    fn new(path: impl AsRef<Path>) -> rusqlite::Result<Self> {
        let conn = Connection::open(path)?;
        conn.execute(include_str!("create_table.sql"), ())?;
        Ok(Self { conn })
    }

    fn insert_molecule(
        &self,
        smiles: String,
        moldata: String,
    ) -> rusqlite::Result<()> {
        let mut stmt =
            self.conn.prepare(include_str!("insert_molecule.sql"))?;
        stmt.execute((smiles, moldata))?;
        Ok(())
    }
}

fn main() {
    env_logger::init();

    let table = Table::new("try.sqlite").unwrap();

    let cli = Cli::parse();
    trace!("initializing mol supplier from {}", cli.molecule_file);
    let m = SDMolSupplier::new(cli.molecule_file).unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let progress = AtomicUsize::new(0);

    const PROGRESS_INTERVAL: usize = 1000;

    let map_op = |mol: Result<ROMol, RDError>| -> Option<(String, String)> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let Ok(mut mol) = mol else {
            return None;
        };
        trace!("calling clean");
        mol.openff_clean();

        Some((mol.to_smiles(), mol.to_json()))
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    for (smiles, moldata) in results {
        table.insert_molecule(smiles, moldata).unwrap();
    }
}
