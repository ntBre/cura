//! this module is based on chembl-search but with a focus on storing molecules
//! from multiple sources (not just ChEMBL) in a SQLite database for later
//! queries instead of reading and processing the source SDF each time.
//! currently, the stored representation is very simple: just the SMILES and
//! JSON representation of the ROMol. later on, I can decorate these entries
//! with their Morgan fingerprints or other useful information

use std::sync::atomic::AtomicUsize;

use clap::{Parser, Subcommand};
use cura::{Table, PROGRESS_INTERVAL};
use log::{trace, warn};
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::{RDError, ROMol, SDMolSupplier};

#[derive(Parser)]
struct Cli {
    #[arg(short, long, default_value = "try.sqlite")]
    database: String,

    /// The number of threads to use. Defaults to the number of logical CPUs as
    /// detected by rayon.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Store molecules from the SDF into the database
    Store { molecule_file: String },
}

fn store(table: &mut Table, molecule_file: String) {
    trace!("initializing mol supplier from {}", molecule_file);
    let m = SDMolSupplier::new(molecule_file).unwrap();

    let progress = AtomicUsize::new(0);

    let map_op = |mol: Result<ROMol, RDError>| -> Option<(String, String)> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let Ok(mut mol) = mol else {
            warn!("error loading molecule, skipping");
            return None;
        };
        mol.openff_clean();

        Some((mol.to_smiles(), mol.to_json()))
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    table.insert_molecules(results).unwrap();
}

fn main() {
    env_logger::init();

    let cli = Cli::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let mut table = Table::create(&cli.database).unwrap();

    match cli.command {
        Commands::Store { molecule_file } => store(&mut table, molecule_file),
    }
}
