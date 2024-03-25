//! this module is based on chembl-search but with a focus on storing molecules
//! from multiple sources (not just ChEMBL) in a SQLite database for later
//! queries instead of reading and processing the source SDF each time.
//! currently, the stored representation is very simple: just the SMILES and
//! JSON representation of the ROMol. later on, I can decorate these entries
//! with their Morgan fingerprints or other useful information

use std::{collections::HashMap, path::Path, sync::atomic::AtomicUsize};

use clap::{Parser, Subcommand};
use cura::{Table, PROGRESS_INTERVAL};
use log::{info, trace, warn};
use openff_toolkit::ForceField;
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};
use rdkit_rs::{RDError, ROMol, SDMolSupplier};
use rsearch::{find_matches, load_want, print_output, write_output};

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

    /// Query existing molecules in the database for matches to OpenFF
    /// parameters in a force field
    Query {
        /// The OpenFF force field to load from the toolkit.
        #[arg(short, long, default_value = "openff-2.1.0.offxml")]
        forcefield: String,

        /// The `Parameter` type for which to extract parameters. Allowed
        /// options are valid arguments to `ForceField.get_parameter_handler`,
        /// such as Bonds, Angles, or ProperTorsions.
        #[arg(short, long, default_value = "ProperTorsions")]
        parameter_type: String,

        /// The path to the file listing the parameter identifiers to match
        /// against, one per line. These must correspond to parameters in the
        /// provided forcefield.
        #[arg(short, long, default_value = "want.params")]
        search_params: String,

        /// Where to write output SMILES files, one for each input parameter. If
        /// false, print the output to stdout.
        #[arg(short, long)]
        output_dir: Option<String>,
    },
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

fn query(
    table: &mut Table,
    forcefield: String,
    parameter_type: String,
    search_params: String,
    output_dir: Option<String>,
) {
    info!("loading moldata from database");
    let data = table.get_moldata().unwrap();

    info!("loading force field and parameters");
    let ff = ForceField::load(&forcefield).unwrap();
    let h = ff.get_parameter_handler(&parameter_type).unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), ROMol::from_smarts(&p.smirks())));
    }

    let want = load_want(&search_params);

    info!("processing data");
    let progress = AtomicUsize::new(0);
    let map_op = |mol: &String| -> Vec<(String, String)> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % PROGRESS_INTERVAL == 0 {
            eprint!("{cur} complete\r");
        }
        let mut mol = ROMol::from_json(mol);

        trace!("calling clean");
        mol.openff_clean(); // avoids pre-condition violation on match

        trace!("calling find_matches");
        let matches = find_matches(&params, &mol);

        let mut res: Vec<(String, String)> = Vec::new();
        let mut smiles = None;
        for pid in matches.intersection(&want) {
            if smiles.is_none() {
                smiles = Some(mol.to_smiles());
            }
            res.push((pid.to_string(), smiles.clone().unwrap()));
        }
        res
    };
    let results: Vec<_> = data.par_iter().flat_map(map_op).collect();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for (pid, mol) in results {
        res.entry(pid.to_string()).or_default().push(mol);
    }

    if let Some(dir) = output_dir {
        let dir = Path::new(&dir);
        if !dir.exists() && std::fs::create_dir_all(dir).is_err() {
            eprintln!("failed to create output dir {dir:?}");
            eprintln!("falling back to stdout");
            print_output(res);
            return;
        }
        write_output(dir, res);
    } else {
        print_output(res);
    }
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
        Commands::Query {
            forcefield,
            parameter_type,
            search_params,
            output_dir,
        } => query(
            &mut table,
            forcefield,
            parameter_type,
            search_params,
            output_dir,
        ),
    }
}
