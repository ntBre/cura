//! this module is based on chembl-search but with a focus on storing molecules
//! from multiple sources (not just ChEMBL) in a SQLite database for later
//! queries instead of reading and processing the source SDF each time.
//! currently, the stored representation is very simple: just the SMILES and
//! JSON representation of the ROMol. later on, I can decorate these entries
//! with their Morgan fingerprints or other useful information

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use cura::{parse::parse, query::query, store::store, table::Table};

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

    /// Parse the output from `cura query`, then fragment, fingerprint, and
    /// cluster the results. Writes the output to `input` with the extension
    /// `html`
    Parse {
        /// The OpenFF force field to load from the toolkit.
        #[arg(short, long, default_value = "openff-2.1.0.offxml")]
        forcefield: String,

        /// The `Parameter` type for which to extract parameters. Allowed
        /// options are valid arguments to `ForceField.get_parameter_handler`,
        /// such as Bonds, Angles, or ProperTorsions.
        #[arg(short, long, default_value = "ProperTorsions")]
        parameter_type: String,

        /// The target parameter. Must correspond to a parameter of type
        /// `parameter_type` in `forcefield`.
        #[arg(short, long)]
        target: String,

        /// The input SMILES file to read. TODO handle a directory recursively
        /// or at least multiple files
        input: PathBuf,
    },
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
        Commands::Store { molecule_file } => store(&mut table, &molecule_file),
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
        Commands::Parse {
            input,
            forcefield,
            parameter_type,
            target,
        } => parse(&mut table, input, forcefield, parameter_type, target),
    }
}
