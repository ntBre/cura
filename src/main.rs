//! this module is based on chembl-search but with a focus on storing molecules
//! from multiple sources (not just ChEMBL) in a SQLite database for later
//! queries instead of reading and processing the source SDF each time.
//! currently, the stored representation is very simple: just the SMILES and
//! JSON representation of the ROMol. later on, I can decorate these entries
//! with their Morgan fingerprints or other useful information

use clap::{Parser, Subcommand};
use cura::{
    query::{query, Filter},
    serve::serve,
    store::store,
    table::Table,
};
use log::info;

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

        /// Filters to apply to the query. Currently supported values are
        /// `inchi:filename` to filter by a sequence of inchi keys stored in
        /// filename; and `natoms:n` to filter by maximum number of atoms n.
        #[arg(short = 'x', long)]
        filters: Vec<String>,

        /// Reset (clear) the stored molecule matches in the database for this
        /// force field before running the query.
        #[arg(short, long, default_value_t = false)]
        reset: bool,
    },

    Serve {
        /// The OpenFF force field to load from the toolkit.
        #[arg(short, long, default_value = "openff-2.1.0.offxml")]
        forcefield: String,
    },
}

/// Parse a sequence of filters like ["inchi:inchis.dat", "natoms:100", ...]
/// into [Filter]s of the appropriate types
fn parse_filters(filters: Vec<String>) -> Vec<Filter> {
    filters
        .into_iter()
        .map(|s| {
            let fields: Vec<&str> = s.trim().split(':').collect();
            match &fields[..] {
                [typ, arg] if *typ == "inchi" => Filter::Inchi(
                    std::fs::read_to_string(*arg)
                        .unwrap()
                        .split_ascii_whitespace()
                        .map(str::trim)
                        .map(str::to_owned)
                        .collect(),
                ),
                [typ, arg] if *typ == "natoms" => {
                    Filter::Natoms(arg.parse().unwrap())
                }
                _ => panic!("unknown filter argument {s}"),
            }
        })
        .collect()
}

#[tokio::main]
async fn main() {
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
            filters,
            reset,
        } => {
            if reset {
                info!("deleting entry for {forcefield} from database");
                table.reset_forcefield(&forcefield).unwrap();
            }
            query(
                &mut table,
                forcefield,
                parameter_type,
                search_params,
                &parse_filters(filters),
            )
        }
        Commands::Serve { forcefield } => serve(table, forcefield).await,
    }
}
