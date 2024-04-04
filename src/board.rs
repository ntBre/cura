use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use axum::{routing::get, Router};
use config::{Config, Parameter};
use openff_toolkit::ForceField;
use templates::Param;
use tokio::net::TcpListener;

use crate::table::Table;

mod config;
mod handlers;
mod templates;

#[derive(Default)]
struct ParamState {
    smiles_list: Option<Vec<String>>,
    param_page: Option<Param>,
    nclusters: usize,
}

struct AppState {
    cli: Config,
    pid_to_smarts: HashMap<String, String>,
    param_states: HashMap<String, ParamState>,
    table: Table,
}

impl AppState {
    #[inline]
    fn param_by_id(&self, id: &str) -> Option<&Parameter> {
        self.cli.parameters.iter().find(|p| p.id == id)
    }

    #[inline]
    fn param_by_id_mut(&mut self, id: &str) -> Option<&mut Parameter> {
        self.cli.parameters.iter_mut().find(|p| p.id == id)
    }
}

pub async fn board(table: Table) {
    let cli = config::Config::load("testfiles/fingerprint.toml");

    let pid_to_smarts: HashMap<String, String> =
        ForceField::load(&cli.forcefield)
            .unwrap()
            .get_parameter_handler("ProperTorsions")
            .unwrap()
            .parameters()
            .into_iter()
            .map(|p| (p.id(), p.smirks()))
            .collect();

    let state = Arc::new(Mutex::new(AppState {
        cli,
        param_states: HashMap::new(),
        pid_to_smarts,
        table,
    }));

    let app = Router::new()
        .route("/", get(handlers::index))
        .route("/param/:pid", get(handlers::param))
        .nest_service("/css", tower_http::services::ServeDir::new("css"))
        .with_state(state);
    let addr = "0.0.0.0:3000";
    let listener = TcpListener::bind(addr).await.unwrap();
    eprintln!("serving on {addr}");
    axum::serve(listener, app).await.unwrap();
}
