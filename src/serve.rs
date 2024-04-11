use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use axum::{
    routing::{get, post},
    Router,
};
use openff_toolkit::ForceField;
use templates::Param;
use tokio::net::TcpListener;

use crate::table::Table;

mod handlers;
mod templates;

#[derive(Default)]
struct ParamState {
    smiles_list: Option<Vec<String>>,
    param_page: Option<Param>,
}

struct AppState {
    forcefield: String,
    pid_to_smarts: HashMap<String, String>,
    param_states: HashMap<String, ParamState>,
    table: Table,
}

pub async fn serve(table: Table, forcefield: String) {
    let pid_to_smarts: HashMap<String, String> = ForceField::load(&forcefield)
        .unwrap()
        .get_parameter_handler("ProperTorsions")
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), p.smirks()))
        .collect();

    let state = Arc::new(Mutex::new(AppState {
        forcefield,
        param_states: HashMap::new(),
        pid_to_smarts,
        table,
    }));

    let app = Router::new()
        .route("/", get(handlers::index))
        .route("/param/:pid", get(handlers::param))
        .route("/cluster/:pid", get(handlers::cluster))
        .route("/add-molecule", post(handlers::add_molecule))
        .route("/reset-dataset", post(handlers::reset_dataset))
        .route("/preview-dataset", get(handlers::preview_dataset))
        .route("/export-dataset", get(handlers::export_dataset))
        .nest_service("/css", tower_http::services::ServeDir::new("css"))
        .nest_service("/js", tower_http::services::ServeDir::new("js"))
        .with_state(state);

    let addr = "0.0.0.0:3000";
    let listener = TcpListener::bind(addr).await.unwrap();
    eprintln!("serving on {addr}");
    axum::serve(listener, app).await.unwrap();
}
