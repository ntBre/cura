use std::path::PathBuf;

use crate::table::Table;

pub fn parse(
    table: &mut Table,
    input: PathBuf,
    forcefield: String,
    parameter_type: String,
    target: String,
) {
    let _ = table;
    todo!("parse {input:?} {forcefield} {parameter_type} {target}");
}
