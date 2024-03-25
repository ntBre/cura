use cura::Table;
use log::info;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rdkit_rs::ROMol;

fn main() {
    env_logger::init();
    let table = Table::open("try.sqlite").unwrap();
    info!("loading moldata from database");
    let data = table.get_moldata().unwrap();
    info!("processing data");
    data.par_iter().for_each(|moldata| {
        let mol = ROMol::from_json(&moldata);
        mol.to_inchi_key();
    });
}
