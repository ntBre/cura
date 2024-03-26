use std::io::Read;
use std::io::Write;
use std::path::Path;

use flate2::{read::ZlibDecoder, write::ZlibEncoder, Compression};
use rusqlite::Connection;

pub const PROGRESS_INTERVAL: usize = 1000;

pub struct Table {
    conn: Connection,
}

impl Table {
    /// Open a new database [Connection] and create a `molecules` table there.
    pub fn create(path: impl AsRef<Path>) -> rusqlite::Result<Self> {
        let conn = Connection::open(path)?;
        conn.execute(include_str!("create_table.sql"), ())?;
        Ok(Self { conn })
    }

    /// Open a [Connection] to an existing database at `path`.
    pub fn open(path: impl AsRef<Path>) -> rusqlite::Result<Self> {
        let conn = Connection::open(path)?;
        Ok(Self { conn })
    }

    /// Insert a single molecule entry into the database.
    pub fn insert_molecule(
        &self,
        smiles: String,
        moldata: String,
    ) -> rusqlite::Result<()> {
        let mut stmt =
            self.conn.prepare(include_str!("insert_molecule.sql"))?;
        stmt.execute((smiles, moldata))?;
        Ok(())
    }

    /// Insert a sequence of SMILES, moldata pairs into the database in a single
    /// transaction.
    pub fn insert_molecules(
        &mut self,
        mols: Vec<(String, String)>,
    ) -> rusqlite::Result<()> {
        let tx = self.conn.transaction()?;
        for (i, (smiles, moldata)) in mols.iter().enumerate() {
            let mut e = ZlibEncoder::new(Vec::new(), Compression::default());
            e.write_all(moldata.as_bytes()).unwrap();
            let moldata = e.finish().unwrap();
            tx.execute(include_str!("insert_molecule.sql"), (smiles, moldata))?;
            if i % PROGRESS_INTERVAL == 0 {
                eprint!("{i} complete\r");
            }
        }
        tx.commit()
    }

    /// Return a flattened vector of JSON molecule data strings from the
    /// database.
    pub fn get_moldata(&self) -> rusqlite::Result<Vec<String>> {
        let mut stmt = self.conn.prepare(include_str!("get_moldata.sql"))?;
        let ret = stmt
            .query_map([], |row| {
                let data: Vec<u8> = row.get(0).unwrap();
                let mut z = ZlibDecoder::new(&data[..]);
                let mut s = String::new();
                z.read_to_string(&mut s).unwrap();
                Ok(s)
            })
            .unwrap()
            .flatten()
            .collect();

        Ok(ret)
    }
}
