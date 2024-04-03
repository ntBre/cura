use std::path::Path;
use std::sync::mpsc::SyncSender;
use std::sync::Mutex;

use log::debug;
use rusqlite::params_from_iter;
use rusqlite::Connection;
use rusqlite::Result as RResult;

use crate::ForceField;
use crate::Molecule;
use crate::Pid;
use crate::PROGRESS_INTERVAL;

pub struct Table {
    pub(crate) conn: Mutex<Connection>,
}

impl Table {
    /// Open a new database [Connection] and create a `molecules` table there.
    pub fn create(path: impl AsRef<Path>) -> RResult<Self> {
        let conn = Connection::open(path)?;
        conn.execute(include_str!("create_molecules.sql"), ())?;
        conn.execute(include_str!("create_forcefields.sql"), ())?;
        Ok(Self {
            conn: Mutex::new(conn),
        })
    }

    /// Insert a single molecule entry into the database.
    pub fn insert_molecule(
        &self,
        smiles: String,
        moldata: String,
    ) -> RResult<()> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("insert_molecule.sql"))?;
        stmt.execute((smiles, moldata))?;
        Ok(())
    }

    fn conn(&self) -> std::sync::MutexGuard<'_, Connection> {
        self.conn.lock().unwrap()
    }

    /// Insert a sequence of SMILES, moldata pairs into the database in a single
    /// transaction.
    pub fn insert_molecules(&mut self, mols: Vec<Molecule>) -> RResult<()> {
        let mut conn = self.conn.lock().unwrap();
        let tx = conn.transaction()?;
        for (
            i,
            Molecule {
                id,
                smiles,
                natoms,
                elements,
            },
        ) in mols.iter().enumerate()
        {
            assert!(
                id.is_none(),
                "attempted to insert existing molecule back into database"
            );
            tx.execute(
                include_str!("insert_molecule.sql"),
                (smiles, natoms, elements),
            )?;
            if i % PROGRESS_INTERVAL == 0 {
                eprint!("{i} complete\r");
            }
        }
        tx.commit()
    }

    pub fn insert_forcefield(&mut self, forcefield: ForceField) -> RResult<()> {
        let conn = self.conn();
        let name = &forcefield.name;
        let blob = postcard::to_stdvec(&forcefield.matches).unwrap();
        let mut stmt = conn.prepare(include_str!("insert_forcefield.sql"))?;
        stmt.execute((name, blob))?;
        Ok(())
    }

    /// return the SMILES from the molecule table matching pid in the force
    /// field table
    pub fn get_smiles_matching(
        &self,
        ffname: &str,
        pid: &Pid,
    ) -> RResult<Vec<String>> {
        let conn = self.conn();
        let mut stmt = conn.prepare(include_str!("get_smiles_matching.sql"))?;
        let ids: Vec<usize> = stmt
            .query_map((ffname,), |row| {
                Ok(ForceField {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    matches: postcard::from_bytes(
                        &row.get::<usize, Vec<u8>>(2)?,
                    )
                    .unwrap(),
                })
            })?
            .flatten() // filter out any errors
            .flat_map(|ff| {
                ff.matches.into_iter().map(|m| {
                    if &m.pid == pid {
                        m.molecules
                    } else {
                        Vec::new()
                    }
                })
            })
            .flatten()
            .collect();

        if ids.is_empty() {
            debug!("no ids found for {pid} in {ffname}");
            return Ok(Vec::new());
        }

        // taken from rusqlite paramsfromiter example
        let mut vars = "?,".repeat(ids.len());
        vars.pop(); // remove trailing comma
        let sql =
            format!("select smiles from molecules where id in ({})", vars);
        let mut stmt = conn.prepare(&sql)?;
        let ret = stmt
            .query_map(params_from_iter(ids), |row| Ok(row.get(0).unwrap()))
            .unwrap()
            .flatten()
            .collect();
        Ok(ret)
    }

    /// Return a flattened vector of SMILES from the database.
    pub fn get_smiles(&self) -> RResult<Vec<String>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("get_smiles.sql"))?;
        let ret = stmt
            .query_map([], |row| Ok(row.get(0).unwrap()))
            .unwrap()
            .flatten()
            .collect();

        Ok(ret)
    }

    pub fn send_smiles(&self, sender: SyncSender<String>) -> RResult<()> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("get_smiles.sql"))?;
        for s in stmt
            .query_map([], |row| Ok(row.get(0).unwrap()))
            .unwrap()
            .flatten()
        {
            sender
                .send(s)
                .map_err(|_| rusqlite::Error::ExecuteReturnedResults)?;
        }
        Ok(())
    }

    pub fn send_molecules(&self, sender: SyncSender<Molecule>) -> RResult<()> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("get_molecules.sql"))?;
        for s in stmt
            .query_map([], |row| {
                Ok(Molecule {
                    id: row.get(0)?,
                    smiles: row.get(1)?,
                    natoms: row.get(2)?,
                    elements: row.get(3)?,
                })
            })?
            .flatten()
        {
            sender
                .send(s)
                .map_err(|_| rusqlite::Error::ExecuteReturnedResults)?;
        }
        Ok(())
    }

    pub fn with_molecules<T, F>(&self, mut f: F) -> RResult<Vec<T>>
    where
        F: FnMut(Molecule) -> T,
    {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("get_molecules.sql"))?;
        let ret = stmt
            .query_map([], |row| {
                Ok(f(Molecule {
                    id: row.get(0)?,
                    smiles: row.get(1)?,
                    natoms: row.get(2)?,
                    elements: row.get(3)?,
                }))
            })?
            .flatten()
            .collect();

        Ok(ret)
    }
}
