use std::path::Path;
use std::sync::mpsc::SyncSender;
use std::sync::Mutex;

use log::debug;
use rusqlite::params_from_iter;
use rusqlite::Connection;

use crate::atomic_num_to_symbol;
use crate::bits_to_elements;
use crate::ForceField;
use crate::Molecule;
use crate::Pid;
use crate::PROGRESS_INTERVAL;

pub(crate) use rusqlite::Result as RResult;

pub(crate) mod dataset;

pub struct Table {
    pub(crate) conn: Mutex<Connection>,
}

impl Table {
    /// Open a new database [Connection] and create a `molecules` table there.
    pub fn create(path: impl AsRef<Path>) -> RResult<Self> {
        let conn = Connection::open(path)?;
        conn.execute(include_str!("sql/create_molecules.sql"), ())?;
        conn.execute(include_str!("sql/create_forcefields.sql"), ())?;
        conn.execute(include_str!("sql/create_dataset.sql"), ())?;
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
        let mut stmt = conn.prepare(include_str!("sql/insert_molecule.sql"))?;
        stmt.execute((smiles, moldata))?;
        Ok(())
    }

    fn conn(&self) -> std::sync::MutexGuard<'_, Connection> {
        self.conn.lock().unwrap()
    }

    /// Insert a sequence of [Molecule]s into the database in a single
    /// transaction.
    pub fn insert_molecules(&mut self, mols: Vec<Molecule>) -> RResult<()> {
        let mut conn = self.conn.lock().unwrap();
        let tx = conn.transaction()?;
        for (
            i,
            Molecule {
                id,
                smiles,
                inchikey,
                natoms,
                elements,
                tag,
            },
        ) in mols.iter().enumerate()
        {
            assert!(
                id.is_none(),
                "attempted to insert existing molecule back into database"
            );
            tx.execute(
                include_str!("sql/insert_molecule.sql"),
                (smiles, inchikey, natoms, elements, tag),
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
        let mut stmt =
            conn.prepare(include_str!("sql/insert_forcefield.sql"))?;
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
        let mut stmt =
            conn.prepare(include_str!("sql/get_smiles_matching.sql"))?;
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

    pub fn get_forcefield(&self, ffname: &str) -> RResult<ForceField> {
        let conn = self.conn();
        let mut stmt =
            conn.prepare(include_str!("sql/get_smiles_matching.sql"))?;
        let res = stmt.query_row((ffname,), |row| {
            Ok(ForceField {
                id: row.get(0)?,
                name: row.get(1)?,
                matches: postcard::from_bytes(&row.get::<usize, Vec<u8>>(2)?)
                    .unwrap(),
            })
        })?;
        Ok(res)
    }

    pub fn get_forcefields(&self) -> RResult<Vec<ForceField>> {
        let conn = self.conn();
        let mut stmt = conn.prepare(include_str!("sql/get_forcefields.sql"))?;
        let res = stmt
            .query_map((), |row| {
                Ok(ForceField {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    matches: postcard::from_bytes(
                        &row.get::<usize, Vec<u8>>(2)?,
                    )
                    .unwrap(),
                })
            })?
            .flatten()
            .collect();
        Ok(res)
    }

    /// Return a flattened vector of SMILES from the database.
    pub fn get_smiles(&self) -> RResult<Vec<String>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("sql/get_smiles.sql"))?;
        let ret = stmt
            .query_map([], |row| Ok(row.get(0).unwrap()))
            .unwrap()
            .flatten()
            .collect();

        Ok(ret)
    }

    pub fn send_smiles(&self, sender: SyncSender<String>) -> RResult<()> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("sql/get_smiles.sql"))?;
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
        let mut stmt = conn.prepare(include_str!("sql/get_molecules.sql"))?;
        for s in stmt
            .query_map([], |row| {
                Ok(Molecule {
                    id: row.get(0)?,
                    smiles: row.get(1)?,
                    inchikey: row.get(2)?,
                    natoms: row.get(3)?,
                    elements: row.get(4)?,
                    tag: row.get(5)?,
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
        let mut stmt = conn.prepare(include_str!("sql/get_molecules.sql"))?;
        let ret = stmt
            .query_map([], |row| {
                Ok(f(Molecule {
                    id: row.get(0)?,
                    smiles: row.get(1)?,
                    inchikey: row.get(2)?,
                    natoms: row.get(3)?,
                    elements: row.get(4)?,
                    tag: row.get(5)?,
                }))
            })?
            .flatten()
            .collect();

        Ok(ret)
    }

    /// Return the number of molecules in the `molecules` table.
    pub fn count_molecules(&self) -> RResult<usize> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(include_str!("sql/get_molecules.sql"))?;
        let res = stmt.query_map([], |_row| Ok(()))?;
        Ok(res.count())
    }

    /// Delete the entry in the forcefields table named `forcefield`.
    pub fn reset_forcefield(&self, forcefield: &str) -> RResult<()> {
        let conn = self.conn();
        let mut stmt =
            conn.prepare(include_str!("sql/delete_forcefield.sql"))?;
        stmt.execute((forcefield,))?;
        Ok(())
    }

    /// Print the status of the database tables to stdout.
    pub fn print_status(&self) {
        let (elements, total_atoms, nmols) = self
            .with_molecules(|mol| (mol.elements, mol.natoms))
            .unwrap()
            .iter()
            .fold((0, 0, 0), |(e, a, n), x| (e | x.0, a + x.1, n + 1));
        println!("Molecules: {nmols}");
        println!(
            "Average size: {:.0} atoms",
            total_atoms as f64 / nmols as f64
        );
        println!(
            "Elements: {:?}",
            atomic_num_to_symbol(bits_to_elements(elements))
        );

        println!();

        let forcefields = self.get_forcefields().unwrap();
        println!("Force fields: {}", forcefields.len());
        for ff in forcefields {
            println!("{}: {} parameters", ff.name, ff.matches.len());
        }
    }
}
