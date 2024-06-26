use crate::Pid;

use super::{RResult, Table};

impl Table {
    pub fn add_to_dataset(&self, smiles: String, pid: String) -> RResult<()> {
        let conn = self.conn();
        let mut stmt =
            conn.prepare(include_str!("../sql/add_to_dataset.sql"))?;
        stmt.execute((smiles, pid))?;
        Ok(())
    }

    /// Delete the entry in the forcefields table named `forcefield`.
    pub fn reset_dataset(&self) -> RResult<()> {
        let conn = self.conn();
        let mut stmt =
            conn.prepare(include_str!("../sql/delete_dataset.sql"))?;
        stmt.execute(())?;
        Ok(())
    }

    pub(crate) fn get_dataset_size(&self) -> RResult<usize> {
        let conn = self.conn();
        let mut stmt = conn.prepare("SELECT COUNT(*) from dataset")?;
        Ok(stmt.query_row((), |row| Ok(row.get(0).unwrap())).unwrap())
    }

    pub(crate) fn get_dataset_entries(&self) -> RResult<Vec<(String, Pid)>> {
        let conn = self.conn();
        let mut stmt = conn.prepare("SELECT smiles, pid from dataset")?;
        Ok(stmt
            .query_map((), |row| Ok((row.get(0).unwrap(), row.get(1).unwrap())))
            .unwrap()
            .flatten()
            .collect())
    }
}
