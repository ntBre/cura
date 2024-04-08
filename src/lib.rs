use std::{
    collections::{HashMap, HashSet},
    path::Path,
    process::exit,
    sync::atomic::{AtomicUsize, Ordering},
};

use log::{debug, info, trace};
use openff_toolkit::ForceField as OFF;
use rayon::iter::{
    IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use rdkit_rs::{
    bitvector::BitVector, find_smarts_matches_mol, fragment::recap_decompose,
    ROMol,
};

pub const PROGRESS_INTERVAL: usize = 1000;

pub mod parse;
pub mod query;
pub mod serve;
pub mod store;
pub mod table;

/// Parameter ID in a force field
type Pid = String;

/// SMIRKS string
type Smirks = String;

/// Database record
pub struct Molecule {
    /// None when inserting into the database but set when retrieving
    id: Option<usize>,
    smiles: String,
    inchikey: String,
    natoms: usize,
    elements: i128,
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Match {
    pub pid: Pid,
    pub smirks: Smirks,
    /// Database row IDs in the [Molecule]s table
    pub molecules: Vec<usize>,
}

/// Force field database record
#[derive(Clone, Debug)]
pub struct ForceField {
    pub id: Option<usize>,
    pub name: String,
    pub matches: Vec<Match>,
}

impl Molecule {
    pub fn new(
        smiles: String,
        inchikey: String,
        natoms: usize,
        elements: i128,
    ) -> Self {
        Self {
            id: None,
            smiles,
            inchikey,
            natoms,
            elements,
        }
    }
}

/// Load an OpenFF [ForceField] from `forcefield` and return a sequence of
/// parameter_id, SMIRKS pattern pairs corresponding to its `parameter_type`
/// [ParameterHandler].
fn load_forcefield(
    forcefield: &str,
    parameter_type: &str,
) -> Vec<(String, ROMol)> {
    OFF::load(forcefield)
        .unwrap()
        .get_parameter_handler(parameter_type)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), ROMol::from_smarts(&p.smirks())))
        .collect()
}

pub struct Report<'a> {
    pub args: Vec<String>,
    pub max: usize,
    pub nfps: usize,
    pub noise: usize,
    pub clusters: Vec<Vec<usize>>,
    pub mols: Vec<ROMol>,
    pub parameter: &'a str,
    pub map: HashMap<String, String>,
    pub mol_map: Vec<(String, ROMol)>,
}

impl Report<'_> {
    pub fn generate(&self, path: impl AsRef<Path>) -> std::io::Result<()> {
        use std::io::Write;
        let mut out = std::fs::File::create(path)?;
        let mut s = String::new();
        self.write(&mut s).unwrap();
        out.write_all(s.as_bytes())?;
        Ok(())
    }

    pub fn write(
        &self,
        mut out: impl std::fmt::Write,
    ) -> Result<(), std::fmt::Error> {
        writeln!(out, "<html>")?;
        writeln!(out, "<pre>args: {:?}</pre>", self.args)?;
        writeln!(
            out,
            "{nfps} molecules, {max} clusters, {noise} noise points, \
        pruned {} empty clusters",
            self.max + 1 - self.clusters.len(),
            nfps = self.nfps,
            max = self.max + 1,
            noise = self.noise
        )?;
        let pid = self.parameter;
        if let Some(smirks) = self.map.get(pid) {
            writeln!(out, "<p>PID: {pid}, SMIRKS: {smirks}</p>")?;
        }
        let mut clusters = self.clusters.clone();
        clusters.sort_by_key(|c| self.mols[c[0]].num_atoms());
        for (i, c) in clusters.iter().enumerate() {
            writeln!(out, "<h1>Cluster {}, {} molecules</h1>", i + 1, c.len())?;
            self.add_svg(&mut out, "Central Molecule", c[0])?;
        }
        writeln!(out, "</html>")?;
        Ok(())
    }

    pub fn add_svg(
        &self,
        out: &mut impl std::fmt::Write,
        msg: &str,
        idx: usize,
    ) -> Result<(), std::fmt::Error> {
        let mol = &self.mols[idx];
        let smile = mol.to_smiles();
        println!("{smile}");
        let svg = self.make_svg(mol);
        writeln!(out, "<p>{msg}</p>")?;
        writeln!(out, "<p>{} atoms</p>", mol.num_atoms())?;
        writeln!(out, "<p>SMILES: {smile}</p>")?;
        writeln!(out, "{svg}")?;
        Ok(())
    }

    pub fn make_svg(&self, mol: &ROMol) -> String {
        let mut hl_atoms = Vec::new();
        let pid = self.parameter;
        if self.map.contains_key(pid) {
            let tmp = find_matches_full(&self.mol_map, mol);
            let got = tmp.iter().find(|(_atoms, param_id)| param_id == &pid);
            if let Some((atoms, _pid)) = got {
                hl_atoms.clone_from(atoms);
            } else {
                panic!("smirks doesn't match any more");
            }
        }
        mol.draw_svg(400, 300, "", &hl_atoms)
    }
}

/// Compute Morgan fingerprints of size `radius` for each of the molecules in
/// `mols`.
pub fn make_fps(mols: &Vec<ROMol>, radius: u32) -> Vec<BitVector> {
    mols.par_iter()
        .map(|mol| mol.morgan_fingerprint_bit_vec::<1024>(radius))
        .collect()
}

pub fn fragment(mols: Vec<ROMol>) -> Vec<ROMol> {
    info!("starting fragment");

    // from Lily's example
    let dummy_replacements = [
        // Handle the special case of -S(=O)(=O)[*] -> -S(=O)(-[O-])
        (
            ROMol::from_smiles("S(=O)(=O)*"),
            ROMol::from_smiles("S(=O)([O-])"),
        ),
        // Handle the general case
        (ROMol::from_smiles("*"), ROMol::from_smiles("[H]")),
    ];

    /// the maximum molecule size to try fragmenting. 100 seems like a good
    /// limit, but maybe only because we're filtering a particular molecule with
    /// 103 atoms
    const MAX_FRAG_ATOMS: usize = 100;

    // apply the replacements above to a molecule and then round-trip through
    // SMILES to prevent radical formation issue
    let replace_fn = |mut m: ROMol| {
        for (inp, out) in &dummy_replacements {
            m = m.replace_substructs(inp, out, true).remove(0);
        }
        let smiles = m.to_smiles();
        ROMol::from_smiles(&smiles)
    };

    let ret = mols
        .into_par_iter()
        .flat_map(|mol| {
            let natoms = mol.num_atoms();
            if natoms > MAX_FRAG_ATOMS {
                debug!("filtered a molecule with {natoms} atoms");
                return vec![mol];
            }
            // we're looking for torsions, so the min_fragment_size is 4. this
            // probably isn't quite true because these might only be heavy atoms
            // at this point
            let leaves =
                recap_decompose(&mol, None, Some(4), None).get_leaves();
            if leaves.is_empty() {
                return vec![mol];
            }
            leaves.into_values().map(replace_fn).collect::<Vec<_>>()
        })
        .collect();

    info!("finished fragment");

    ret
}

/// Load a sequence of molecules from `smiles`, optionally fragmenting them with
/// the RECAP algorithm, and filtering out any larger than `max_atoms` or not
/// matching `mol_map[pid]`.
pub fn load_mols(
    smiles: Vec<&str>,
    max_atoms: usize,
    do_fragment: bool,
    pid: &str,
    mol_map: &[(String, ROMol)],
) -> Vec<ROMol> {
    // this gives each of the "fragments" from the original smiles
    let mut mols: Vec<_> = smiles.into_iter().map(ROMol::from_smiles).collect();

    info!("collected {} initial molecules", mols.len());

    if do_fragment {
        mols = fragment(mols);
    }

    let too_big = AtomicUsize::new(0);
    let no_match = AtomicUsize::new(0);
    let fragments = AtomicUsize::new(0);

    let mut ret: Vec<_> = mols
        .into_par_iter()
        .flat_map(|mut mol| {
            mol.openff_clean();
            fragments.fetch_add(1, Ordering::Relaxed);
            if mol.num_atoms() > max_atoms {
                too_big.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            let matches = find_matches_full(mol_map, &mol);
            if !matches.values().any(|p| p.as_str() == pid) {
                trace!(
                    "pid {pid} not found in {:?} for {}",
                    matches.into_values().collect::<HashSet<String>>(),
                    mol.to_smiles(),
                );
                no_match.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            Some((mol.to_smiles(), mol))
        })
        .collect();

    info!(
        "expanded initial smiles to {} fragments",
        fragments.into_inner()
    );

    let presort = ret.len();

    ret.sort_by_key(|(smiles, _mol)| smiles.clone());
    ret.dedup_by_key(|(smiles, _mol)| smiles.clone());

    info!(
        "filtered {} for size, {} for smirks, {} for duplicates",
        too_big.into_inner(),
        no_match.into_inner(),
        presort - ret.len(),
    );

    let (_, ret): (Vec<String>, _) = ret.into_iter().unzip();

    ret
}

/// returns a map of chemical environments in `mol` to their matching parameter
/// ids. matching starts with the first parameter and proceeds through the
/// whole sequence of parameters, so this should follow the SMIRNOFF typing
/// rules
pub fn find_matches_full(
    params: &[(String, ROMol)],
    mol: &ROMol,
) -> HashMap<Vec<usize>, String> {
    let mut matches = HashMap::new();
    for (id, smirks) in params {
        let env_matches = find_smarts_matches_mol(mol, smirks);
        for mut mat in env_matches {
            if mat.first().unwrap() > mat.last().unwrap() {
                mat.reverse();
            }
            matches.insert(mat, id.clone());
        }
    }
    matches
}

/// returns the set of parameter ids matching `mol`. matching starts with the
/// first parameter and proceeds through the whole sequence of parameters, so
/// this should follow the SMIRNOFF typing rules
pub fn find_matches(
    params: &[(String, ROMol)],
    mol: &ROMol,
) -> HashSet<String> {
    find_matches_full(params, mol).into_values().collect()
}

const PTABLE: [&str; 119] = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds ", "Rg ", "Cn ", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// Encode the elements in `mol` as a 128-bit set
pub fn get_elements(mol: &ROMol) -> i128 {
    let mut ret = 0;
    for i in mol.elements() {
        ret |= 1 << i;
    }
    ret
}

pub fn bits_to_elements(bits: i128) -> Vec<usize> {
    let mut ret = Vec::new();
    for i in 0..128 {
        if (bits & (1 << i)) != 0 {
            ret.push(i);
        }
    }
    ret
}

pub fn elements_to_bits(symbols: &[String]) -> i128 {
    let mut ret = 0;
    for s in symbols {
        let Some(pos) = PTABLE.iter().position(|sym| *sym == s.as_str()) else {
            eprintln!("unrecognized symbol {s}");
            exit(1);
        };
        ret |= 1 << pos;
    }
    ret
}

pub fn atomic_num_to_symbol(v: Vec<usize>) -> Vec<&'static str> {
    v.into_iter().map(|n| PTABLE[n]).collect()
}

#[cfg(test)]
mod tests {
    use crate::query::query;

    use self::{store::store, table::Table};

    use super::*;

    #[test]
    fn test_get_elements() {
        let mol = ROMol::from_smiles("CCO");
        let got = bits_to_elements(get_elements(&mol));
        let want = vec![6, 8];
        assert_eq!(got, want);
    }

    #[test]
    fn test_store() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut table = Table::create(tmp.path()).unwrap();
        // this file has multiple entries, but only one SMILES
        store(&mut table, "testfiles/small.sdf");
        let got = table.get_smiles().unwrap().len();
        let want = 1;
        assert_eq!(got, want);

        query(
            &mut table,
            "openff-2.1.0.offxml".to_string(),
            "ProperTorsions".to_string(),
            "testfiles/want.params".to_string(),
            &[],
        );

        let res = table.with_molecules(|mol| mol.natoms).unwrap();
        assert_eq!(res.len(), 1);
    }
}
