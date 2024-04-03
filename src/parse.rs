use clustrs::{dbscan, Label};
use openff_toolkit::ForceField;
use rayon::iter::{
    IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use rdkit_rs::{
    bitvector::BitVector, fingerprint::tanimoto, fragment::recap_decompose,
    ROMol,
};
use std::{
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};

use log::{debug, info, trace};

use crate::{find_matches_full, table::Table};

type Pid = String;
type Smirks = String;

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
    let overlap = AtomicUsize::new(0);
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
            if !matches.iter().any(|(_, p)| p.as_str() == pid) {
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
        "filtered {} for size, {} for smirks, {} for inchi, {} for duplicates",
        too_big.into_inner(),
        no_match.into_inner(),
        overlap.into_inner(),
        presort - ret.len(),
    );

    let (_, ret): (Vec<String>, _) = ret.into_iter().unzip();

    ret
}

pub fn parse(
    _table: &mut Table, // TODO use this
    input: PathBuf,
    forcefield: String,
    parameter_type: String,
    target: String,
) {
    let s = std::fs::read_to_string(&input)
        .unwrap_or_else(|e| panic!("failed to read {:?} for {e}", input));
    let smiles: Vec<_> = s.lines().collect();

    info!("{} initial smiles", smiles.len());

    let parameters = ForceField::load(&forcefield)
        .unwrap()
        .get_parameter_handler(&parameter_type)
        .unwrap()
        .parameters();

    let map: HashMap<Pid, Smirks> =
        parameters.iter().map(|p| (p.id(), p.smirks())).collect();

    let mol_map: Vec<(Pid, ROMol)> = parameters
        .into_iter()
        .map(|p| (p.id(), ROMol::from_smarts(&p.smirks())))
        .collect();

    const MAX_ATOMS: usize = 80;
    const DO_FRAGMENT: bool = false;
    const MORGAN_RADIUS: u32 = 4;

    let mols = load_mols(smiles, MAX_ATOMS, DO_FRAGMENT, &target, &mol_map);

    info!("making fingerprints");

    let fps: Vec<_> = make_fps(&mols, MORGAN_RADIUS);
    let nfps = fps.len();
    let distance_fn = |i, j| {
        if i == j {
            0.0
        } else {
            1.0 - tanimoto(&fps[i], &fps[j])
        }
    };

    info!("running DBSCAN");

    const DBSCAN_EPS: f64 = 0.5;
    const DBSCAN_MIN_PTS: usize = 1;

    let labels = dbscan(nfps, nfps, distance_fn, DBSCAN_EPS, DBSCAN_MIN_PTS);

    let max = match labels
        .iter()
        .filter_map(|l| match l {
            Label::Cluster(n) => Some(n),
            _ => None,
        })
        .max()
    {
        Some(n) => *n,
        None => {
            dbg!(labels);
            eprintln!("error: all noise points, exiting");
            std::process::exit(1);
        }
    };

    let mut clusters: Vec<Vec<usize>> = vec![vec![]; max + 1];
    let mut noise = 0;
    for (i, l) in labels.iter().enumerate() {
        match l {
            Label::Cluster(n) => clusters[*n].push(i),
            _ => noise += 1,
        }
    }

    info!("generating report");

    let output = input.with_extension("html");
    Report {
        args: std::env::args().collect::<Vec<_>>(),
        max,
        nfps,
        noise,
        clusters,
        mols,
        parameter: &target,
        map,
        mol_map,
    }
    .generate(output)
    .unwrap();
}
