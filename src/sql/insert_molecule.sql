-- insert a new molecule if the SMILES is not already present
insert or ignore into molecules (smiles, inchikey, natoms, elements, tag)
values (?1, ?2, ?3, ?4, ?5)
