-- insert a new molecule if the SMILES is not already present
insert or ignore into molecules (smiles, natoms, elements)
values (?1, ?2, ?3)
