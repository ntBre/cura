-- insert a new molecule if the SMILES is not already present
insert or ignore into molecules (smiles, moldata) values (?1, ?2)
