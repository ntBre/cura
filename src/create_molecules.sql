create table if not exists molecules (
id integer primary key,
smiles text unique,
inchikey text unique,
natoms integer,
elements blob
)
