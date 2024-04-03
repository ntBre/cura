create table if not exists molecules (
id integer primary key,
smiles text unique,
natoms integer,
elements blob
)
