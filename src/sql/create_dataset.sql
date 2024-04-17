create table if not exists dataset (
id integer primary key,
smiles text,
pid text,
CONSTRAINT unq UNIQUE (smiles, pid)
)
