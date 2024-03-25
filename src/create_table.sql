create table if not exists molecules (
id integer primary key,
smiles text unique,
moldata text
)
