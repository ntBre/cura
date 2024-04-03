create table if not exists forcefields (
id integer primary key,
name text unique,
matches blob
)
