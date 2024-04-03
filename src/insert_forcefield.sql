-- insert a new forcefield if the name is not already present
insert or ignore into forcefields (name, matches)
values (?1, ?2)
