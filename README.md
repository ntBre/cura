# cura
Chemical dataset curation

## Usage

### Initialize the database

First, you must initialize the molecule database by running the following
command on one or more SDF files. For testing, I've been using the [ChEMBL
33][chembl] database.

``` shell
cura -d DATABASE store MOLECULE_FILE
```

### Query for OpenFF parameters

Next, you can search for molecules in the existing `DATABASE` with the `query`
subcommand:

``` shell
cura -d DATABASE query -f FORCEFIELD -s PARAM_FILE
```

`FORCEFIELD` should be the name of an OpenFF-compatible force field (either
built-in like `openff-2.1.0.offxml` or a local file), and `PARAM_FILE` should be
a file containing a sequence of parameter IDs, one per line. `cura` defaults to
searching for `ProperTorsions`, but this can be changed to another type of
`ParameterHandler`with the `--parameter-type/-p` flag. You can also apply
filters at this stage, using the `--filters/-x` flag. Currently supported
filters are:

| Flag             | Description                                               |
|------------------|-----------------------------------------------------------|
| `inchi:FILENAME` | Filter out molecules matching the InchiKeys in `FILENAME` |
| `natoms:N`       | Filter out molecules with more than `N` atoms             |


For example, to filter out InchiKeys in `inchi.dat` and molecules with more than
100 atoms:

``` shell
cura -d DATABASE query -f FORCEFIELD -s PARAM_FILE -x inchi:inchi.dat -x natoms:100
```

The query results are also stored back in the database, with separate results
tables for each unique force field name. If you want to re-run a query for a
given force field, overwriting the previous results, you can use the
`--reset/-r` flag.

[chembl]: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_33/
