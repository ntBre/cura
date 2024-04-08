CHEMBL_BASE := /home/brent/omsf/chembl
CHEMBL := $(CHEMBL_BASE)/data/chembl_33.sdf

clippy:
	cargo clippy --workspace --tests

test:
	cargo test -- $(ARGS)

clean:
	rm try.sqlite

doc:
	cargo doc --open

# Usage:
# $(call run, SUBCOMMAND, ARGS...)
run = cargo run --release -- $1 $2

define memprof
	RUSTFLAGS='-g' cargo build --release
	valgrind -- target/release/cura $1 $2
endef

store:
	$(call run, store, $(CHEMBL))

query_args := -s $(CHEMBL_BASE)/input/want.params -x 'inchi:work/inchis.dat' -x	\
'natoms:100' -x 'elements:Cl, P, Br, I, H, C, B, Si, O, N, F, S' --reset

query:
	$(call run, query, $(query_args))

parse:
	$(call run, parse, output/t18b.smiles -t t18b)

status:
	$(call run, status)

serve:
	$(call run, serve)

memprof.query:
	$(call memprof, query, $(query_args))
