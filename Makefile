CHEMBL_BASE := /home/brent/omsf/chembl
CHEMBL := $(CHEMBL_BASE)/data/chembl_33.sdf
bench := /home/brent/omsf/projects/benchmarking

ff := $(bench)/forcefields/tm-2.2.offxml

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

query_args := -s work/uncovered.dat -x 'inchi:work/inchis.dat' -x 'natoms:100'	\
-x 'elements:Cl, P, Br, I, H, C, B, Si, O, N, F, S' --reset -f $(ff)

query:
	$(call run, query, $(query_args))

status:
	$(call run, status)

serve:
	$(call run, serve, -f $(ff))

memprof.query:
	$(call memprof, query, $(query_args))
