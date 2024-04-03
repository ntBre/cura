CHEMBL_BASE := /home/brent/omsf/chembl
CHEMBL := $(CHEMBL_BASE)/data/chembl_33.sdf

clippy:
	cargo clippy --workspace --tests

test:
	cargo test

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

query_args := -s $(CHEMBL_BASE)/input/want.params -o output

query:
	$(call run, query, $(query_args))

parse:
	$(call run, parse, output/t18b.smiles -t t18b)

board:
	$(call run, board)

memprof.query:
	$(call memprof, query, $(query_args))
