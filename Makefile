CHEMBL_BASE := /home/brent/omsf/chembl
CHEMBL := $(CHEMBL_BASE)/data/chembl_33.sdf

clippy:
	cargo clippy --workspace --tests

test:
	cargo test

# Usage:
# $(call run, SUBCOMMAND, ARGS...)
run = cargo run --release -- $1 $2

store:
	$(call run, store, $(CHEMBL))

query:
	$(call run, query, -s $(CHEMBL_BASE)/input/want.params -o output)
