CHEMBL := /home/brent/omsf/chembl/data/chembl_33.sdf

clippy:
	cargo clippy

run:
	cargo run --release --bin cura -- store $(CHEMBL)

query:
	cargo run --release --bin query
