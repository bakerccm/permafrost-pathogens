# run these commands from repo root to create the DAG, rulegraph and filegraph saved here

# the default rule 'all', plus the three rules 'descriptive', 'missingdata' and 'coinertia_etc'
# (which do not have output files), should be sufficient to re-run the code in its entirety:
snakemake -c 1 -n -F all

# generate DAG etc using:
snakemake -c 1 -n -F all --dag | dot -Tpdf >docs/snakemake_dag.pdf
snakemake -c 1 -n -F all --rulegraph | dot -Tpdf >docs/snakemake_rulegraph.pdf
snakemake -c 1 -n -F all --filegraph | dot -Tpdf >docs/snakemake_filegraph.pdf
