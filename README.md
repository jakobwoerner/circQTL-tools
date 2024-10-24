Analysis tools for circadian QTL discovery

To run the circQTL pipeline, use `circQTL.R` with the following arguments in order

* known_chr
* start_interval
* end_interval
* tissue_name
* covariates_file
* time_file
* ge_file
* snp_file (with ref/alt)
* disttogene
* out_dir
Note that without 10 arguments, the script will not run.

Use `run_circQTL.sh` as an example, and `run_adipose_chr22.sh` as a specific Voight lab example.
