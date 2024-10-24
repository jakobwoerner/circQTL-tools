#!/bin/bash

start_time=$(date +%s)

Rscript path/to/circQTL.R chr start_interval end_interval tissue_name covariates_file time_file ge_file snp_file disttogene out_dir

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))

eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
