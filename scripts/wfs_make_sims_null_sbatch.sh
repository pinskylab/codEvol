#!/bin/bash

# Job name:
#SBATCH --job-name=MLPmakesims
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit: 10 hrs
#SBATCH --time=10:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=20G
#
# Number of cores:
#SBATCH --cpus-per-task=1


## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors
module load R/3.3.2

## Run the p-value calculations
# arguments: pop c1 c2 yr1 yr2 trimlowsampsize rerunlow repnum nsims
Rscript scripts/wfs_make_sims_null_function.r $1 $2 $3 $4 $5 $6 $7 $8 $9
