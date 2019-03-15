#!/bin/bash

# Job name:
#SBATCH --job-name=Malindefault
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=120:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=2G
#
# Number of cores:
#SBATCH --cpus-per-task=16


## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors
module load R/3.3.2

## Run the p-value calculations
# Arguments: myalcnt1 myalcnt2 pop myyr1 myyr2 maxcores
Rscript scripts/wfs_nullmodel_function.r $1 $2 $3 $4 $5 $6
