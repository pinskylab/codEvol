#!/bin/bash

# specify two arguments: 
# First: 0 for all loci, 1 for gatk loci only
# 	only 0 works right now
# Second: 1 for Can, 2 for Lof0711, 3 for Lof0714, 4 for Lof1114

# Job name: Shuffle Theta windows
#SBATCH --job-name=thetaShuffle
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-20:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=40G
#
# Number of cores:
#SBATCH --cpus-per-task=1


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load R/3.6.2-foss-2019b

## Run the calculations
Rscript --vanilla scripts/angsd_theta_siteshuffle_null.r $1 $2
