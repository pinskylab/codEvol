#!/bin/bash

# Job name: Theta by pop
#SBATCH --job-name=thetaPop
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-25:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=50G
#
# Number of cores:
#SBATCH --cpus-per-task=1


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load R/3.6.2-foss-2019b

## Run the calculations
Rscript --vanilla scripts/angsd_theta_bypop.R
