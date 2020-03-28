#!/bin/bash

####################################################################
## NOT FINISHED
## need *.thetas.idx files calculated by ANGSD before running this
## see steps 1 and 2 in http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests
####################################################################

# Job name: Calculate pi from ANGSD output
#SBATCH --job-name=calcPi
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-01:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=2G
#
# Number of cores:
#SBATCH --cpus-per-task=16

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load angsd/0.931-GCC-8.2.0-2.31.1

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/


#get the global estimate
thetaStat do_stat out.thetas.idx
					
#to get windowed
# -type 2 to start at pos=1
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -type 2 -outnames theta.thetasWindow.gz
