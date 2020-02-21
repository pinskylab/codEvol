#!/bin/bash

# Job name: Calculate pFst by population
#SBATCH --job-name=pFst
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-12:00:00
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

# get conda set up: https://documentation.sigma2.no/apps/userinstallsw.html
module load Anaconda3/2019.03 # for conda
export PS1=\$ # Set the ${PS1} (needed in the source of the Anaconda environment)
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh # Source the conda environment setup. The variable ${EBROOTANACONDA3} comes with the module load command

# Activate the environment by using the full path (not name)
# to the environment. The full path is listed if you do
# conda info --envs at the command prompt.
conda activate /cluster/home/mlpinsky/.conda/envs/testenv


# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# run pFst for Can
pFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 --background 24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44 --file data_31_01_20/Can.vcf.gz --deltaaf 0 --type GL >analysis/pfst_can.out

# pFst for Lof 07-11

# pFst for Lof 07-14

# pFst for Lof 11-14