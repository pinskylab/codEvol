#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

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
# These are standard, but they cause conda activate to fail
#set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

# get conda set up: https://documentation.sigma2.no/apps/userinstallsw.html
module load Anaconda3/2019.03 # for conda
export PS1=\$ # Set the ${PS1} (needed in the source of the Anaconda environment)
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh # Source the conda environment setup. The variable ${EBROOTANACONDA3} comes with the module load command

# Activate the environment by using the full path (not name)
# to the environment. The full path is listed if you do
# conda info --envs at the command prompt.
conda activate /cluster/home/mlpinsky/.conda/envs/testenv


# run pFst for Can
pFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 --background 24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44 --file data_31_01_20/Can2.vcf.gz --deltaaf 0 --type GP >analysis/pfst_can.out

# pFst for Lof 07-11
pFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21 --background 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45 --file data_31_01_20/Lof2.vcf.gz --deltaaf 0 --type GP >analysis/pfst_lof0711.out

# pFst for Lof 07-14
pFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21 --background 46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67 --file data_31_01_20/Lof2.vcf.gz --deltaaf 0 --type GP >analysis/pfst_lof0714.out

# pFst for Lof 11-14
pFst --target 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45 --background 46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67 --file data_31_01_20/Lof2.vcf.gz --deltaaf 0 --type GP >analysis/pfst_lof1114.out
