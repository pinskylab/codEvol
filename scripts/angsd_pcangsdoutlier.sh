#!/bin/bash
# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Run pcangsd selection scan by population
#SBATCH --job-name=pca_sel
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-06:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=16G
#
# Number of cores:
#SBATCH --cpus-per-task=2

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load PCAngsd/200115-foss-2019a-Python-2.7.15 # 0.982

# Canada
pcangsd.py -beagle data_31_01_20/Can_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can -sites_save

# Lof 07-11
pcangsd.py -beagle data_31_01_20/Lof0711_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711 -sites_save

# Lof 07-14
pcangsd.py -beagle data_31_01_20/Lof0714_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714 -sites_save

# Lof 11-14
pcangsd.py -beagle data_31_01_20/Lof1114_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114 -sites_save