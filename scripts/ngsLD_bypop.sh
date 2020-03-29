#!/bin/bash

# Job name: ngsLD calcs
#SBATCH --job-name=ngsLD
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=01-00:00:00
#
# Max memory usage:
#SBATCH --mem=10G
#
# Number of cores:
#SBATCH --cpus-per-task=40


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load ngsLD/191108-GCC-8.2.0-2.31.1

# create position files w/out headers
zcat data_31_01_20/Can_40_freq.pos.gz | tail -n +2 | gzip > tmp/can40.pos.gz
zcat data_31_01_20/Can_14_freq.pos.gz | tail -n +2 | gzip > tmp/can14.pos.gz
zcat data_31_01_20/Lof_07_freq.pos.gz | tail -n +2 | gzip > tmp/lof07.pos.gz
zcat data_31_01_20/Lof_11_freq.pos.gz | tail -n +2 | gzip > tmp/lof11.pos.gz
zcat data_31_01_20/Lof_14_freq.pos.gz | tail -n +2 | gzip > tmp/lof14.pos.gz

# calculations
# Can 40
ngsLD --geno data_31_01_20/Can_40_freq.beagle.gz \
	--probs \ # Beagle GP format for input
	--n_ind 20 \
	--n_sites 18535474 \
	--pos tmp/can40.pos.gz \
	--max_kb_dist 5 \ # out to 5kb
	--n_threads 40 \
	--out analysis/ld.Can_40

# Can 14
ngsLD --geno data_31_01_20/Can_14_freq.beagle.gz \
	--probs \ 
	--n_ind 23 \
	--n_sites 18932496 \
	--pos tmp/can14.pos.gz \
	--max_kb_dist 5 \ 
	--n_threads 40 \
	--out analysis/ld.Can_14

# Lof 07
ngsLD --geno data_31_01_20/Lof_07_freq.beagle.gz \
	--probs \ 
	--n_ind 21 \
	--n_sites 18543439 \
	--pos tmp/lof07.pos.gz \
	--max_kb_dist 5 \ 
	--n_threads 40 \
	--out analysis/ld.Lof_07

# Lof 11
ngsLD --geno data_31_01_20/Lof_11_freq.beagle.gz \
	--probs \ 
	--n_ind 23 \
	--n_sites 18954865 \
	--pos tmp/lof11.pos.gz \
	--max_kb_dist 5 \ 
	--n_threads 40 \
	--out analysis/ld.Lof_11

# Lof 14
ngsLD --geno data_31_01_20/Lof_14_freq.beagle.gz \
	--probs \ 
	--n_ind 21 \
	--n_sites 18957426 \
	--pos tmp/lof14.pos.gz \
	--max_kb_dist 5 \ 
	--n_threads 40 \
	--out analysis/ld.Lof_14
	
	
# clean up
rm tmp/*.pos.gz
