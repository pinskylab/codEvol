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
#SBATCH --mem=40G
#
# Number of cores:
#SBATCH --cpus-per-task=16


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load ngsLD/191108-GCC-8.2.0-2.31.1

# filter to gatk positions
cat data_31_01_20/GATK_filtered_SNP_set.tab | tail -n +2 | awk '{print $1"_"$2}' > tmp/gatk.beagle.pos # make the list of GATK positions in the right format
zcat data_31_01_20/Can_40_freq.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.beagle.pos | gzip > tmp/Can_40_freq.gatk.beagle.gz # filter to pos we want. Also removes header
zcat data_31_01_20/Can_14_freq.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.beagle.pos | gzip > tmp/Can_14_freq.gatk.beagle.gz 
zcat data_31_01_20/Lof_07_freq.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.beagle.pos | gzip > tmp/Lof_07_freq.gatk.beagle.gz 
zcat data_31_01_20/Lof_11_freq.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.beagle.pos | gzip > tmp/Lof_11_freq.gatk.beagle.gz 
zcat data_31_01_20/Lof_14_freq.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.beagle.pos | gzip > tmp/Lof_14_freq.gatk.beagle.gz 

# create position files w/out headers from the beagle files
# cuts off the first column, splits into two based on _, and separates with a tab
zcat tmp/Can_40_freq.gatk.beagle.gz | cut -f 1 | cut -d "_" -f 1,2 --output-delimiter $'\t' | gzip > tmp/can40.pos.gz
zcat tmp/Can_14_freq.gatk.beagle.gz | cut -f 1 | cut -d "_" -f 1,2 --output-delimiter $'\t' | gzip > tmp/can14.pos.gz
zcat tmp/Lof_07_freq.gatk.beagle.gz | cut -f 1 | cut -d "_" -f 1,2 --output-delimiter $'\t' | gzip > tmp/lof07.pos.gz
zcat tmp/Lof_11_freq.gatk.beagle.gz | cut -f 1 | cut -d "_" -f 1,2 --output-delimiter $'\t' | gzip > tmp/lof11.pos.gz
zcat tmp/Lof_14_freq.gatk.beagle.gz | cut -f 1 | cut -d "_" -f 1,2 --output-delimiter $'\t' | gzip > tmp/lof14.pos.gz

# find the number of sites for each population
sitescan40=$(zcat tmp/can40.pos.gz | wc -l)
sitescan14=$(zcat tmp/can14.pos.gz | wc -l)
siteslof07=$(zcat tmp/lof07.pos.gz | wc -l)
siteslof11=$(zcat tmp/lof11.pos.gz | wc -l)
siteslof14=$(zcat tmp/lof14.pos.gz | wc -l)

# calculations
# Can 40
ngsLD --geno tmp/Can_40_freq.gatk.beagle.gz \
	--probs \
	--n_ind 20 \
	--n_sites $sitescan40 \
	--pos tmp/can40.pos.gz \
	--max_kb_dist 10 \
	--n_threads 16 \
	--out analysis/ld.Can_40.gatk

# Can 14
ngsLD --geno tmp/Can_14_freq.gatk.beagle.gz \
	--probs \
	--n_ind 23 \
	--n_sites $sitescan14 \
	--pos tmp/can14.pos.gz \
	--max_kb_dist 10 \
	--n_threads 16 \
	--out analysis/ld.Can_14.gatk

# Lof 07
ngsLD --geno tmp/Lof_07_freq.gatk.beagle.gz \
	--probs \
	--n_ind 21 \
	--n_sites $siteslof07 \
	--pos tmp/lof07.pos.gz \
	--max_kb_dist 10 \
	--n_threads 16 \
	--out analysis/ld.Lof_07.gatk

# Lof 11
ngsLD --geno tmp/Lof_11_freq.gatk.beagle.gz \
	--probs \
	--n_ind 23 \
	--n_sites $siteslof11 \
	--pos tmp/lof11.pos.gz \
	--max_kb_dist 10 \
	--n_threads 16 \
	--out analysis/ld.Lof_11.gatk

# Lof 14
ngsLD --geno tmp/Lof_14_freq.gatk.beagle.gz \
	--probs \
	--n_ind 21 \
	--n_sites $siteslof14 \
	--pos tmp/lof14.pos.gz \
	--max_kb_dist 10 \
	--n_threads 16 \
	--out analysis/ld.Lof_14.gatk
	
	
# compress the output to save space
gzip analysis/ld.Can_40.gatk
gzip analysis/ld.Can_14.gatk
gzip analysis/ld.Lof_07.gatk
gzip analysis/ld.Lof_11.gatk
gzip analysis/ld.Lof_14.gatk

# clean up
rm tmp/*.pos.gz
rm tmp/*.gatk.beagle.gz
