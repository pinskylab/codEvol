#!/bin/bash

# Job name: Merge angsd vcf files into populations
#SBATCH --job-name=mergeVCF
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
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load angsd/0.931-GCC-8.2.0-2.31.1

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

# merge the population-level VCF files into one for Can and one for Lof
# useful for running pFst later
vcf-merge data_31_01_20/Can_14_freq.vcf.gz data_31_01_20/Can_40_freq.vcf.gz | bgzip -c > data_31_01_20/Can.vcf.gz

vcf-merge data_31_01_20/Lof_07_freq.vcf.gz data_31_01_20/Lof_11_freq.vcf.gz data_31_01_20/Lof_14_freq.vcf.gz | bgzip -c > data_31_01_20/Lof.vcf.gz