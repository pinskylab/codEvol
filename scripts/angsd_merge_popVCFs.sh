#!/bin/bash

# Job name: Merge angsd vcf files into populations
#SBATCH --job-name=mergeVCF
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-08:00:00
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
module load BCFtools/1.9-intel-2018b

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# first have to index the vcf files if not already done
# tabix -p vcf data_31_01_20/Can_14_freq.vcf.gz
# tabix -p vcf data_31_01_20/Can_40_freq.vcf.gz
# tabix -p vcf data_31_01_20/Lof_07_freq.vcf.gz
# tabix -p vcf data_31_01_20/Lof_11_freq.vcf.gz
# tabix -p vcf data_31_01_20/Lof_14_freq.vcf.gz

# merge the population-level VCF files into one for Can and one for Lof
# also exclude sites with any missing data
# needed for running pFst later
bcftools merge --output-type v --force-samples data_31_01_20/Can_14_freq.vcf.gz data_31_01_20/Can_40_freq.vcf.gz | bcftools view --exclude 'GP="."' -o data_31_01_20/Can.vcf

bcftools merge --output-type v --force-samples data_31_01_20/Lof_07_freq.vcf.gz data_31_01_20/Lof_11_freq.vcf.gz data_31_01_20/Lof_14_freq.vcf.gz | bcftools view --exclude 'GP="."' -o data_31_01_20/Lof.vcf

#add GT columns (needed by pFst)
scripts/dummyGT.pl data_31_01_20/Can.vcf >data_31_01_20/Can2.vcf
scripts/dummyGT.pl data_31_01_20/Lof.vcf >data_31_01_20/Lof2.vcf

# compress
bgzip -c data_31_01_20/Can2.vcf >data_31_01_20/Can2.vcf.gz
bgzip -c data_31_01_20/Lof2.vcf >data_31_01_20/Lof2.vcf.gz

# index the merged files
tabix -p vcf data_31_01_20/Can2.vcf.gz
tabix -p vcf data_31_01_20/Lof2.vcf.gz
