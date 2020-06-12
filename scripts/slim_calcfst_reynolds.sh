#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Calculate Reynolds FST components for slim output
#SBATCH --job-name=fstslim
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-12:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=1G
#
# Number of cores:
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error


# gzip any raw vcf files
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-intel-2018b # for merging the vcf files

shopt -s nullglob # so that a pattern that matches nothing "disappears", rather than treated as a literal string. good if there are no .vcf files left in the directory
for f in analysis/slim_sim/slim_sim_n*.vcf
do
	bgzip -f $f # have to bgzip. -f to overwrite existing files
done
shopt -u nullglob # unset this option


# calc for each simulation
module --quiet purge  # Reset the modules to the system default
module load R/3.6.2-foss-2019b # to convert files to beagle format (using R)

for f in analysis/slim_sim/slim_sim_n*_1.vcf.gz # for all gen 1 files
do
	f2=${f/%_1.vcf.gz/_11.vcf.gz} # create filename for the later generation file
	out=${f/%_1.vcf.gz/.fst.csv.gz} # set up output file prefix

	Rscript --vanilla scripts/fst_reynolds_fromvcf.R $f $f2 $out
done

