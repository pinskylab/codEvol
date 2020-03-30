#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Calculate sliding window FST for slim output with vcftools
#SBATCH --job-name=fstslim
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-3:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=1G
#
# Number of cores:
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-intel-2018b # for merging the vcf files
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0 # load vcftools


# first have to index the vcf files
for f in analysis/slim_sim/slim_sim_n*.vcf
do
	bgzip -f $f # have to bgzip to get tabix to work. -f to overwrite existing files
	tabix -p vcf $f.gz;
done

# merge the population-level VCF files into one for each sim (step through all gen 1 files)
for f in analysis/slim_sim/slim_sim_n*_1.vcf.gz # for all gen 1 files
do
	f2=${f/%_1.vcf.gz/_11.vcf.gz} # create filename for the later generation file
	out=${f/%_1.vcf.gz/_comb.bcf.gz} # create out filename
	#echo $f # a test
	#echo $f2 # a test
	#echo $out # a test
	bcftools merge --output-type b --force-samples -o $out $f $f2
done


# calc fst for output files
for f in analysis/slim_sim/slim_sim_n*_comb.bcf.gz
do
	# write file with list of individuals in early pop (first 21 for Lof 07)
	bcftools query -l $f | head -n 21 >tmp/early.txt

	# write file with list of individuals in late pop (second 23 for Lof 07)
	bcftools query -l $f | tail -n 23 >tmp/late.txt

	# sliding window fst. output file has suffix .weir.fst
	vcftools --bcf $f --weir-fst-pop tmp/early.txt --weir-fst-pop tmp/late.txt --fst-window-size 50000 --fst-window-step 10000 --out ${f%_comb.bcf.gz}
done

# clean up the intermediate files
rm analysis/slim_sim/slim_sim_n*_comb.bcf.gz
rm tmp/early.txt
rm tmp/late.txt