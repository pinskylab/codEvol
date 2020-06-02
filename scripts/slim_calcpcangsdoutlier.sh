#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Calculate pcangsd outlier test using slim output
#SBATCH --job-name=pcaslim
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
#SBATCH --cpus-per-task=2

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

# first combine the vcf files from the two timepoints and convert to beagle format
module --quiet purge  # Reset the modules to the system default
module load R/3.6.2-foss-2019b # to convert files to beagle format (using R)

for f in analysis/slim_sim/slim_sim_n*_1.vcf.gz # for all gen 1 files
do
	f2=${f/%_1.vcf.gz/_11.vcf.gz} # create filename for the later generation file
	out=${f:18} # strip off directory to create base of out filename
	out=tmp/${out/%_1.vcf.gz/.beagle.gz} # create out filename

	#echo $f2 # a test
	#echo $out # a test
	
	Rscript --vanilla scripts/vcftobeagle.R 0.005 $f $f2 $out
done


# run pcangsd
module --quiet purge
module load PCAngsd/200115-foss-2019a-Python-2.7.15 # 0.982

for f in tmp/slim_sim_n*.beagle.gz
do
	pref=${f:4} #set up output file prefix
	pref=analysis/slim_sim/${pref%.beagle.gz} # remove suffix
	
	# echo $pref # a test
	
	pcangsd.py -beagle $f -selection -minMaf 0.05 -threads 2 -o $pref -sites_save
done

# clean up the intermediate files
rm tmp/slim_sim_n*.beagle.gz
