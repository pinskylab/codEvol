#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/
# this assumes there are 20 replicates of the s0 simulations

# Job name: Calculate pcangsd outlier test using slim output
#SBATCH --job-name=pcaslim
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-21:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=1G
#
# Number of cores:
#SBATCH --cpus-per-task=2

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

# gzip any raw vcf files
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-intel-2018b # for merging the vcf files

shopt -s nullglob # so that a pattern that matches nothing "disappears", rather than treated as a literal string. good if there are no .vcf files left in the directory
for f in analysis/slim_sim/slim_sim_n*.vcf # for all vcf files
do
	bgzip -f $f # -f to overwrite existing files
done

shopt -u nullglob # unset this option

# first combine the vcf files from the two timepoints and convert to beagle format
# only run if final file does not exist
module --quiet purge  # Reset the modules to the system default
module load R/3.6.2-foss-2019b # to convert files to beagle format (using R)

for f in analysis/slim_sim/slim_sim_n*_1.vcf.gz # for all gen 1 files
do
	f2=${f/%_1.vcf.gz/_11.vcf.gz} # create filename for the later generation file
	outbase=${f:18} # strip off directory to create base of out filename
	out=tmp/${outbase/%_1.vcf.gz/.beagle.gz} # create out filename

	#echo $f2 # a test
	#echo $out # a test

	finalfile=analysis/slim_sim/${outbase/%_1.vcf.gz/.selection.npy} # add directory and suffix for the pcangsd file eventually to be created

	if [ ! -f "$finalfile" ]; then # only run if final file doesn't exist
		Rscript --vanilla scripts/vcftobeagle.R 0.005 $f $f2 $out # use genotype error of 0.005
	fi
done

#merge each sim with 20 other chromosomes with s=0
# only run if final file does not exist
for f in tmp/slim_sim_n*[!b].beagle.gz # [!b] so that _comb files aren't also included
do
	szero=$(echo "$f" | sed -e 's/_s0.[1-9]_\|_s1_\|_s1.[1-9]_\+/_s0_/g') # get name of the corresponding s0 files
	szero=${szero%_i*.beagle.gz} # remove iteration number and to get base name for all 20 iterations
	
	f2=${f/.beagle.gz/_comb.beagle.gz} # create output filename
	
	finalbase=${f2:4} # strip off directory to create base of out filename
	finalfile=analysis/slim_sim/${finalbase/%.beagle.gz/.selection.npy} # add directory and suffix for the pcangsd file eventually to be created

	if [ ! -f "$finalfile" ]; then # only run if final file doesn't exist
		#echo $szero
		echo $f2 # a simple progress indicator
	
		# cat together the s0 files after removing headers
		# also replace the chromosome number with a made-up new one
		# really only need to do this once for each combination of ne and f, so could streamline
		rm -f tmp/szero.beagle.gz # force remove if it exists
	
		chrom=2
		for s in $szero*[!b].beagle.gz # [!b] so that _comb files aren't also included
		do
			zcat $s | tail -n+2 | sed -e "s/^1/$chrom/g" | gzip >> tmp/szero.beagle.gz
			chrom=$((++chrom)) # move on to the next chromosome
		done
	
		zcat $f tmp/szero.beagle.gz | gzip > $f2 # cat file and the s0 files together together
	fi
done


# run pcangsd
# only run if final file does not exist
module --quiet purge
module load PCAngsd/200115-foss-2019a-Python-2.7.15 # 0.982

for f in tmp/slim_sim_n*.beagle.gz # will get both the single-chromosome sims and the _comb sims with s=0 chromosomes added
do
	pref=${f:4} # remove tmp/ to set up output file prefix
	pref=analysis/slim_sim/${pref%.beagle.gz} # remove suffix from output file prefix
	
	# echo $pref # a test

	finalfile=${pref}.selection.npy # add suffix for the pcangsd file eventually to be created

	if [ ! -f "$finalfile" ]; then # only run if final file doesn't exist
		pcangsd.py -beagle $f -selection -minMaf 0.05 -threads 2 -o $pref -sites_save
	fi
done

# clean up the intermediate files
rm tmp/slim_sim_n*.beagle.gz
rm tmp/szero.beagle.gz
