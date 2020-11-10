#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/
# this assumes there are 20 replicates of the s0 simulations
# this assumes that analysis/slim_sim/slim_sim_n*.fst.csv.gz files exist from running slim_calcfst_reynolds.sh

# Job name: Calculate FST siteshuffle windowed outlier test using slim output
#SBATCH --job-name=fstslim
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=07-00:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=1G
#
# Number of cores:
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error



#merge each sim with 20 other chromosomes with s=0
for f in analysis/slim_sim/slim_sim_n*.fst.csv.gz #
do
	szero=$(echo "$f" | sed -e 's/_s0.[1-9]_\|_s1_\|_s1.[1-9]_\+/_s0_/g') # get name of the corresponding s0 files
	szero=${szero%_i*.fst.csv.gz} # remove iteration number and to get base name for all 20 iterations
	
	out=${f:18} # remove analysis/slim_sim/
	out=tmp/${out/.fst.csv.gz/_comb.fst.csv.gz} # create output filename

	outshuff=${out:4} # output file name from site shuffling (to test if already exists). cut off tmp/
	outshuff=analysis/slim_sim/${outshuff/.fst.csv.gz/.fst.siteshuffle.csv.gz} # set up output file name

	if [ ! -f "$outshuff" ]; then # only run if final shuffled file does not exist
		#echo $szero
		echo $out # a simple progress indicator
	
		# cat together the s0 files after removing headers
		# also replace the chromosome number with a made-up new one
		rm -f tmp/szero.fst.csv.gz # force remove if it exists
	
		chrom=2
		for s in $szero*[!b].fst.csv.gz # [!b] so that _comb files aren't also included
		do
			zcat $s | tail -n+2 | sed -e "s/^1/$chrom/g" | gzip >> tmp/szero.fst.csv.gz # replace CHROM number and append
			chrom=$((++chrom)) # move on to the next chromosome
		done
	
		zcat $f tmp/szero.fst.csv.gz | gzip > $out # cat file and the s0 files together together
	fi
done


# run fst site shuffle for single chromosomes
module --quiet purge
module load R/3.6.2-foss-2019b 

for f in analysis/slim_sim/slim_sim_n*.fst.csv.gz
do
	out=${f/.fst.csv.gz/.fst.siteshuffle.csv.gz} # set up output file name
	
	if [ ! -f "$out" ]; then # only run if final shuffled file does not exist
		echo $out
		
		Rscript --vanilla scripts/slim_fst_siteshuffle_null.r $f $out
	fi
done

# run fst site shuffle for full genomes
module --quiet purge
module load R/3.6.2-foss-2019b 

for f in tmp/slim_sim_n*_comb.fst.csv.gz
do
	out=${f:4} # cut off tmp/
	out=analysis/slim_sim/${out/.fst.csv.gz/.fst.siteshuffle.csv.gz} # set up output file name
	
	if [ ! -f "$out" ]; then # only run if final shuffled file does not exist
		echo $out
		
		Rscript --vanilla scripts/slim_fst_siteshuffle_null.r $f $out
	fi
done


# clean up the intermediate files
rm tmp/slim_sim_n*.fst.csv.gz
rm tmp/szero.fst.csv.gz
