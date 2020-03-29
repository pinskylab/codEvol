#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/
# Takes one argument to designate the population comparison to make: CAN, LOF0711, LOF0714, or LOF1114
# expects a directory called tmp/

# Job name: Calculate bootstrapped Fst from ANGSD output
#SBATCH --job-name=bootFST
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=01-00:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=2G
#
# Number of cores:
#SBATCH --cpus-per-task=40

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load angsd/0.931-GCC-8.2.0-2.31.1


if [ $# -eq 0 ]
then
	echo No arguments supplied

elif [ $1 = "CAN" ]
then
	echo can
	# bootstrap and calculate the 2D SFS prior
	# -P: nthreads
	# -bootstrap: number of bootstrap replicates
#	realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -P 40 -bootstrap 100 >analysis/Can_40.Can_14.boot.ml

	#calculate fst
	# have to cut out each bootstrap separately
	for i in {1..100} 
	do
		cat analysis/Can_40.Can_14.boot.ml | head -n ${i} | tail -n 1 >tmp/thisbootcan.ml # cut out just one line
		realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -sfs tmp/thisbootcan.ml -P 40 -fstout tmp/thisbootcan # calculate FSTs
		realSFS fst stats tmp/thisbootcan.fst.idx -P 40 >>analysis/Can_40.Can_14.fsts.boot.txt # save the global estimate
	done
	
	rm tmp/thisbootcan.*

elif [ $1 == "LOF0711" ]
then
	echo lof0711

	realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 40 -bootstrap 100 >analysis/Lof_07.Lof_11.boot.ml

	for i in {1..100} 
	do
		cat analysis/Lof_07.Lof_11.boot.ml | head -n ${i} | tail -n 1 >tmp/thisbootlof0711.ml
		realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs tmp/thisbootlof0711.ml -P 40 -fstout tmp/thisbootlof0711
		realSFS fst stats tmp/thisbootlof0711.fst.idx -P 40 >>analysis/Lof_07.Lof_11.fsts.boot.txt
	done

	rm tmp/thisbootlof0711.*

elif [ $1 == "LOF0714" ]
then
	echo lof0714

	realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 40 -bootstrap 100 >analysis/Lof_07.Lof_14.boot.ml

	for i in {1..100} 
	do
		cat analysis/Lof_07.Lof_14.boot.ml | head -n ${i} | tail -n 1 >tmp/thisbootlof0714.ml
		realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs tmp/thisbootlof0714.ml -P 40 -fstout tmp/thisbootlof0714
		realSFS fst stats tmp/thisbootlof0714.fst.idx -P 40 >>analysis/Lof_07.Lof_14.fsts.boot.txt
	done

	rm tmp/thisbootlof0714.*

elif [ $1 == "LOF1114" ]
then
	echo lof1114

	realSFS data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 40 -bootstrap 100 >analysis/Lof_11.Lof_14.boot.ml

	for i in {1..100} 
	do
		cat analysis/Lof_11.Lof_14.boot.ml | head -n ${i} | tail -n 1 >tmp/thisbootlof1114.ml
		realSFS fst index data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs tmp/thisbootlof1114.ml -P 40 -fstout tmp/thisbootlof1114
		realSFS fst stats tmp/thisbootlof1114.fst.idx -P 40 >>analysis/Lof_11.Lof_14.fsts.boot.txt
	done

	rm tmp/thisbootlof1114.*

else
	echo Incorrect population specified. Options are CAN, LOF0711, LOF0714, or LOF1114
	
fi
