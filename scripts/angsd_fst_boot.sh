#!/bin/bash

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
#SBATCH --cpus-per-task=16

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load angsd/0.931-GCC-8.2.0-2.31.1

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# bootstrap and calculate the 2D SFS prior
# -P: nthreads
# -bootstrap: number of bootstrap replicates
realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -P 16 -bootstrap 100 >analysis/Can_40.Can_14.boot.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 -bootstrap 100 >analysis/Lof_07.Lof_11.boot.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -bootstrap 100 >analysis/Lof_07.Lof_14.boot.ml
realSFS data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -bootstrap 100 >analysis/Lof_11.Lof_14.boot.ml

#calculate fst
# have to cut out each bootstrap separately
# CAN first
#for i in {1..100} 
#do
#	cat analysis/Can_40.Can_14.boot.ml | head -n ${i} | tail -n 1 >tmp/thisboot.ml # cut out just one line
#	realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -sfs tmp/thisboot.ml -P 16 -fstout tmp/thisboot # calculate FSTs
#	realSFS fst stats tmp/thisboot.fst.idx -P 16 >>analysis/Can_40.Can_14.fsts.boot.txt # save the global estimate
#done

# LOF 07-11
#for i in {1..100} 
#do
#	cat analysis/Lof_07.Lof_11.boot.ml | head -n ${i} | tail -n 1 >tmp/thisboot.ml # just one line
#	realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs tmp/thisboot.ml -P 16 -fstout tmp/thisboot
#	realSFS fst stats tmp/thisboot.fst.idx -P 16 >>analysis/Lof_07.Lof_11.fsts.boot.txt
#done

# LOF 07-14
#for i in {1..100} 
#do
#	cat analysis/Lof_07.Lof_14.boot.ml | head -n ${i} | tail -n 1 >tmp/thisboot.ml # just one line
#	realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs tmp/thisboot.ml -P 16 -fstout tmp/thisboot
#	realSFS fst stats tmp/thisboot.fst.idx -P 16 >>analysis/Lof_07.Lof_14.fsts.boot.txt
#done

# LOF 11-14
#for i in {1..100} 
#do
#	cat analysis/Lof_11.Lof_14.boot.ml | head -n ${i} | tail -n 1 >tmp/thisboot.ml # just one line
#	realSFS fst index data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs tmp/thisboot.ml -P 16 -fstout tmp/thisboot
#	realSFS fst stats tmp/thisboot.fst.idx -P 16 >>analysis/Lof_11.Lof_14.fsts.boot.txt
#done
