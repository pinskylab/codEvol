#!/bin/bash

# Job name: Calculate Fst from ANGSD output
#SBATCH --job-name=calcFST
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-01:00:00
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

# calculate the 2D SFS prior
# -P: nthreads
realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -P 16 >analysis/Can_40.Can_14.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 >analysis/Lof_07.Lof_11.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 >analysis/Lof_07.Lof_14.ml
realSFS data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 >analysis/Lof_11.Lof_14.ml

#calculate fst by site for all sites
realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -sfs analysis/Can_40.Can_14.ml -fstout analysis/Can_40.Can_14
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Lof_07.Lof_11.ml -fstout analysis/Lof_07.Lof_11
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_07.Lof_14.ml -fstout analysis/Lof_07.Lof_14
realSFS fst index data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_11.Lof_14.ml -fstout analysis/Lof_11.Lof_14

#get the global estimate from all sites
realSFS fst stats analysis/Can_40.Can_14.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_11.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_14.fst.idx 
realSFS fst stats analysis/Lof_11.Lof_14.fst.idx 
					
#to get windowed fst for all sites, see instructions: http://www.popgen.dk/angsd/index.php/Fst
# -type 2 to start at pos=1
realSFS fst stats2 analysis/Can_40.Can_14.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Can_14.slide
realSFS fst stats2 analysis/Lof_07.Lof_11.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_11.slide
realSFS fst stats2 analysis/Lof_07.Lof_14.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_14.slide
realSFS fst stats2 analysis/Lof_11.Lof_14.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_11.Lof_14.slide

#calculate fst by site for GATK sites
cat data_31_01_20/GATK_filtered_SNP_set.tab | tail -n +2 > tmp/gatksites.tab # trim off the header from the list of good sites so that ANGSD can index it
angsd sites index tmp/gatksites.tab # index so that we can trim by site

realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -sfs analysis/Can_40.Can_14.ml -fstout analysis/Can_40.Can_14.gatk -sites tmp/gatksites.tab
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Lof_07.Lof_11.ml -fstout analysis/Lof_07.Lof_11.gatk -sites tmp/gatksites.tab
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_07.Lof_14.ml -fstout analysis/Lof_07.Lof_14.gatk -sites tmp/gatksites.tab
realSFS fst index data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_11.Lof_14.ml -fstout analysis/Lof_11.Lof_14.gatk -sites tmp/gatksites.tab

#get the global estimate from GATK sites
realSFS fst stats analysis/Can_40.Can_14.gatk.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_11.gatk.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_14.gatk.fst.idx 
realSFS fst stats analysis/Lof_11.Lof_14.gatk.fst.idx 

#to get windowed fst for GATK sites
realSFS fst stats2 analysis/Can_40.Can_14.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Can_14.gatk.slide
realSFS fst stats2 analysis/Lof_07.Lof_11.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_11.gatk.slide
realSFS fst stats2 analysis/Lof_07.Lof_14.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_14.gatk.slide
realSFS fst stats2 analysis/Lof_11.Lof_14.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_11.Lof_14.gatk.slide


# output A and B (numerator and denominator) by site for designing a null model
realSFS fst print analysis/Can_40.Can_14.fst.idx | gzip > analysis/Can_40.Can_14.fst.AB.gz
realSFS fst print analysis/Lof_07.Lof_11.fst.idx | gzip > analysis/Lof_07.Lof_11.fst.AB.gz
realSFS fst print analysis/Lof_07.Lof_14.fst.idx | gzip > analysis/Lof_07.Lof_14.fst.AB.gz
realSFS fst print analysis/Lof_11.Lof_14.fst.idx | gzip > analysis/Lof_11.Lof_14.fst.AB.gz


# clean up
rm tmp/gatksites.*