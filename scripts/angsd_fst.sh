#!/bin/bash

# Job name: Calculate Fst from ANGSD output (global, by site, and windowed)
#SBATCH --job-name=calcFST
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
# needs output from ngsLD_find_unlinked.r to run GATK unlinked calcs

#############################################
#calculate fst for all sites
#############################################
# calculate the 2D SFS prior
# -P: nthreads
realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -P 16 >analysis/Can_40.Can_14.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 >analysis/Lof_07.Lof_11.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 >analysis/Lof_07.Lof_14.ml
realSFS data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 >analysis/Lof_11.Lof_14.ml

realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Lof_07_freq.saf.idx -P 16 >analysis/Can_40.Lof_07.ml
realSFS data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 >analysis/Can_14.Lof_11.ml
realSFS data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 >analysis/Can_14.Lof_14.ml


#calculate fst by site
realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -sfs analysis/Can_40.Can_14.ml -fstout analysis/Can_40.Can_14
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Lof_07.Lof_11.ml -fstout analysis/Lof_07.Lof_11
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_07.Lof_14.ml -fstout analysis/Lof_07.Lof_14
realSFS fst index data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_11.Lof_14.ml -fstout analysis/Lof_11.Lof_14

realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Lof_07_freq.saf.idx -sfs analysis/Can_40.Lof_07.ml -fstout analysis/Can_40.Lof_07
realSFS fst index data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Can_14.Lof_11.ml -fstout analysis/Can_14.Lof_11
realSFS fst index data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Can_14.Lof_14.ml -fstout analysis/Can_14.Lof_14


#get the global estimate from all sites
realSFS fst stats analysis/Can_40.Can_14.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_11.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_14.fst.idx 
realSFS fst stats analysis/Lof_11.Lof_14.fst.idx 

realSFS fst stats analysis/Can_40.Lof_07.fst.idx 
realSFS fst stats analysis/Can_14.Lof_11.fst.idx 
realSFS fst stats analysis/Can_14.Lof_14.fst.idx 

					
#to get windowed fst for all sites, see instructions: http://www.popgen.dk/angsd/index.php/Fst
# -type 2 to start at pos=1
realSFS fst stats2 analysis/Can_40.Can_14.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Can_14.slide
realSFS fst stats2 analysis/Lof_07.Lof_11.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_11.slide
realSFS fst stats2 analysis/Lof_07.Lof_14.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_14.slide
realSFS fst stats2 analysis/Lof_11.Lof_14.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_11.Lof_14.slide

realSFS fst stats2 analysis/Can_40.Lof_07.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Lof_07.slide
realSFS fst stats2 analysis/Can_14.Lof_11.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_14.Lof_11.slide
realSFS fst stats2 analysis/Can_14.Lof_14.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_14.Lof_14.slide


# output A and B (numerator and denominator) by site for designing a site-shuffle null model
realSFS fst print analysis/Can_40.Can_14.fst.idx | gzip > analysis/Can_40.Can_14.fst.AB.gz
realSFS fst print analysis/Lof_07.Lof_11.fst.idx | gzip > analysis/Lof_07.Lof_11.fst.AB.gz
realSFS fst print analysis/Lof_07.Lof_14.fst.idx | gzip > analysis/Lof_07.Lof_14.fst.AB.gz
realSFS fst print analysis/Lof_11.Lof_14.fst.idx | gzip > analysis/Lof_11.Lof_14.fst.AB.gz


###################################################
#calculate fst for GATK nodam2 sites (inc. linked)
###################################################
# prep list of nodam2 sites
# needs to be tab-separated, no header. use only chr and pos, since otherwise angsd interprets cols 3 and 4 as ref and alt
cat data_2020.05.07/GATK_filtered_SNP_no_dam2.tab > tmp/gatk.tab # trim off the header from the list of good sites so that ANGSD can index it

# index so that we can trim by site
angsd sites index tmp/gatk.tab

# calculate the 2D SFS prior
realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -P 16 -sites tmp/gatk.tab  >analysis/Can_40.Can_14.gatk.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 -sites tmp/gatk.tab >analysis/Lof_07.Lof_11.gatk.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -sites tmp/gatk.tab >analysis/Lof_07.Lof_14.gatk.ml
realSFS data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -sites tmp/gatk.tab >analysis/Lof_11.Lof_14.gatk.ml

realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Lof_07_freq.saf.idx -P 16 -sites tmp/gatk.tab >analysis/Can_40.Lof_07.gatk.ml
realSFS data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 -sites tmp/gatk.tab >analysis/Can_14.Lof_11.gatk.ml
realSFS data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -sites tmp/gatk.tab >analysis/Can_14.Lof_14.gatk.ml


#calculate fst by site
realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -sfs analysis/Can_40.Can_14.gatk.ml -fstout analysis/Can_40.Can_14.gatk -sites tmp/gatk.tab
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Lof_07.Lof_11.gatk.ml -fstout analysis/Lof_07.Lof_11.gatk -sites tmp/gatk.tab
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_07.Lof_14.gatk.ml -fstout analysis/Lof_07.Lof_14.gatk -sites tmp/gatk.tab
realSFS fst index data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_11.Lof_14.gatk.ml -fstout analysis/Lof_11.Lof_14.gatk -sites tmp/gatk.tab

realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Lof_07_freq.saf.idx -sfs analysis/Can_40.Lof_07.gatk.ml -fstout analysis/Can_40.Lof_07.gatk -sites tmp/gatk.tab
realSFS fst index data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Can_14.Lof_11.gatk.ml -fstout analysis/Can_14.Lof_11.gatk -sites tmp/gatk.tab
realSFS fst index data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Can_14.Lof_14.gatk.ml -fstout analysis/Can_14.Lof_14.gatk -sites tmp/gatk.tab


#get the global estimate
realSFS fst stats analysis/Can_40.Can_14.gatk.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_11.gatk.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_14.gatk.fst.idx 
realSFS fst stats analysis/Lof_11.Lof_14.gatk.fst.idx 

realSFS fst stats analysis/Can_40.Lof_07.gatk.fst.idx 
realSFS fst stats analysis/Can_14.Lof_11.gatk.fst.idx 
realSFS fst stats analysis/Can_14.Lof_14.gatk.fst.idx 


#windowed fst
realSFS fst stats2 analysis/Can_40.Can_14.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Can_14.gatk.slide
realSFS fst stats2 analysis/Lof_07.Lof_11.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_11.gatk.slide
realSFS fst stats2 analysis/Lof_07.Lof_14.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_14.gatk.slide
realSFS fst stats2 analysis/Lof_11.Lof_14.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_11.Lof_14.gatk.slide

realSFS fst stats2 analysis/Can_40.Lof_07.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Lof_07.gatk.slide
realSFS fst stats2 analysis/Can_14.Lof_11.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_14.Lof_11.gatk.slide
realSFS fst stats2 analysis/Can_14.Lof_14.gatk.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_14.Lof_14.gatk.slide



#############################################
#calculate fst for GATK nodam2 unlinked sites
#############################################
# prep list of unlinked nodam2 sites
# needs to be tab-separated, no header. use only chr and pos, since otherwise angsd interprets cols 3 and 4 as ref and alt
zcat analysis/ld.unlinked.Can.gatk.nodam.csv.gz | tail -n +2 | sed 's/\,/\t/g' | cut -f 1-2 > tmp/ld.unlinked.Can.gatk.nodam.tab
zcat analysis/ld.unlinked.Lof.gatk.nodam.csv.gz | tail -n +2 | sed 's/\,/\t/g' | cut -f 1-2 > tmp/ld.unlinked.Lof.gatk.nodam.tab
zcat analysis/ld.unlinked.gatk.nodam.csv.gz | tail -n +2 | sed 's/\,/\t/g' | cut -f 1-2 > tmp/ld.unlinked.gatk.nodam.tab

# index so that we can trim by site
angsd sites index tmp/ld.unlinked.Can.gatk.nodam.tab
angsd sites index tmp/ld.unlinked.Lof.gatk.nodam.tab
angsd sites index tmp/ld.unlinked.gatk.nodam.tab

# calculate the 2D SFS prior
realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -P 16 -sites tmp/ld.unlinked.Can.gatk.nodam.tab  >analysis/Can_40.Can_14.gatk.unlinked.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 -sites tmp/ld.unlinked.Lof.gatk.nodam.tab >analysis/Lof_07.Lof_11.gatk.unlinked.ml
realSFS data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -sites tmp/ld.unlinked.Lof.gatk.nodam.tab >analysis/Lof_07.Lof_14.gatk.unlinked.ml
realSFS data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -sites tmp/ld.unlinked.Lof.gatk.nodam.tab >analysis/Lof_11.Lof_14.gatk.unlinked.ml

realSFS data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Lof_07_freq.saf.idx -P 16 -sites tmp/ld.unlinked.gatk.nodam.tab >analysis/Can_40.Lof_07.gatk.unlinked.ml
realSFS data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -P 16 -sites tmp/ld.unlinked.gatk.nodam.tab >analysis/Can_14.Lof_11.gatk.unlinked.ml
realSFS data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -P 16 -sites tmp/ld.unlinked.gatk.nodam.tab >analysis/Can_14.Lof_14.gatk.unlinked.ml


#calculate fst by site
realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Can_14_freq.saf.idx -sfs analysis/Can_40.Can_14.gatk.unlinked.ml -fstout analysis/Can_40.Can_14.gatk.unlinked -sites tmp/ld.unlinked.Can.gatk.nodam.tab
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Lof_07.Lof_11.gatk.unlinked.ml -fstout analysis/Lof_07.Lof_11.gatk.unlinked -sites tmp/ld.unlinked.Lof.gatk.nodam.tab
realSFS fst index data_31_01_20/Lof_07_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_07.Lof_14.gatk.unlinked.ml -fstout analysis/Lof_07.Lof_14.gatk.unlinked -sites tmp/ld.unlinked.Lof.gatk.nodam.tab
realSFS fst index data_31_01_20/Lof_11_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Lof_11.Lof_14.gatk.unlinked.ml -fstout analysis/Lof_11.Lof_14.gatk.unlinked -sites tmp/ld.unlinked.Lof.gatk.nodam.tab

realSFS fst index data_31_01_20/Can_40_freq.saf.idx data_31_01_20/Lof_07_freq.saf.idx -sfs analysis/Can_40.Lof_07.gatk.unlinked.ml -fstout analysis/Can_40.Lof_07.gatk.unlinked -sites tmp/ld.unlinked.gatk.unlinked.tab
realSFS fst index data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_11_freq.saf.idx -sfs analysis/Can_14.Lof_11.gatk.unlinked.ml -fstout analysis/Can_14.Lof_11.gatk.unlinked -sites tmp/ld.unlinked.gatk.unlinked.tab
realSFS fst index data_31_01_20/Can_14_freq.saf.idx data_31_01_20/Lof_14_freq.saf.idx -sfs analysis/Can_14.Lof_14.gatk.unlinked.ml -fstout analysis/Can_14.Lof_14.gatk.unlinked -sites tmp/ld.unlinked.gatk.unlinked.tab


#get the global estimate
realSFS fst stats analysis/Can_40.Can_14.gatk.unlinked.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_11.gatk.unlinked.fst.idx 
realSFS fst stats analysis/Lof_07.Lof_14.gatk.unlinked.fst.idx 
realSFS fst stats analysis/Lof_11.Lof_14.gatk.unlinked.fst.idx 

realSFS fst stats analysis/Can_40.Lof_07.gatk.unlinked.fst.idx 
realSFS fst stats analysis/Can_14.Lof_11.gatk.unlinked.fst.idx 
realSFS fst stats analysis/Can_14.Lof_14.gatk.unlinked.fst.idx 


#windowed fst
realSFS fst stats2 analysis/Can_40.Can_14.gatk.unlinked.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Can_14.gatk.unlinked.slide
realSFS fst stats2 analysis/Lof_07.Lof_11.gatk.unlinked.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_11.gatk.unlinked.slide
realSFS fst stats2 analysis/Lof_07.Lof_14.gatk.unlinked.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_07.Lof_14.gatk.unlinked.slide
realSFS fst stats2 analysis/Lof_11.Lof_14.gatk.unlinked.fst.idx -win 50000 -step 10000 -type 2 >analysis/Lof_11.Lof_14.gatk.unlinked.slide

realSFS fst stats2 analysis/Can_40.Lof_07.gatk.unlinked.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_40.Lof_07.gatk.unlinked.slide
realSFS fst stats2 analysis/Can_14.Lof_11.gatk.unlinked.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_14.Lof_11.gatk.unlinked.slide
realSFS fst stats2 analysis/Can_14.Lof_14.gatk.unlinked.fst.idx -win 50000 -step 10000 -type 2 >analysis/Can_14.Lof_14.gatk.unlinked.slide


###############
# clean up
###############
rm tmp/gatk.*
rm tmp/ld.unlinked.*