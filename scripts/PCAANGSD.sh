#!/bin/bash

# Job name:
#SBATCH --job-name=PCAANGSD
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=02:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=5G
#
#SBATCH --ntasks-per-node=2


module purge

module load PCAngsd/200115-foss-2019a-Python-2.7.15


pcangsd.py -beagle Ind_GATK_DP5_no_dam2.ex_inv.beagle.gz -sites_save \
-o Ind_PCA_DP5_no_dam2 \
-threads 2 > Ind_ANGSD_DP5_no_dam2.log


pcangsd.py -beagle Ind_GATK_DP9_no_dam2.ex_inv.beagle.gz -sites_save \
-o Ind_PCA_DP9_no_dam2 \
-threads 2 > Ind_ANGSD_DP9_no_dam2.log
