#!/bin/bash

# Job name:
#SBATCH --job-name=ANGSD_SNP
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=15:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=4G
#
#SBATCH--ntasks-per-node=14

#if you need long time
# --partition=long

#If you need hugemem...

# --partition=hugemem

module load angsd/0.931-GCC-8.2.0-2.31.1

##HOME_DIR


angsd -b bam.list -anc /cluster/projects/nn9244k/databases/gadus_morhua/gadMor2/sep_mt_gen_nu_gen/gadMor2.nu.fasta \
-out Freq/$1 \
-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 \
-doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 16 \
-minInd 3 -minQ 20 -sites SNP.list.txt\

