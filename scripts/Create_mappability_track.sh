#!/bin/bash

# Job name:
#SBATCH --job-name=GEM_MAP
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=144:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=20G
#
#SBATCH--ntasks-per-node=1

#if you need long time
# --partition=long

#If you need hugemem...

# --partition=hugemem

## Set up the job environment
source /cluster/bin/jobsetup


#### Argument $1 gives the window size in base-pair. 


echo "Get reference"

ln -s /projects/cees/databases/gadus_morhua/gadMor2/sep_mt_gen_nu_gen/gadMor2.nu.fasta .

echo "Done"

echo "module load GEM"

module load GEM

echo "Done"


echo "Index reference"
gem-indexer -i gadMor2.nu.fasta -o gadMor2.nu.fasta


echo "Done"

echo ${1}

echo "Calculate Mappability"

gem-mappability -T 1 -I gadMor2.nu.fasta.gem -l $1 -m 0 -o Map_$1


echo "Done"


echo "Create wig track"
gem-2-wig -I gadMor2.nu.fasta.gem -i Map_${1}.mappability -o Map_${1}.wig


module load BEDops

echo "Create BED track"core_makefile.yaml


wig2bed < Map_${1}.wig.wig > Map_${1}.bed

echo "Done"
