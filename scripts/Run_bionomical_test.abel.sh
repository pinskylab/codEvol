#!/bin/bash

# Job name:
#SBATCH --job-name=Binomial_test
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G
#
#SBATCH--ntasks-per-node=1

#if you need long time
# --partition=long

#If you need hugemem...

# --partition=hugemem

# Bastiaan's instructions:
# The script of running the heterozygosity on the before genotype depth filtered files is here
# Its an adaptation to run the script per chromosome and on Abel (to be faster, it takes some time). The script I call with a bash loop e.g.:  (for f in $(ls All_*HF);do sbatch Run_bionomical_test.abel.sh $f;done)
# This would submit all chromosome files simultaneously and works great when the script is in your path. I think its a fantastic ID to place all scripts in the Github area, I’ll go look in how committing works. The scripts may need some additional “readme” padding as most scripts I run are likely to be environment specific (with hard coded referrals to certain scripts or files) for some parts.

## Set up the job environment
source /cluster/bin/jobsetup


module purge
module load vcftools/0.1.14.zlib.1.2.8
module load R

#####Script to run a bionomial test of allele contribution on heterozygote genotype calls. A true mendelian SNP will fall within the 0.5 frequency
#### R script is made by Malin Pinsky, with input from Bastiaan Star
#### Use at own risk




echo "Extracting allele data..."

##Select output from VCF (allele data)
vcftools --gzvcf $1 \
--extract-FORMAT-info AD \
--out $1

echo "Extracting genotype data..."

##Select output from VCF (genotype data)
vcftools --gzvcf $1 \
--extract-FORMAT-info GT \
--out $1


echo "Running bionomical test on heterozygote genotype calls..."


##run bionomical test
Rscript /projects/cees/bin/filter_allele_balance.r $1.GT.FORMAT $1.AD.FORMAT $1

## Select the postions to EXCLUDE
awk '$5 < 0.05 {print $0}' $1.tsv > $1.exclude_binom_05.tsv


echo "Done!"
