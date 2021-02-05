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



module purge
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load R/3.5.1-foss-2018b

#####Script to run a bionomial test of allele contribution (read depth of alleles) on heterozygote genotype calls. A true mendelian SNP will fall within the 0.5 frequency
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


echo "Running bionomial test on heterozygote genotype calls..."


##run bionomial test
Rscript /cluster/projects/nn9244k/databases/old_bin_CEES/filter_allele_balance.r $1.GT.FORMAT $1.AD.FORMAT $1

## Select the postions to EXCLUDE (significant test)
awk '$5 < 0.05 {print $0}' $1.tsv > $1.exclude_binom_05.tsv


echo "Done!"
