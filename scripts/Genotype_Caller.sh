#!/bin/bash

# Job name:
#SBATCH --job-name=2N_GNTP
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=144:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=30G
#
#if you need long time
#   --partition=long

#If you need hugemem...

#S --partition=hugemem

#SBATCH --ntasks-per-node=1

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi



###$1 points to chromosome (we run this per chromosome in parallel). 

#cd $SCRATCH

/usr/bin/java -Xmx30g -jar /projects/cees/bin/GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R gadMor2.nu.fasta \
-L ${1} \
-V all_historic.list \
--max_alternate_alleles 3 \
-o All_${1}_no_clip.vcf.gz 2> ALL_${1}_no_clip_rerun.out

echo "Finished Genotype_GVCF ${1}....."


