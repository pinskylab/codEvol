#!/bin/bash

# Job name:
#SBATCH --job-name=HAPLO
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=50:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=30G
#
#SBATCH--ntasks-per-node=1

#if you need long time
# --partition=long

#If you need hugemem...

# --partition=hugemem

## Set up the job environment
source /cluster/bin/jobsetup
source /projects/cees/bin/paleomix/loadModules.sh

##HOME_DIR


##Copy data to $SCRATCH
cp $1* $SCRATCH
cd $SCRATCH


#HAPLOTYPE_CALLER
cp /projects/cees/databases/gadus_morhua/gadMor2/sep_mt_gen_nu_gen/gadMor2.nu.fasta $SCRATCH
cp /projects/cees/databases/gadus_morhua/gadMor2/sep_mt_gen_nu_gen/gadMor2.nu.dict $SCRATCH
cp /projects/cees/databases/gadus_morhua/gadMor2/sep_mt_gen_nu_gen/gadMor2.nu.fasta.fai $SCRATCH

ls -hal

java -version

which java

#echo "Initiating GATK on ${1}"
### $1 points to BAM file

/usr/bin/java -jar /projects/cees/bin/GATK/GenomeAnalysisTK.jar \
-T HaplotypeCaller -R gadMor2.nu.fasta \
-I ${1}.bam \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-o ${1}.g.vcf.gz 2> Haplotype_caller.${1}.out 



mkdir -p $SUBMITDIR/Report
mkdir -p $SUBMITDIR/GVCF


cp *.err $SUBMITDIR/Report/
cp *.txt $SUBMITDIR/Report/
cp *.out $SUBMITDIR/Report/
cp *.vcf.gz $SUBMITDIR/GVCF/
cp *.vcf.gz.tbi $SUBMITDIR/GVCF/





#DONE
echo "DONE!...."






