#!/bin/bash#vcftools --gzvcf data_2018.09.13/All_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND.vcf.gz --bed Gadmor_unique_25k-mer_150bp_from_edge.bed --recode --recode-INFO-all --stdout | bgzip > Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz



#vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz --exclude-positions Top5_perc_All_variance.tab --recode --recode-INFO-all --stdout \
#| bgzip > Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz

#Log:
#VCFtools - 0.1.16
#(C) Adam Auton and Anthony Marcketta 2009

#Parameters as interpreted:
#	--gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz
#	--exclude-positions Top5_perc_All_variance.tab
#	--recode-INFO-all
#	--recode
#	--stdout
#After filtering, kept 113 out of 113 Individuals
#Outputting VCF file...
#After filtering, kept 784237 out of a possible 825514 Sites
#Run Time = 194.00 seconds


#vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz \
#--exclude-positions Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz.exclude_binom_05.tsv \
#--recode --recode-INFO-all --stdout \
# | bgzip > Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom.vcf.gz

#Log:
#VCFtools - 0.1.16
#(C) Adam Auton and Anthony Marcketta 2009

#Parameters as interpreted:
#	--gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz
#	--exclude-positions Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz.exclude_binom_05.tsv
#	--recode-INFO-all
#	--recode
#	--stdout

#After filtering, kept 113 out of 113 Individuals
#Outputting VCF file...
#After filtering, kept 720987 out of a possible 784237 Sites
#Run Time = 211.00 seconds

#zcat Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom.vcf.gz |awk '($4!="C" && $5 != "T" ) {print $0}' | \
#awk '($4!="G" && $5 != "A" ) {print $0}' | bcftools convert -o Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam.vcf.gz -O z

zcat Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom.vcf.gz | awk '{if ($4== "C" && $5 == "T" || $4== "T" && $5 == "C"  \
|| $4== "G" && $5 == "A" || $4== "A" && $5 == "G" ) {} else print $0}' | \
bcftools convert -o Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz -O z

# Job name:
#SBATCH --job-name=2N_GNTP
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=5G
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


module load bcftools/1.6

bcftools filter -i 'FS<60.0 && SOR<4 && MQRankSum>-12.5 && MQRankSum<12.5 && ReadPosRankSum>-4.0 && QD>2.0 && MQ>30.0' --SnpGap 10 -O z -o $1_HF $1 &> $1_HF.out

module load vcftools

vcftools --gzvcf $1_HF --remove-indels --min-alleles 2 --max-alleles 2 \
--maf 0.05 --max-maf 0.95 --exclude-bed Atlantic_cod_repeats.tab --recode --recode-INFO-all --stdout | bgzip > $1_HF_GQ.vcf.gz

mkdir -p per_pop

for f in $(ls pop*)


do
#Calculate HWE
vcftools --gzvcf ${1}_HF_GQ.vcf.gz --keep $f --hardy --out per_pop/$1_$f 2> per_pop/$1_$f.log

done


#rm keep_hwe_max_missing_list.txt


#Collect all SNPS with HWE > 0.001, and only take those in HWE 

cat per_pop/$1_*.hwe |grep -v POS |awk '$6 > 0.001 {print $1"\t"$2}' |sort -n | uniq -c |awk '$1>1 {print $2"\t"$3}' > keep_$1_hwe_max_missing_list.txt



vcftools --gzvcf $1_HF_GQ.vcf.gz --positions keep_hwe_max_missing_list.txt --recode --recode-INFO-all --stdout | bgzip > $1_HF_GQ_HWE_MISS.vcf.gz


vcftools --gzvcf $1_HF_GQ_HWE_MISS.vcf.gz --remove exclude_due_pop_missinging_MT --recode --recode-INFO-all --stdout | bgzip > $1_HF_GQ_HWE_MISS_IND.vcf.gz




vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND.vcf.gz --bed Gadmor_unique_25k-mer_150bp_from_edge.bed --recode --recode-INFO-all --stdout \
| bgzip > Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz


### Here postions first need to be obtained by the variance calculation script
vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz --exclude-positions Top5_perc_All_variance.tab --recode --recode-INFO-all --stdout \
| bgzip > Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz

### Here positions first need to be obtained by the bionomial distribution script.

vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz \
--exclude-positions Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR.vcf.gz.exclude_binom_05.tsv \
--recode --recode-INFO-all --stdout \
 | bgzip > Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom.vcf.gz

### excludes size potentially affected by aDNA degradation. 
zcat Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom.vcf.gz | awk '{if ($4== "C" && $5 == "T" || $4== "T" && $5 == "C"  \
|| $4== "G" && $5 == "A" || $4== "A" && $5 == "G" ) {} else print $0}' | \
bcftools convert -o Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz -O z
