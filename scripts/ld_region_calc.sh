#!/bin/bash

# calculate LD for regions in 1e4 windows
# to be run on a cod node
# not Unplaced, using VCF filtered for kmer25, individual read depth variance, depth, and 80% of individuals

module load vcftools

windsz=30000
windnm=30kb

# each of these vcftools runs takes ~5 min
# using 2014 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --geno-r2 --ld-window-bp $windsz --keep data_2019_03_18/pop_LOF_S_14.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map --out analysis/LOF_S_14_$windnm

# using 2011 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --geno-r2 --ld-window-bp $windsz --keep data_2019_03_18/pop_LOF_S_11.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map --out analysis/LOF_S_11_$windnm

# using 1907 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --geno-r2 --ld-window-bp $windsz --keep data_2019_03_18/pop_LOF_07.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map --out analysis/LOF_07_$windnm

# using CAN 1940 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --geno-r2 --ld-window-bp $windsz --keep data_2019_03_18/pop_CAN40.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map --out analysis/CAN40_$windnm

# using CAN Modern samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --geno-r2 --ld-window-bp $windsz --keep data_2019_03_18/pop_CANMod.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map --out analysis/CANMod_$windnm

# zip up the files to save space
gzip analysis/LOF_S_14_$windnm.geno.ld
gzip analysis/LOF_S_11_$windnm.geno.ld
gzip analysis/LOF_07_$windnm.geno.ld
gzip analysis/CAN40_$windnm.geno.ld
gzip analysis/CANMod_$windnm.geno.ld
