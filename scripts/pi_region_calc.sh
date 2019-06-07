#!/bin/bash

# calculate pi for regions in 1e4 windows
# to be run on a cod node
# not Unplaced, using VCF filtered for kmer25, individual read depth variance, depth, and 80% of individuals

module load vcftools

# using 2014 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --window-pi 10000 --keep data_2019_03_18/pop_LOF_S_14.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map -c | grep -v "nan" > analysis/LOF_S_14_10kb.windowed.pi

# using 2011 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --window-pi 10000 --keep data_2019_03_18/pop_LOF_S_11.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map -c | grep -v "nan" > analysis/LOF_S_11_10kb.windowed.pi

# using 1907 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --window-pi 10000 --keep data_2019_03_18/pop_LOF_07.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map -c | grep -v "nan" > analysis/LOF_07_10kb.windowed.pi

# using CAN 1940 samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --window-pi 10000 --keep data_2019_03_18/pop_CAN40.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map -c | grep -v "nan" > analysis/CAN40_10kb.windowed.pi

# using CAN Modern samples
vcftools --gzvcf data_2019_06_06/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_Readdepthvar5perc_Hetbinom_minDP7_MM.vcf.gz --window-pi 10000 --keep data_2019_03_18/pop_CANMod.txt --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map -c | grep -v "nan" > analysis/CANMod_10kb.windowed.pi
