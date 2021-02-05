# to run on saga
# e.g.,  Rscript scripts/calc_allele_balance.r data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.GT.FORMAT data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.AD.FORMAT  data_2020.05.07/popCan13.txt analysis/allele_balance_popCan13
# run after making AD and GT files, with, e.g., 
# module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
# vcftools --gzvcf data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz --extract-FORMAT-info AD --out data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2
# vcftools --gzvcf data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz --extract-FORMAT-info GT --out data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2

require(data.table)


#accept arguments
args<-commandArgs(TRUE)

gt <- fread(args[1])
ad <- fread(args[2])
pop <- fread(args[3], header = FALSE) # the population file. a list of individual names.
outnm <- paste0(args[4],".csv.gz")


# simple check
dim(gt)
dim(ad)

# trim out individuals not in this pop
reminds <- setdiff(names(gt), c('CHROM', 'POS', pop$V1))
set(gt, j = reminds, value = NULL)  
set(ad, j = reminds, value = NULL)  

# merge GT and AD
setnames(gt, 3:ncol(gt), paste(names(gt)[3:ncol(gt)], '_gt', sep=''))
setkey(gt, CHROM, POS)
setkey(ad, CHROM, POS)
ad2 <- merge(ad, gt)

# mask uncalled or homozygote genotypes (set to NA)
indivs <- setdiff(names(ad), c('CHROM', 'POS'))
for(i in indivs){
	ad2[get(paste(i, '_gt', sep='')) %in% c('0/0', './.', '1/1'), eval(i):=NA]
}

# check that it worked
#	ad2[1:10, .(BM_026_gt, BM_026)]
#	ad2[BM_026_gt=='0/1', .(BM_026_gt, BM_026)]
#
#	ad2[1:10, .(BM_027_gt, BM_027)]
#	ad2[BM_027_gt=='0/1', .(BM_027_gt, BM_027)]

# remove gt columns (to reduce data in memory)
for(i in indivs){
	ad2[, eval(paste(i, '_gt', sep='')):=NULL]
}

# extract ref and alt alleles for het individuals (a couple seconds per individual)
for(i in indivs){
	cat(paste(i, ' ', sep=''))
	ad2[, eval(paste(i, '_ref', sep='')) := as.numeric(sapply(strsplit(get(i), split=','), '[', 1))]
	ad2[, eval(paste(i, '_alt', sep='')) := as.numeric(sapply(strsplit(get(i), split=','), '[', 2))]
}

# calc REF proportion
for(i in indivs){
  ad2[, eval(paste0(i, '_refprop')) := get(paste0(i, '_ref'))/(get(paste0(i, '_ref')) + get(paste0(i, '_alt')))]
}

# remove all non-refprop columns
for(i in indivs){
  ad2[, eval(i):=NULL]
  ad2[, eval(paste0(i, '_ref')):=NULL]
  ad2[, eval(paste0(i, '_alt')):=NULL]
}
set(ad2, j = c('CHROM', 'POS'), value = NULL)

# combine all non-NA values together
refprop <- numeric(0)
for(i in names(ad2)){
  theserefprops <- ad2[[i]]  
  refprop <- c(refprop, theserefprops[!is.na(theserefprops)])
}

# write out
write.table(refprop, file=gzfile(outnm), quote = FALSE, col.names = FALSE, row.names = FALSE)
