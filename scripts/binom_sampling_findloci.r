# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
	require(plyr)
	require(parallel)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
}

# read in data
dat1 <- fread('analysis/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat1, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF07_11')) # for 1907 vs. 2011. sample sizes
dat2 <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat2, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF07_14')) # for 1907 vs. 2014. sample sizes
dat3 <- fread('analysis/Frequency_table_Lof11_Lof14.txt', header=TRUE); setnames(dat3, 3:7, c('N_CHR_11', 'Freq_11', 'N_CHR_14', 'Freq_14', 'ABS_DIFF11_14')) # for 2011 vs. 2014

# make sure rows line up
all(dat1$CHROM == dat2$CHROM) 
all(dat1$POS == dat2$POS) 
all(dat1$CHROM == dat3$CHROM) 
all(dat1$POS == dat3$POS) 

# calculate MAF
dat1[,maf07:=Freq_07]
dat1[maf07>0.5, maf07:=1-maf07]
dat1[,maf11:=Freq_11]
dat1[maf11>0.5, maf11:=1-maf11]

dat2[,maf07:=Freq_07]
dat2[maf07>0.5, maf07:=1-maf07]
dat2[,maf14:=Freq_14]
dat2[maf14>0.5, maf14:=1-maf14]

dat3[,maf11:=Freq_11]
dat3[maf11>0.5, maf11:=1-maf11]
dat3[,maf14:=Freq_14]
dat3[maf14>0.5, maf14:=1-maf14]

# find loci
inds07to11or14 <- dat1[,maf07 <= 0.05 & ABS_DIFF07_11>=0.25] & dat2[, maf07 <= 0.05 & ABS_DIFF07_14>=0.25] # fall in over-represented categories for 1907-2011 and 1907-2014
	sum(inds07to11or14) # 296
	
inds11to14 <- dat3[,(maf11 <= 0.05 & ABS_DIFF11_14>=0.05) | (maf11 <= 0.1 & maf11 > 0.05 & ABS_DIFF11_14>=0.15) | (maf11 <= 0.15 & maf11 > 0.1 & ABS_DIFF11_14 >= 0.3)] # fall in over-represented categories for 2011-2014 (doesn't get all of them, but gets the major ones)
	sum(inds11to14) # 19259

	sum(inds07to11or14 & !inds11to14) # still 296: no overlap

dat3[inds07to11or14, summary(maf1)] # allele freq in 2011 for candidate outliers
dat3[inds07to11or14, summary(ABS_DIFF)] # change from 2011 to 2014

	# trim to candidate loci
setkey(dat1, CHROM, POS, N_CHR_07, Freq_07, maf07)
setkey(dat2, CHROM, POS, N_CHR_07, Freq_07, maf07)
cands <- dat1[dat2,.(CHROM, POS, N_CHR_07, maf07, N_CHR_11, maf11, N_CHR_14, maf14, ABS_DIFF07_11, ABS_DIFF07_14)]
cands <- cands[inds07to11or14,]

t(t(cands[, table(CHROM)])) # distribution of candidates across linkage groups


# write out candidate loci
write.table(cands, file='analysis/Frequency_table_Lof07_Lof11_Lof14_candidates.tsv', quote=FALSE, row.names=FALSE, sep='\t')


# find distance among loci
setkey(cands, CHROM, POS) # sort by position
cands[,ndist := NA] # nearest neighbor at an earlier position (measure only to left)
for(i in 1:nrow(cands)){
	j <- which(cands$CHROM == cands$CHROM[i]) # other loci on same chromosome
	j <- j[j<i] # remove focal locus and any later loci
	if(length(j)>0) cands$ndist[i] <- min(abs(cands$POS[i] - cands$POS[j]))
} 