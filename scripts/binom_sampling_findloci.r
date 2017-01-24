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
dat1 <- fread('analysis/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat1, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # for 1907 vs. 2011. sample sizes
dat2 <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat2, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # for 1907 vs. 2014. sample sizes
dat3 <- fread('analysis/Frequency_table_Lof11_Lof14.txt', header=TRUE); setnames(dat3, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # for 2011 vs. 2014

# make sure rows line up
all(dat1$CHROM == dat2$CHROM) 
all(dat1$POS == dat2$POS) 
all(dat1$CHROM == dat3$CHROM) 
all(dat1$POS == dat3$POS) 

# calculate MAF
dat3[,maf1:=Freq_1]
dat3[maf1>0.5, maf1:=1-maf1]
dat3[,maf2:=Freq_2]
dat3[maf2>0.5, maf2:=1-maf2]

# find loci
inds07to11or14 <- dat1[,(Freq_1 <= 0.05 | Freq_1 >= 0.95) & ABS_DIFF>=0.25] & dat2[, (Freq_1 <= 0.05 | Freq_1 >= 0.95) & ABS_DIFF>=0.25] # fall in over-represented categories for 1907-2011 and 1907-2014
	sum(inds07to11or14) # 297, including LG18 POS 18099619
	
inds11to14 <- dat3[,((Freq_1 <= 0.05 | Freq_1 > 0.95) & ABS_DIFF>=0.05) | (((Freq_1 <= 0.1 & Freq_1 > 0.05) | (Freq_1 > 0.9 & Freq_1 <= 0.95)) & ABS_DIFF>=0.15) | (((Freq_1 <= 0.15 & Freq_1 > 0.1) | (Freq_1 > 0.85 & Freq_1 <= 0.9)) & ABS_DIFF>=0.3)] # fall in over-represented categories for 2011-2014 (doesn't get all of them, but gets the major ones)
	sum(inds11to14) # 19089

	sum(inds07to11or14 & !inds11to14) # still 297: no overlap

dat3[inds07to11or14, summary(maf1)] # allele freq in 2011 for candidate outliers
dat3[inds07to11or14, summary(ABS_DIFF)] # change from 2011 to 2014