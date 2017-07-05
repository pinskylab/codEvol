# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){ # if not on cod node
	require(data.table)
	require(plyr)
	require(parallel)
	ncores=3
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){ # if on a cod node
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	ncores=20
}

# read in data (choose one)
	# 1907-2011
dat <- fread('data_29.06.17/Frequency_table_Lof07_Lof11_25k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=50; c2=46; nm='1907-2011'; gen=11 # for 1907 vs. 2011. sample sizes. Use 25 kmer since most conservative (least likely biased by read mapping errors)

dat <- fread('data_29.06.17/Frequency_table_Lof07_Lof11_150k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=50; c2=46; nm='1907-2011'; gen=11 # for 1907 vs. 2011. sample sizes.


	# 1907-2014
dat <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=50; c2=48; nm='1907-2014'; gen=11 # for 1907 vs. 2014. sample sizes. Use 25 kmer since most conservative (least likely biased by read mapping errors)

dat <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_50k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=50; c2=48; nm='1907-2014'; gen=11 # for 1907 vs. 2014. sample sizes

dat <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_100k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=50; c2=48; nm='1907-2014'; gen=11 # for 1907 vs. 2014. sample sizes

dat <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=50; c2=48; nm='1907-2014'; gen=11 # for 1907 vs. 2014. sample sizes

summary(dat)

# calculate Fs' (Foll et al. 2015)
gen=11
dat[,maf1:=Freq_1] # minor allele freq for first time point
dat[,maf2:=Freq_2]
	switch <- dat[,maf1>0.5] # which alleles should convert to other allele freq (so < 0.5)
	dat[switch,maf1:=1-maf1] # convert from >0.5 to <0.5
	dat[switch,maf2:=1-maf2] # also convert maf2
dat[,z:=(maf1+maf2)/2]
dat[,enye:=2/(1/N_CHR_1 + 1/N_CHR_2)] # harmonic mean sample size (enye in Foll et al. 2014, s_l in NeEstimator)
dat[,Fs:=(maf1-maf2)^2/(z*(1-z))] # for the minor allele
dat[,Fsprime := 1/gen * (Fs*(1-1/(2*enye))-2/enye)/((1+Fs/4)*(1-1/N_CHR_2))]

	# understand NA values
	dat[,summary(Fsprime)]
	summary(dat[is.na(Fsprime),]) # maf1 and maf2 are 0 where Fsprime is NA

	# choose a set of loci
inds <- !(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) # not inversions or unplaced (use this)
#inds <- !(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12'))
#inds <- !(dat$CHROM %in% c('Unplaced')) # not Unplaced
#inds <- (dat$CHROM %in% c('Unplaced')) # only Unplaced
#inds <- (dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')) # only inversions
#inds <- rep(TRUE, nrow(dat)) # all loci
sum(inds)

dat[inds,mean(Fsprime, na.rm=TRUE)]
dat[inds,1/2/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (diploid individuals)
dat[inds,1/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (# chromosomes: for comparison to wfabc_1 output)

	# Jorde & Ryman/NeEstimator approach
	dat[,z2:=((1-maf1)+(1-maf2))/2] # z for the 2nd allele
	dat[,FsJRnum := (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2] # sum across the alleles for the numerator of Fs
	dat[,FsJRdenom := z*(1-z) + z2*(1-z2)] # sum across the alleles for denominator of Fs
	dat[,sl := 2/(2/N_CHR_1 + 2/N_CHR_2)] # harmonic mean sample size for each locus, in # individuals

	FsJR <- dat[inds,sum(FsJRnum)]/dat[inds,sum(FsJRdenom)] # from NeEstimator calculations manual
	S <- sum(inds)/sum(1/dat[inds,2/(2/N_CHR_1 + 2/N_CHR_2)]) # harmonic mean sample size in # individuals, across loci and across both times. don't need to worry about weights since all 2 alleles.
	
	S <- 1/mean(1/dat[inds,sl]) # harmonic mean sample size in # individuals, across loci. don't need to worry about weights since all 2 alleles.
	S2 <- 1/mean(dat[inds,2/N_CHR_2]) # harmonic mean sample size of 2nd sample in # individuals, across loci. don't need to worry about weights since all 2 alleles.
	FsJRprime <- (FsJR*(1-1/(4*S))-1/S)/((1+FsJR/4)*(1-1/(2*S2)))
	gen/2/FsJRprime # calculation of Ne in # diploid individuals
	gen/FsJRprime # calculation of Ne in # chromosomes
