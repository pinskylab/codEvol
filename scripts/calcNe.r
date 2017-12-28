# Calculate Ne from frequency change at many loci

####################
# load functions
####################
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

################################
# read in data 
################################
## Frequency change and 25 kmer filter: choose ONE
	# 1907-2011
dat <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); nm='1907-2011'; gen=11 # for 1907 vs. 2011.
loc25 <- fread('data_2017.11.24/Norway_25K_mer_positions.txt')

	# 1907-2014
dat <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); nm='1907-2014'; gen=11 # for 1907 vs. 2014
loc25 <- fread('data_2017.11.24/Norway_25K_mer_positions.txt')

	# Canada
dat <- fread('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); nm='Can'; gen=8 # for Canada
loc25 <- fread('data_2017.11.24/Canada_25K_mer_positions.txt')


##############
# Trim loci
##############

# trim to loci that meet 25kmer filter
setkey(dat, CHROM, POS)
setkey(loc25, CHROM, POS)
	nrow(dat)
dat <- dat[loc25, nomatch=0] # nomatch=0 so that non-matching rows are dropped
	nrow(dat)

# trim out inversions and Unplaced
dat <- dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

	dim(dat)

######################
# Run calculations
######################

summary(dat)

### calculate Fs' (Foll et al. 2015)
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


dat[,mean(Fsprime, na.rm=TRUE)]
dat[,1/2/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (diploid individuals)
dat[,1/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (# chromosomes: for comparison to wfabc_1 output)


### Jorde & Ryman/NeEstimator approach
	dat[,z2:=((1-maf1)+(1-maf2))/2] # z for the 2nd allele
	dat[,FsJRnum := (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2] # sum across the alleles for the numerator of Fs
	dat[,FsJRdenom := z*(1-z) + z2*(1-z2)] # sum across the alleles for denominator of Fs
	dat[,sl := 2/(2/N_CHR_1 + 2/N_CHR_2)] # harmonic mean sample size for each locus, in # individuals

	FsJR <- dat[,sum(FsJRnum)]/dat[,sum(FsJRdenom)] # from NeEstimator calculations manual
	S <- nrow(dat)/sum(1/dat[,2/(2/N_CHR_1 + 2/N_CHR_2)]) # harmonic mean sample size in # individuals, across loci and across both times. don't need to worry about weights since all 2 alleles.
	
	S <- 1/mean(1/dat[,sl]) # harmonic mean sample size in # individuals, across loci. don't need to worry about weights since all 2 alleles.
	S2 <- 1/mean(dat[,2/N_CHR_2]) # harmonic mean sample size of 2nd sample in # individuals, across loci. don't need to worry about weights since all 2 alleles.
	FsJRprime <- (FsJR*(1-1/(4*S))-1/S)/((1+FsJR/4)*(1-1/(2*S2)))
	gen/2/FsJRprime # calculation of Ne in # diploid individuals
	gen/FsJRprime # calculation of Ne in # chromosomes
