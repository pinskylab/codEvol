# Calculate Ne from frequency change at many loci
# start from ANGSD output
# set up to run on Saga

####################
# load functions
####################
require(data.table)
require(boot) # for bootstrap CIs

################################
# read in data 
################################
# MAFS
dat07 <- fread('data_2020.01.31/Lof_07_freq.mafs.gz', header=TRUE) # 1907 Lof
dat11 <- fread('data_2020.01.31/Lof_11_freq.mafs.gz', header=TRUE) # 2011 Lof
dat14 <- fread('data_2020.01.31/Lof_14_freq.mafs.gz', header=TRUE) # 2014 Lof
datCan40 <- fread('data_2020.01.31/Can_40_freq.mafs.gz', header=TRUE) # 1940 Can
datCan14 <- fread('data_2020.01.31/Can_14_freq.mafs.gz', header=TRUE) # 2014 Can

# high quality loci
#gatk <- fread('data_2020.01.31/GATK_filtered_SNP_set.tab')
gatk <- fread('data_2020.01.31/GATK_filtered_SNP_set_no_Dam.tab')
setnames(gatk, c('CHROM', 'POS'), c('chromo', 'position'))

##############
# Trim loci
##############

# trim to loci that meet GATK filter
setkey(dat07, chromo, position)
setkey(dat11, chromo, position)
setkey(dat14, chromo, position)
setkey(datCan40, chromo, position)
setkey(datCan14, chromo, position)
setkey(gatk, chromo, position)

nrow(dat07)
nrow(dat11)
nrow(dat14)
nrow(datCan40)
nrow(datCan14)

dat07 <- dat07[gatk, nomatch=0] # nomatch=0 so that non-matching rows are dropped
dat11 <- dat11[gatk, nomatch=0] # nomatch=0 so that non-matching rows are dropped
dat14 <- dat14[gatk, nomatch=0] # nomatch=0 so that non-matching rows are dropped
datCan40 <- datCan40[gatk, nomatch=0] # nomatch=0 so that non-matching rows are dropped
datCan14 <- datCan14[gatk, nomatch=0] # nomatch=0 so that non-matching rows are dropped

nrow(dat07)
nrow(dat11)
nrow(dat14)
nrow(datCan40)
nrow(datCan14)

# trim out inversions and Unplaced
dat07 <- dat07[!(chromo %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
dat11 <- dat11[!(chromo %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
dat14 <- dat14[!(chromo %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
datCan40 <- datCan40[!(chromo %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
datCan14 <- datCan14[!(chromo %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

nrow(dat07)
nrow(dat11)
nrow(dat14)
nrow(datCan40)
nrow(datCan14)


################################
# Merge population comparisons
# also trim to loci genotyped in both populations and >0.1
################################
setnames(dat07, c("knownEM", 'nInd'), c("freq1", 'nInd1'))
setnames(dat11, c("knownEM", 'nInd'), c("freq2", 'nInd2'))
setnames(dat14, c("knownEM", 'nInd'), c("freq2", 'nInd2'))
setnames(datCan40, c("knownEM", 'nInd'), c("freq1", 'nInd1'))
setnames(datCan14, c("knownEM", 'nInd'), c("freq2", 'nInd2'))

dat0711 <- dat07[dat11, .(chromo, position, freq1, freq2, nInd1, nInd2)][!is.na(freq1) & !is.na(freq2) & freq1 > 0.1 & freq2 > 0.1, ]
dat0714 <- dat07[dat14, .(chromo, position, freq1, freq2, nInd1, nInd2)][!is.na(freq1) & !is.na(freq2) & freq1 > 0.1 & freq2 > 0.1, ]
datCan <- datCan40[datCan14, .(chromo, position, freq1, freq2, nInd1, nInd2)][!is.na(freq1) & !is.na(freq2) & freq1 > 0.1 & freq2 > 0.1, ]

nrow(dat0711) # 101481
nrow(dat0714) # 99231
nrow(datCan) # 106905

######################
# Run calculations
######################

### calculate Fs' (Foll et al. 2015)

# maf1 and maf2 are allele freqs at each time point
# n1 and n2 are sample sizes
# gen: number of generations between
fsprime <- function(maf1, maf2, n1, n2, gen){
	swtch <- maf1 > 0.5 # which alleles should convert to other allele freq (so < 0.5)
	maf1[swtch] <- 1-maf1[swtch] # convert from >0.5 to <0.5
	maf2[swtch] <- 1-maf2[swtch] # also convert maf2
	z <- (maf1+maf2)/2
	enye <- 2/(1/n1 + 1/n2) # harmonic mean sample size (enye in Foll et al. 2014, s_l in NeEstimator)
	Fs <- (maf1-maf2)^2/(z*(1-z)) # for the minor allele
	return(1/gen * (Fs*(1-1/(2*enye))-2/enye)/((1+Fs/4)*(1-1/n2))) # Fsprime
}

dat0711[, Fsprime := fsprime(freq1, freq2, nInd1, nInd2, 11)]
dat0714[, Fsprime := fsprime(freq1, freq2, nInd1, nInd2, 11)]
datCan[, Fsprime := fsprime(freq1, freq2, nInd1, nInd2, 8)]


dat0711[,mean(Fsprime, na.rm=TRUE)]
dat0714[,mean(Fsprime, na.rm=TRUE)]
datCan[,mean(Fsprime, na.rm=TRUE)]

dat0711[,1/2/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (diploid individuals)
dat0714[,1/2/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (diploid individuals)
datCan[,1/2/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (diploid individuals)

dat0711[,1/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (# chromosomes: for comparison to wfabc_1 output)
dat0714[,1/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (# chromosomes: for comparison to wfabc_1 output)
datCan[,1/mean(Fsprime, na.rm=TRUE)] # calculation of Ne (# chromosomes: for comparison to wfabc_1 output)


### Jorde & Ryman/NeEstimator approach
# Jorde & Ryman 2007

# Ne in # diploid individuals
# as I wrote it first in 2017
jrNe1 <- function(maf1, maf2, n1, n2, gen){
	FsJRnum <- (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2 # the numerator

	z <- (maf1+maf2)/2 # for the first allele
	z2 <- ((1-maf1)+(1-maf2))/2 # z for the 2nd allele
	FsJRdenom <- z*(1-z) + z2*(1-z2) # the denominator of Fs

	sl <- 2/(2/n1 + 2/n2) # harmonic mean sample size for each locus, in # individuals

	FsJR <- sum(FsJRnum)/sum(FsJRdenom) # from NeEstimator calculations manual eq 4.9
	S <- length(FsJRnum)/sum(1/(2/(2/n1 + 2/n2))) # harmonic mean sample size in # individuals, across loci and across both times. don't need to worry about weights since all 2 alleles.
	S <- 1/mean(1/sl) # harmonic mean sample size in # individuals, across loci. don't need to worry about weights since all 2 alleles.
	S2 <- 1/mean(2/n2) # harmonic mean sample size of 2nd sample in # individuals, across loci. don't need to worry about weights since all 2 alleles.
	FsJRprime <- (FsJR*(1 - 1/(4*S)) - 1/S)/((1 + FsJR/4)*(1 - 1/(2*S2)))
	return(gen/2/FsJRprime) # calculation of Ne in # diploid individuals
}

# based on NeEstimator manual v2.1 
# written in 2020
jrNe2 <- function(maf1, maf2, n1, n2, gen){
	Fsnum <- (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2 # the numerator, summing across the two alleles

	z <- (maf1+maf2)/2 # for the first allele
	z2 <- ((1-maf1)+(1-maf2))/2 # z for the 2nd allele
	Fsdenom <- z*(1-z) + z2*(1-z2) # the denominator of Fs, summing across the 2 alleles
	Fs <- sum(Fsnum)/sum(Fsdenom) # from NeEstimator calculations manual

	sl <- 2/(1/n1 + 1/n2) # harmonic mean sample size for each locus, in # individuals

	S <- length(maf1)*2/sum(2/sl) # harmonic mean sample size in # individuals, across loci and across both times. 2 alleles. Eq. 4.10 in NeEstimator v2.1 manual
	S2 <- length(maf2)*2/sum(2/n2) # harmonic mean sample size of 2nd sample in # individuals, across loci. all 2 alleles. See NeEstimator v2.1 below Eq. 4.13
	Fsprime <- (Fs*(1 - 1/(4*S)) - 1/S)/((1 + Fs/4)*(1 - 1/(2*S2))) # Eq. 4.13 in NeEstimator v2.1
	return(gen/(2*Fsprime)) # calculation of Ne in # diploid individuals
}


# old JR estimator
dat0711[, jrNe1(freq1, freq2, nInd1, nInd2, 11)]
dat0714[, jrNe1(freq1, freq2, nInd1, nInd2, 11)]
datCan[, jrNe1(freq1, freq2, nInd1, nInd2, 8)]

# new JR estimator
dat0711[, jrNe2(freq1, freq2, nInd1, nInd2, 11)]
dat0714[, jrNe2(freq1, freq2, nInd1, nInd2, 11)]
datCan[, jrNe2(freq1, freq2, nInd1, nInd2, 8)]



####################################
# bootstrap over loci to get CIs
####################################
# Jorde & Ryman Ne estimator, for boot() to use
jrNe2boot <- function(data, gen, indices){
	maf1 <- data$freq1[indices]
	maf2 <- data$freq2[indices]
	n1 <- data$nInd1[indices]
	n2 <- data$nInd2[indices]

	Fsnum <- (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2 # the numerator, summing across the two alleles

	z <- (maf1+maf2)/2 # for the first allele
	z2 <- ((1-maf1)+(1-maf2))/2 # z for the 2nd allele
	Fsdenom <- z*(1-z) + z2*(1-z2) # the denominator of Fs, summing across the 2 alleles
	Fs <- sum(Fsnum)/sum(Fsdenom) # from NeEstimator calculations manual

	sl <- 2/(1/n1 + 1/n2) # harmonic mean sample size for each locus, in # individuals

	S <- length(maf1)*2/sum(2/sl) # harmonic mean sample size in # individuals, across loci and across both times. 2 alleles. Eq. 4.10 in NeEstimator v2.1 manual
	S2 <- length(maf2)*2/sum(2/n2) # harmonic mean sample size of 2nd sample in # individuals, across loci. all 2 alleles. See NeEstimator v2.1 below Eq. 4.13
	Fsprime <- (Fs*(1 - 1/(4*S)) - 1/S)/((1 + Fs/4)*(1 - 1/(2*S2))) # Eq. 4.13 in NeEstimator v2.1
	Ne <- gen/(2*Fsprime)
	if(Ne < 0) Ne <- Inf
	return(Ne) # calculation of Ne in # diploid individuals
}

# for block bootstrapping across LGs
# frin https://stackoverflow.com/questions/11919808/block-bootstrap-from-subject-list
jrNe2block <- function(lgs, gen, alldata, indices){
	mydata <- do.call("rbind", lapply(indices, function(n) subset(alldata, chromo==lgs[n])))
	return(jrNe2boot(mydata, gen, indices = 1:nrow(mydata)))
}

# regular bootstrap calculations
# not enough memory to do BCa CIs
boot0711 <- boot(data = dat0711, statistic = jrNe2boot, R = 1000, gen = 11)
boot.ci(boot0711, type = 'perc')

boot0714 <- boot(data = dat0714, statistic = jrNe2boot, R = 1000, gen = 11)
boot.ci(boot0714, type = 'perc')

bootCan <- boot(data = datCan, statistic = jrNe2boot, R = 1000, gen = 8)
boot.ci(bootCan, type = 'perc')


# block bootstrapping across LGs
lgs <- dat0711[, sort(unique(chromo))]

boot0711lg <- boot(lgs, jrNe2block, 1000, gen = 11, alldata = dat0711)
boot.ci(boot0711lg, type = c('norm', 'basic', 'perc'))

boot0714lg <- boot(lgs, jrNe2block, 1000, gen = 11, alldata = dat0714)
boot.ci(boot0714lg, type = c('norm', 'basic', 'perc'))

bootCanlg <- boot(lgs, jrNe2block, 1000, gen = 8, alldata = datCan)
boot.ci(bootCanlg, type = c('norm', 'basic', 'perc'))

