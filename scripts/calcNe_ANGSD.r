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
gatk <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab')
setnames(gatk, c('chromo', 'position', 'ref', 'alt'))

# unlinked loci
unlnkCan <- fread('analysis/ld.unlinked.Can.gatk.nodam.csv.gz', drop=1)
unlnkLof <- fread('analysis/ld.unlinked.Lof.gatk.nodam.csv.gz', drop=1)
setnames(unlnkCan, c('chromo', 'position', 'cluster', 'nloci'))
setnames(unlnkLof, c('chromo', 'position', 'cluster', 'nloci'))

##############
# Trim loci
##############

# trim to loci that meet GATK and unlinked filter
setkey(dat07, chromo, position)
setkey(dat11, chromo, position)
setkey(dat14, chromo, position)
setkey(datCan40, chromo, position)
setkey(datCan14, chromo, position)
setkey(gatk, chromo, position)
setkey(unlnkCan, chromo, position)
setkey(unlnkLof, chromo, position)

nrow(dat07)
nrow(dat11)
nrow(dat14)
nrow(datCan40)
nrow(datCan14)

dat07 <- dat07[gatk, nomatch=0] # nomatch=0 so that non-matching rows are dropped
dat11 <- dat11[gatk, nomatch=0]
dat14 <- dat14[gatk, nomatch=0]
datCan40 <- datCan40[gatk, nomatch=0]
datCan14 <- datCan14[gatk, nomatch=0]

nrow(dat07)
nrow(dat11)
nrow(dat14)
nrow(datCan40)
nrow(datCan14)

dat07 <- dat07[unlnkLof, nomatch=0] # nomatch=0 so that non-matching rows are dropped
dat11 <- dat11[unlnkLof, nomatch=0]
dat14 <- dat14[unlnkLof, nomatch=0]
datCan40 <- datCan40[unlnkCan, nomatch=0]
datCan14 <- datCan14[unlnkCan, nomatch=0]

nrow(dat07)
nrow(dat11)
nrow(dat14)
nrow(datCan40)
nrow(datCan14)

# trim out Unplaced
dat07 <- dat07[!(chromo == 'Unplaced'),]
dat11 <- dat11[!(chromo == 'Unplaced'),]
dat14 <- dat14[!(chromo == 'Unplaced'),]
datCan40 <- datCan40[!(chromo == 'Unplaced'),]
datCan14 <- datCan14[!(chromo == 'Unplaced'),]

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

nrow(dat0711) # 94976 (w/out unplaced), 95426 (w/ unplaced)
nrow(dat0714) # 92730, 93154
nrow(datCan) # 93728, 94263

######################
# Run calculations
######################

### Jorde & Ryman/NeEstimator approach
# Jorde & Ryman 2007

# Ne in # diploid individuals
# based on NeEstimator manual v2.1 
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
# boot0711 <- boot(data = dat0711, statistic = jrNe2boot, R = 1000, gen = 11)
# boot.ci(boot0711, type = 'perc')
# 
# boot0714 <- boot(data = dat0714, statistic = jrNe2boot, R = 1000, gen = 11)
# boot.ci(boot0714, type = 'perc')
# 
# bootCan <- boot(data = datCan, statistic = jrNe2boot, R = 1000, gen = 8)
# boot.ci(bootCan, type = 'perc')


# block bootstrapping across LGs
lgs <- dat0711[, sort(unique(chromo))]

boot0711lg <- boot(lgs, jrNe2block, 4000, gen = 11, alldata = dat0711)
print(boot0711lg)
median(boot0711lg$t[is.finite(boot0711lg$t)]) # median bootstrap
boot.ci(boot0711lg, type = c('perc'))

boot0714lg <- boot(lgs, jrNe2block, 10000, gen = 11, alldata = dat0714)
print(boot0714lg)
median(boot0714lg$t[is.finite(boot0714lg$t)]) # median bootstrap
boot.ci(boot0714lg, type = c('perc'))

bootCanlg <- boot(lgs, jrNe2block, 4000, gen = 8, alldata = datCan)
print(bootCanlg)
median(bootCanlg$t[is.finite(bootCanlg$t)]) # median bootstrap
boot.ci(bootCanlg, type = c('perc'))
