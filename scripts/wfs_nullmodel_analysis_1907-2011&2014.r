# Compare 1907-2011 and 1907-2014 outlier loci from null model analysis

locnms <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_07-14_25k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

	# make sure it matches
	nrow(dat)
	nrow(locnms)

	# merge
	dat <- as.data.table(dat)
	dat14 <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_14)], by='locusnum')
	setnames(dat14, 'p', 'p14')
	setnames(dat14, 'pmax', 'pmax14')
	setnames(dat14, 'p.adj', 'p.adj14')
	setnames(dat14, 'cnt1', 'cnt07')
	setnames(dat14, 'cnt2', 'cnt14')

locnms <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof11_25k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_07-11_25k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

	# make sure it matches
	nrow(dat)
	nrow(locnms)

	# merge
	dat <- as.data.table(dat)
	dat11 <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_11, ABS_DIFF)], by='locusnum')
	setnames(dat11, 'p', 'p11')
	setnames(dat11, 'pmax', 'pmax11')
	setnames(dat11, 'p.adj', 'p.adj11')
	setnames(dat11, 'cnt1', 'cnt07')
	setnames(dat11, 'cnt2', 'cnt11')
	
# merge
dat <- dat11[dat14, .(CHROM, POS, cnt07, cnt11, cnt14, Freq_07, Freq_11, Freq_14, p11, p14, pmax11, pmax14, p.adj11, p.adj14)]
nrow(dat)
nrow(dat11)
nrow(dat14)

# combine p-values
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE) # https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/

dat[,p.comb := pchisq(-2 * sum(log(c(pmax11, pmax14))),df=2*2,lower=FALSE), by=1:nrow(dat)]

# examine
hist(dat[,-log10(p.comb)], breaks=40, col='grey')


# FDR-correct
dat[,p.comb.adj := p.adjust(dat$p.comb, method='fdr')]


# examine
hist(dat[,-log10(p.comb.adj)], breaks=40, col='grey')

dat[, sum(p.comb.adj < 0.01)]
dat[, sum(p.comb.adj < 0.0001)]

dat[p.comb.adj<0.0001,]