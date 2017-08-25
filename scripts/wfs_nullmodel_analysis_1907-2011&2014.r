# Compare 1907-2011 and 1907-2014 outlier loci from null model analysis (25k and 150k)
require(data.table)

# read in data
	# 1907-2014 25k
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

	# 1907-2011 25k
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

	# 1907-2014 150k
locnms <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF0714_150k')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_07-14_150k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

	# make sure it matches
	nrow(dat)
	nrow(locnms)

	# merge
	dat <- as.data.table(dat)
	dat14_150k <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_14, ABS_DIFF0714_150k)], by='locusnum')
	setnames(dat14_150k, 'p', 'p14_150k')
	setnames(dat14_150k, 'pmax', 'pmax14_150k')
	setnames(dat14_150k, 'p.adj', 'p.adj14_150k')
	setnames(dat14_150k, 'cnt1', 'cnt07')
	setnames(dat14_150k, 'cnt2', 'cnt14')

	
# merge
dat <- dat11[dat14, .(CHROM, POS, cnt11, Freq_11, p11, p14, pmax11, pmax14, p.adj11, p.adj14)]
	nrow(dat)
	nrow(dat11)
	nrow(dat14)

setkey(dat, CHROM, POS)
setkey(dat14_150k, CHROM, POS)
dat2 <- dat[dat14_150k, .(CHROM, POS, cnt07, cnt11, cnt14, Freq_07, Freq_11, Freq_14, p11, p14, p14_150k, pmax11, pmax14, pmax14_150k, p.adj11, p.adj14, p.adj14_150k)]
	nrow(dat14_150k)
	nrow(dat2)

# combine p-values
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE) # https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/

dat2[,p.comb11_14 := pchisq(-2 * sum(log(c(pmax11, pmax14))),df=2*2,lower=FALSE), by=1:nrow(dat2)]

# examine
hist(dat2[,-log10(p.comb11_14)], breaks=40, col='grey')


# FDR-correct
dat2[,p.comb11_14.adj := p.adjust(p.comb11_14, method='fdr')]


# examine
hist(dat2[,-log10(p.comb11_14.adj)], breaks=40, col='grey')

dat2[, sum(p.comb11_14.adj < 0.01)]
dat2[, sum(p.comb11_14.adj < 0.0001)]

dat2[p.comb11_14.adj<0.0001,]


# mark outliers based on fdr-corrected p-values (not distance to other loci)
# don't include the inversions
	# 1907-2011 25k
dat2[,outlier07_11_25k := 0]
dat2[p.adj11<0.1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_11_25k := 1]

	# 1907-2014 25k
dat2[,outlier07_14_25k := 0]
dat2[p.adj14<0.1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_14_25k := 1]

	# 1907-2014 150k
dat2[,outlier07_14_150k := 0]
dat2[p.adj14_150k<0.05 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_14_150k := 1]

	# combined 1907-2011 and 1907-2014
dat2[,outlier07_11_14_25k := 0]
dat2[p.comb11_14.adj<0.001 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_11_14_25k := 1]


# examine outliers
	# 1907-2011 25k
dat2[outlier07_11_25k == 1,] # only 3

	# 1907-2014 25k
print(dat2[outlier07_14_25k == 1,], nrow=dat2[,sum(outlier07_14_25k == 1)]) # LG08 mapping problem region

	# 1907-2014 150k
print(dat2[outlier07_14_150k == 1,.(CHROM, POS, cnt07, cnt11, cnt14, Freq_07, Freq_11, Freq_14, p14_150k, pmax14_150k, p.adj14_150k, outlier07_14_150k)], nrow=dat2[,sum(outlier07_14_150k == 1)]) # 

	# combined 1907-2011 and 1907-2014
dat2[outlier07_11_14_25k == 1,]


# write out tab-separated for Bastiaan
write.table(dat2, 'analysis/wfs_nullmodel_outliers_07-11&14_25k&150k.tsv', sep='\t', row.names=FALSE, quote=FALSE)