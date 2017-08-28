# Compare Canada, 1907-2011, 1907-2014, and 2011-2014 outlier loci from null model analysis (25k or 150k)
# to be run on my mac

require(data.table)

# read in data (choose 25kmer or 150kmer)
	##########
	## 25k
	kmer='25'
	# 1907-2014
locnms <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_07-14_25k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		dat14 <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_14)], by='locusnum')
		setnames(dat14, c('p', 'pmax', 'p.adj', 'cnt1', 'cnt2'), c('p14', 'pmax14', 'p.adj14', 'cnt07', 'cnt14'))

	# 1907-2011 25k
locnms <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof11_25k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_07-11_25k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		dat11 <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_11)], by='locusnum')
		setnames(dat11, c('p', 'pmax', 'p.adj', 'cnt1', 'cnt2'), c('p11', 'pmax11', 'p.adj11', 'cnt07', 'cnt11'))

	# 2011-2014 25k
locnms <- fread('data/data_29.06.17/Frequency_table_Lof11_Lof14_25k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_11', 'Freq_11', 'N_CHR_14', 'Freq_14', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_11-14_25k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		dat1114 <- merge(dat, locnms[,.(locusnum, CHROM, POS)], by='locusnum')
		setnames(dat1114, c('p', 'pmax', 'p.adj', 'cnt1', 'cnt2'), c('p1114', 'pmax1114', 'p.adj1114', 'cnt11', 'cnt14'))

	# Canada 25k
locnms <- fread('data/data_11.07.17/Frequency_table_Can_40_Can_25k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_Can40', 'Freq_Can40', 'N_CHR_CanMod', 'Freq_CanMod', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_Can_25k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		datCan <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_Can40, Freq_CanMod)], by='locusnum')
		setnames(datCan, c('p', 'pmax', 'p.adj', 'cnt1', 'cnt2'), c('pCan', 'pmaxCan', 'p.adjCan', 'cntCan40', 'cntCanMod'))


	##########
	## 150k
	kmer='150'
	# 1907-2014 150k
locnms <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF0714_150k')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_07-14_150k.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

	# make sure it matches
	nrow(dat)
	nrow(locnms)

	# merge
	dat <- as.data.table(dat)
	dat14 <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_14, ABS_DIFF0714_150k)], by='locusnum')
	setnames(datC14, c('p', 'pmax', 'p.adj', 'cnt1', 'cnt2'), c('p14', 'pmax14', 'p.adj14', 'cnt07', 'cnt14'))

	# 1907-2011 150k
	
	# 2011-2014 150k
	
	# Canada 150k

##################
# merge
##################

dat <- dat11[dat14, .(CHROM, POS, cnt07, cnt11, cnt14, Freq_07, Freq_11, Freq_14, p11, p14, pmax11, pmax14, p.adj11, p.adj14)]
	nrow(dat)
	nrow(dat11)
	nrow(dat14)

setkey(dat, CHROM, POS)
setkey(dat1114, CHROM, POS)
dat <- dat[dat1114, .(CHROM, POS, cnt07, cnt11, cnt14, Freq_07, Freq_11, Freq_14, p11, p14, p1114, pmax11, pmax14, pmax1114, p.adj11, p.adj14, p.adj1114)]
	nrow(dat)
	nrow(dat1114)

setkey(dat, CHROM, POS)
setkey(datCan, CHROM, POS)
dat <- datCan[dat, .(CHROM, POS, cnt07, cnt11, cnt14, cntCan40, cntCanMod, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, p11, p14, p1114, pCan, pmax11, pmax14, pmax1114, pmaxCan, p.adj11, p.adj14, p.adj1114, p.adjCan)]
	nrow(dat)
	nrow(dat1114)
	nrow(datCan)


###############################################
# combine p-values for individual loci
# https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/
###############################################

dat[,p.comb11_14 := pchisq(-2 * sum(log(c(pmax11, pmax14))),df=2*2,lower=FALSE), by=1:nrow(dat)]
dat[,p.comb11_Can := pchisq(-2 * sum(log(c(pmax11, pmaxCan))),df=2*2,lower=FALSE), by=1:nrow(dat)]
dat[,p.comb14_Can := pchisq(-2 * sum(log(c(pmax14, pmaxCan))),df=2*2,lower=FALSE), by=1:nrow(dat)]

# FDR-correct
dat[,p.comb11_14.adj := p.adjust(p.comb11_14, method='fdr')]
dat[,p.comb11_Can.adj := p.adjust(p.comb11_Can, method='fdr')]
dat[,p.comb14_Can.adj := p.adjust(p.comb14_Can, method='fdr')]


#########################
# examine
#########################
# examine combined p-values
hist(dat[,-log10(p.comb11_14)], breaks=40, col='grey')
hist(dat[,-log10(p.comb11_Can)], breaks=40, col='grey')
hist(dat[,-log10(p.comb14_Can)], breaks=40, col='grey')

# examine histograms for adjusted combined p-values
hist(dat[,-log10(p.comb11_14.adj)], breaks=40, col='grey')
hist(dat[,-log10(p.comb11_Can.adj)], breaks=40, col='grey')
hist(dat[,-log10(p.comb14_Can.adj)], breaks=40, col='grey')

	# zoom
hist(dat[,-log10(p.comb11_14.adj)], breaks=200, col='grey', xlim=c(2,6), ylim=c(0,50))
hist(dat[,-log10(p.comb11_Can.adj)], breaks=200, col='grey', xlim=c(2,6), ylim=c(0,50))
hist(dat[,-log10(p.comb14_Can.adj)], breaks=200, col='grey', xlim=c(2,6), ylim=c(0,50))

# how many
	# basic: p<X
dat[, sum(p11 < 0.05, na.rm=TRUE)]
dat[, sum(p14 < 0.05, na.rm=TRUE)]
dat[, sum(pCan < 0.05, na.rm=TRUE)]

dat[, sum(p11 < 0.0005, na.rm=TRUE)]
dat[, sum(p14 < 0.0005, na.rm=TRUE)]
dat[, sum(pCan < 0.0005, na.rm=TRUE)]

	# p.comb.adj < X
dat[, sum(p.comb11_14.adj < 0.01, na.rm=TRUE)]
dat[, sum(p.comb11_14.adj < 0.0001, na.rm=TRUE)]

dat[, sum(p.comb11_Can.adj < 0.01, na.rm=TRUE)]
dat[, sum(p.comb11_Can.adj < 0.001, na.rm=TRUE)]
dat[, sum(p.comb11_Can.adj < 0.0001, na.rm=TRUE)]

dat[, sum(p.comb14_Can.adj < 0.01, na.rm=TRUE)]
dat[, sum(p.comb14_Can.adj < 0.001, na.rm=TRUE)]
dat[, sum(p.comb14_Can.adj < 0.0001, na.rm=TRUE)]

# examine specific outliers
dat[p.comb11_14.adj<0.001,]
dat[p.comb11_Can.adj<0.001,]
dat[p.comb14_Can.adj<0.001,]


# mark outliers based on fdr-corrected p-values (not distance to other loci)
# don't include the inversions or unplaced
	# 1907-2011
dat[,outlier07_11 := 0]
dat[p.adj11<0.1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_11 := 1]

	# 1907-2014
dat[,outlier07_14 := 0]
dat[p.adj14<0.1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_14 := 1]

	# 1907-2014
dat[,outlierCan := 0]
dat[p.adjCan<0.1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlierCan := 1]

	# combined 1907-2011 and 1907-2014
dat[,outlier07_11_14 := 0]
dat[p.comb11_14.adj<0.001 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_11_14 := 1]

	# combined 1907-2011 and Canada
dat[,outlier07_11_Can := 0]
dat[p.comb11_Can.adj<0.001 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_11_Can := 1]

	# combined 1907-2014 and Canada
dat[,outlier07_14_Can := 0]
dat[p.comb14_Can.adj<0.001 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_14_Can := 1]



# examine outliers
	# 1907-2011
dat[outlier07_11 == 1,] # only 3

	# 1907-2014
print(dat[outlier07_14 == 1,], nrow=dat[,sum(outlier07_14 == 1)]) # LG04 locus and LG08 mapping problem region

	# Canada
dat[outlierCan == 1,] # only 4

	# combined 1907-2011 and 1907-2014
dat[outlier07_11_14 == 1,]

	# combined 1907-2011 and Canada
dat[outlier07_11_Can == 1,]

	# combined 1907-2014 and Canada
dat[outlier07_14_Can == 1,]

# write out tab-separated for Bastiaan
outfile <- paste('analysis/wfs_nullmodel_outliers_07-11-14_Can_', kmer, 'k.tsv.gz', sep='')
outfile
write.table(dat, file=gzfile(outfile), sep='\t', row.names=FALSE, quote=FALSE)