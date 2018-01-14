# Compare Canada, 1907-2011, 1907-2014, and 2011-2014 outlier loci from null model analysis (25k)
# to be run on my mac

require(data.table)

# read in data
	##########
	## 25k
	kmer='25'
	# 1907-2011-2014
locnms <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
locnms2 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms2, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	setkey(locnms, CHROM, POS, N_CHR_07, Freq_07)
	setkey(locnms2, CHROM, POS, N_CHR_07, Freq_07)
	locnms <- locnms[locnms2,]
load('analysis/wfs_nullmodel_pvals_07-11-14.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		datLof <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_11, Freq_14)], by='locusnum')
		setnames(datLof, c('p', 'p.adj', 'cnt1', 'cnt2', 'cnt3'), c('pLof', 'p.adjLof', 'cnt07', 'cnt11', 'cnt14'))

	# Canada
locnms <- fread('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_Can40', 'Freq_Can40', 'N_CHR_CanMod', 'Freq_CanMod', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_Can.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		datCan <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_Can40, Freq_CanMod)], by='locusnum')
		setnames(datCan, c('p', 'p.adj', 'cnt1', 'cnt2'), c('pCan', 'p.adjCan', 'cntCan40', 'cntCanMod'))




##################
# merge
##################

setkey(datLof, CHROM, POS)
setkey(datCan, CHROM, POS)
dat <- merge(datCan, datLof, all=TRUE)[, .(CHROM, POS, cnt07, cnt11, cnt14, cntCan40, cntCanMod, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof, pCan)]
	nrow(dat)
	nrow(datLof)
	nrow(datCan)
	
################
# add kmer25
################
kmer25Lof <- fread('data_2017.11.24/Norway_25K_mer_positions.txt')
kmer25Can <- fread('data_2017.11.24/Canada_25K_mer_positions.txt')

setkey(kmer25Lof, CHROM, POS)
setkey(kmer25Can, CHROM, POS)
kmer25 <- merge(kmer25Lof, kmer25Can, all=TRUE) # merge and keep all
	nrow(kmer25Lof)
	nrow(kmer25Can)
	nrow(kmer25)

kmer25[,kmer25 := 1] # add a flag
dat <- merge(dat, kmer25, by=c('CHROM', 'POS'), all.x=TRUE)
dat[is.na(kmer25), kmer25:=0] # set NAs to 0
	nrow(dat)
	dat[,sum(kmer25)]

###############################################
# combine p-values for individual loci
# https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/
###############################################

dat[,p.comb := pchisq(-2 * sum(log(c(pLof, pCan))),df=2*2,lower=FALSE), by=1:nrow(dat)]

# FDR-correct
dat[,p.comb.adj := p.adjust(p.comb, method='fdr')]
dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),p.comb.adj2 := p.adjust(p.comb, method='fdr')] # after masking out inversions and unplaced


#########################
# examine
#########################
# examine combined p-values
hist(dat[,-log10(p.comb)], breaks=40, col='grey')

# examine histograms for adjusted combined p-values
hist(dat[,-log10(p.comb.adj)], breaks=40, col='grey')
hist(dat[,-log10(p.comb.adj2)], breaks=40, col='grey')

	# zoom
hist(dat[,-log10(p.comb.adj)], breaks=200, col='grey', xlim=c(2,6), ylim=c(0,50))
hist(dat[,-log10(p.comb.adj2)], breaks=200, col='grey', xlim=c(2,6), ylim=c(0,50))

# how many
	# basic: p<X
dat[, sum(pLof < 0.05, na.rm=TRUE)]
dat[, sum(pCan < 0.05, na.rm=TRUE)]

dat[kmer25==1, sum(pLof < 0.05, na.rm=TRUE)]
dat[kmer25==1, sum(pCan < 0.05, na.rm=TRUE)]

dat[, sum(pLof < 0.00001, na.rm=TRUE)]
dat[, sum(pCan < 0.00001, na.rm=TRUE)]

dat[kmer25==1, sum(pLof < 0.00001, na.rm=TRUE)]
dat[kmer25==1, sum(pCan < 0.00001, na.rm=TRUE)]


	# p.comb.adj < X
dat[, sum(p.comb.adj < 0.01, na.rm=TRUE)]
dat[, sum(p.comb.adj < 0.0001, na.rm=TRUE)]

dat[kmer25==1, sum(p.comb.adj < 0.01, na.rm=TRUE)]
dat[kmer25==1, sum(p.comb.adj < 0.0001, na.rm=TRUE)]


	# p.comb.adj2 < X
dat[, sum(p.comb.adj2 < 0.01, na.rm=TRUE)]
dat[, sum(p.comb.adj2 < 0.0001, na.rm=TRUE)]

dat[kmer25==1, sum(p.comb.adj2 < 0.01, na.rm=TRUE)]
dat[kmer25==1, sum(p.comb.adj2 < 0.0001, na.rm=TRUE)]


# examine specific outliers
dat[p.comb.adj<0.00001,]


# mark outliers based on fdr-corrected p-values (not distance to other loci)
# don't include the inversions or unplaced or those outside kmer25
	# combined 1907-2011 and 1907-2014
dat[,outlier07_11_14_Can := 0]
dat[p.comb.adj2<0.05 & kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), outlier07_11_14_Can := 1]
	dat[,sum(outlier07_11_14_Can, na.rm=TRUE)]


# examine outliers
	# combined 1907-2011-2014 and Canada
print(dat[outlier07_11_14_Can == 1,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof, pCan, p.comb.adj2)], nrow=300)

# write out tab-separated for Bastiaan
outfile <- paste('analysis/wfs_nullmodel_outliers_07-11-14_Can_', kmer, 'k.tsv.gz', sep='')
outfile
write.table(dat, file=gzfile(outfile), sep='\t', row.names=FALSE, quote=FALSE)