# Compare Canada, 1907-2011, 1907-2014, and 2011-2014 outlier loci from null model analysis (25k)
# to be run on my mac

require(data.table)

# read in data
	# 1907-2011-2014
locnms <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
locnms2 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms2, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	setkey(locnms, CHROM, POS, N_CHR_07, Freq_07)
	setkey(locnms2, CHROM, POS, N_CHR_07, Freq_07)
	locnms <- locnms[locnms2,]
load('analysis/wfs_nullmodel_pvals_07-11-14.rdata') # dat. has the p-values.
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
# merge Lof and Can
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
kmer='25' # flag for file names

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

#########################
# add depth statistic
#########################
dpLof <- fread('analysis/Outlier_sequencing_depth_statistic_All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz.txt')
dpCan <- fread('analysis/Outlier_sequencing_depth_statistic_All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz.txt')
	setnames(dpLof, c('CHROM', 'POS', 'dpstatLof'))
	setnames(dpCan, c('CHROM', 'POS', 'dpstatCan'))
depthstat <- merge(dpLof, dpCan, all=TRUE)
	nrow(dpLof)
	nrow(dpCan)
	nrow(depthstat)

dat <- merge(dat, depthstat, by=c('CHROM', 'POS'), all.x=TRUE)
	nrow(dat)

dat[,dpLofFlag := dpstatLof < quantile(dpstatLof, na.rm=TRUE, probs=0.95)]
	dat[is.na(dpLofFlag), dpLofFlag := FALSE]
dat[,dpCanFlag := dpstatCan < quantile(dpstatCan, na.rm=TRUE, probs=0.95)]
	dat[is.na(dpCanFlag), dpCanFlag := FALSE]

	dat[,sum(dpLofFlag)]
	dat[,sum(dpLofFlag)]/dat[,sum(!is.na(dpstatLof))] # should be 0.95
	dat[,sum(dpCanFlag)]
	dat[,sum(dpCanFlag)]/dat[,sum(!is.na(dpstatCan))]

dat[,dpFlag := dpLofFlag & dpCanFlag]
	dat[,sum(dpFlag)]

###############################################
# combine p-values for individual loci
# https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/
###############################################

dat[,p.comb := pchisq(-2 * sum(log(c(pLof, pCan))),df=2*2,lower=FALSE), by=1:nrow(dat)]

# FDR-correct
dat[,p.comb.adj := p.adjust(p.comb, method='fdr')]
dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),p.comb.adj2 := p.adjust(p.comb, method='fdr')] # after masking out inversions and unplaced
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),p.comb.adj3 := p.adjust(p.comb, method='fdr')] # after masking out inversions and unplaced

# FDR-correct each population separately
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),p.Lof.adj3 := p.adjust(pLof, method='fdr')] # after masking out inversions and unplaced
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),p.Can.adj3 := p.adjust(pCan, method='fdr')] # after masking out inversions and unplaced


#########################
# examine
#########################
# number of loci
nrow(dat) # total
dat[,sum(!is.na(cnt07))] # Lof
dat[,sum(!is.na(cntCan40))] # Can
dat[,sum(!is.na(cnt07) & !is.na(cntCan40))] # both

nrow(dat[kmer25==1,])
dat[kmer25==1,sum(!is.na(cnt07))] # Lof
dat[kmer25==1,sum(!is.na(cntCan40))] # Can
dat[kmer25==1,sum(!is.na(cnt07) & !is.na(cntCan40))] # both

dat[kmer25==1 & dpLofFlag==1,sum(!is.na(cnt07))] # Lof
dat[kmer25==1 & dpCanFlag==1,sum(!is.na(cntCan40))] # Can
dat[kmer25==1 & dpFlag==1,sum(!is.na(cnt07) & !is.na(cntCan40))] # both

dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(!is.na(cnt07) & !is.na(cntCan40))] # both


# do the same SNPs appear as outliers in both?
ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan),.N]
nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan),sum(pLof<1e-3)]
nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pLof),sum(pCan<1e-3)]
ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(pCan<1e-3 & pLof<1e-3)]
nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
nLof
nCan
ncomb # number in both
nexp

binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot)

# examine combined p-values in various classes
dat[,range(p.comb, na.rm=TRUE)]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb, na.rm=TRUE)]

hist(dat[,-log10(p.comb)], breaks=40, col='grey')

	require(RColorBrewer)
	cols <- brewer.pal(3, 'Set1')
	bks <- seq(0,1,by=0.01)
	hist14_25k <- dat14_25k[,hist(ABS_DIFF_0714, breaks=bks, plot=FALSE)]
	hist11_25k <- dat11_25k[,hist(ABS_DIFF_0711, breaks=bks, plot=FALSE)]
	hist11_14_25k <- dat11_14_25k[,hist(ABS_DIFF_1114, breaks=bks, plot=FALSE)]

	quartz(width=9, height=3)
	# pdf(width=9, height=3, file='analysis/figures/abs_diff_hist_25k&150k.pdf')
	par(mfrow=c(1,3), las=1, cex.axis=0.8, tcl=-0.2, mgp=c(2.5,0.5, 0), mai=c(0.5, 0.5, 0.2, 0.1))

	plot(hist11_25k$mids, hist11_25k$density/sum(hist11_25k$density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='Proportion of loci', log='y', main='1907-2011', ylim=c(1e-6, 1e-1), xlim=c(0,0.6))
	lines(hist11_150k$mids, hist11_150k$density/sum(hist11_150k$density), type='o', col=cols[2], cex=0.5, main='1907-2011')

	plot(hist14_25k$mids, hist14_25k$density/sum(hist14_25k$density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='', log='y', main='1907-2014', ylim=c(1e-6, 1e-1), xlim=c(0,0.6))
	lines(hist14_150k$mids, hist14_150k$density/sum(hist14_150k$density), type='o', col=cols[2], cex=0.5)

	plot(hist11_14_25k$mids, hist11_14_25k$density/sum(hist11_14_25k$density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='', log='y', main='2011-2014', ylim=c(1e-6, 1e-1), xlim=c(0,0.6))
	lines(hist11_14_150k$mids, hist11_14_150k$density/sum(hist11_14_150k$density), type='o', col=cols[2], cex=0.5)

	legend('topright', legend=c('25 kmer', '150 kmer'), lwd=1, pch=1, col=cols[1:2], cex=0.8)


# examine adjusted combined p-values
dat[,range(p.comb.adj, na.rm=TRUE)]
dat[,range(p.comb.adj2, na.rm=TRUE)]

dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb.adj, na.rm=TRUE)]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb.adj2, na.rm=TRUE)]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb.adj3, na.rm=TRUE)]

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
dat[, sum(p.comb.adj2 < 0.1, na.rm=TRUE)]
dat[, sum(p.comb.adj2 < 0.01, na.rm=TRUE)]
dat[, sum(p.comb.adj2 < 0.001, na.rm=TRUE)]
dat[, sum(p.comb.adj2 < 0.0001, na.rm=TRUE)]
dat[, sum(p.comb.adj2 < 0.00001, na.rm=TRUE)]

dat[kmer25==1, sum(p.comb.adj2 < 0.2, na.rm=TRUE)]
dat[kmer25==1, sum(p.comb.adj2 < 0.01, na.rm=TRUE)]
dat[kmer25==1, sum(p.comb.adj2 < 0.0001, na.rm=TRUE)]

dat[kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), sum(p.comb.adj2 < 0.2, na.rm=TRUE)]

	# p.comb.adj3 < X (filter on kmer25, depth, chromosome)
dat[, sum(p.comb.adj3 < 0.2, na.rm=TRUE)]
dat[, sum(p.comb.adj3 < 0.1, na.rm=TRUE)]
dat[, sum(p.comb.adj3 < 0.01, na.rm=TRUE)]
dat[, sum(p.comb.adj3 < 0.001, na.rm=TRUE)]

# examine specific outliers
dat[p.comb.adj3<0.2,]


# mark outliers based on fdr-corrected p-values (not distance to other loci)
# don't include the inversions or unplaced or those outside kmer25
	# combined 1907-2011 and 1907-2014
dat[,outlier07_11_14_Can := 0]
dat[p.comb.adj3<0.2, outlier07_11_14_Can := 1]
	dat[,sum(outlier07_11_14_Can, na.rm=TRUE)]


# plot outliers based on allele frequency change in both populations
par(mfrow=c(1,2))
dat[kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),][sample(.N,5000), plot(abs((Freq_11 +Freq_14)/2 - Freq_07), abs(Freq_CanMod - Freq_Can40), main='all that pass filters\n(sample of 5000)', xlim=c(0,0.5), ylim=c(0,0.5))]
dat[p.comb.adj3<0.2 & kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), plot(abs((Freq_11 +Freq_14)/2 - Freq_07), abs(Freq_CanMod - Freq_Can40), main='outliers p.comb.adj2<0.05', xlim=c(0,0.5), ylim=c(0,0.5))]


# write out tab-separated for Bastiaan
outfile <- paste('analysis/wfs_nullmodel_outliers_07-11-14_Can_', kmer, 'k.tsv.gz', sep='')
outfile
write.table(dat, file=gzfile(outfile), sep='\t', row.names=FALSE, quote=FALSE)