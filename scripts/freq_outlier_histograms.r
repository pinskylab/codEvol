## Plot histograms of allele frequency change and of p-values
## Compare across different data sets

# load functions
require(data.table)

# read in data
dat <- fread("gunzip -c analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz")
	dat


#################################
# plot histograms of frequency change
#################################
require(RColorBrewer)
cols <- brewer.pal(7, 'Set1')
bks <- seq(0,1,by=0.01)

# histogram
hist14_allk <- dat[, hist(abs(Freq_14-Freq_07), breaks=bks, plot=FALSE)] # all loci
hist11_allk <- dat[, hist(abs(Freq_11-Freq_07), breaks=bks, plot=FALSE)]
hist11_14_allk <- dat[, hist(abs(Freq_14-Freq_11), breaks=bks, plot=FALSE)]
histCan_allk <- dat[, hist(abs(Freq_CanMod-Freq_Can40), breaks=bks, plot=FALSE)]

hist14_25k <- dat[kmer25==1, hist(abs(Freq_14-Freq_07), breaks=bks, plot=FALSE)] # 25kmer filter
hist11_25k <- dat[kmer25==1, hist(abs(Freq_11-Freq_07), breaks=bks, plot=FALSE)]
hist11_14_25k <- dat[kmer25==1, hist(abs(Freq_14-Freq_11), breaks=bks, plot=FALSE)]
histCan_25k <- dat[kmer25==1, hist(abs(Freq_CanMod-Freq_Can40), breaks=bks, plot=FALSE)]

hist14_25k_dp <- dat[kmer25==1 & dpLofFlag==1, hist(abs(Freq_14-Freq_07), breaks=bks, plot=FALSE)] # 25kmer and depth filter
hist11_25k_dp <- dat[kmer25==1 & dpLofFlag==1, hist(abs(Freq_11-Freq_07), breaks=bks, plot=FALSE)]
hist11_14_25k_dp <- dat[kmer25==1 & dpLofFlag==1, hist(abs(Freq_14-Freq_11), breaks=bks, plot=FALSE)]
histCan_25k_dp <- dat[kmer25==1 & dpCanFlag==1, hist(abs(Freq_CanMod-Freq_Can40), breaks=bks, plot=FALSE)]

hist14_25k_dp_chr <- dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), hist(abs(Freq_14-Freq_07), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter
hist11_25k_dp_chr <- dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), hist(abs(Freq_11-Freq_07), breaks=bks, plot=FALSE)]
hist11_14_25k_dp_chr <- dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), hist(abs(Freq_14-Freq_11), breaks=bks, plot=FALSE)]
histCan_25k_dp_chr <- dat[kmer25==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), hist(abs(Freq_CanMod-Freq_Can40), breaks=bks, plot=FALSE)]

hist14_n25k <- dat[kmer25==0, hist(abs(Freq_14-Freq_07), breaks=bks, plot=FALSE)] # fail 25kmer filter
hist11_n25k <- dat[kmer25==0, hist(abs(Freq_11-Freq_07), breaks=bks, plot=FALSE)]
hist11_14_n25k <- dat[kmer25==0, hist(abs(Freq_14-Freq_11), breaks=bks, plot=FALSE)]
histCan_n25k <- dat[kmer25==0, hist(abs(Freq_CanMod-Freq_Can40), breaks=bks, plot=FALSE)]

hist14_ndp <- dat[dpLofFlag==0, hist(abs(Freq_14-Freq_07), breaks=bks, plot=FALSE)] # fail depth filter
hist11_ndp <- dat[dpLofFlag==0, hist(abs(Freq_11-Freq_07), breaks=bks, plot=FALSE)]
hist11_14_ndp <- dat[dpLofFlag==0, hist(abs(Freq_14-Freq_11), breaks=bks, plot=FALSE)]
histCan_ndp <- dat[dpLofFlag==0, hist(abs(Freq_CanMod-Freq_Can40), breaks=bks, plot=FALSE)]

hist14_nchr <- dat[CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'), hist(abs(Freq_14-Freq_07), breaks=bks, plot=FALSE)] # fail chromosome filter
hist11_nchr <- dat[CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'), hist(abs(Freq_11-Freq_07), breaks=bks, plot=FALSE)]
hist11_14_nchr <- dat[CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'), hist(abs(Freq_14-Freq_11), breaks=bks, plot=FALSE)]
histCan_nchr <- dat[CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'), hist(abs(Freq_CanMod-Freq_Can40), breaks=bks, plot=FALSE)]

quartz(width=6, height=6)
# pdf(width=6, height=6, file='figures/abs_diff_hist.pdf')
par(mfrow=c(2,2), las=1, cex.axis=0.8, tcl=-0.2, mgp=c(2,0.5, 0), mai=c(0.5, 0.5, 0.2, 0.1))

with(hist11_allk, plot(mids, density/sum(density), type='o', col=cols[1], cex=0.5, xlab='', ylab='Proportion of loci', log='y', main='1907-2011', ylim=c(1e-6, 1e-1), xlim=c(0,1)))
with(hist11_25k, lines(mids, density/sum(density), type='o', col=cols[2], cex=0.5))
with(hist11_25k_dp, lines(mids, density/sum(density), type='o', col=cols[3], cex=0.5))
with(hist11_25k_dp_chr, lines(mids, density/sum(density), type='o', col=cols[4], cex=0.5))
with(hist11_n25k, lines(mids, density/sum(density), type='o', col=cols[5], cex=0.5))
with(hist11_ndp, lines(mids, density/sum(density), type='o', col=cols[6], cex=0.5))
with(hist11_nchr, lines(mids, density/sum(density), type='o', col=cols[7], cex=0.5))
#with(hist11_14_allk, lines(mids, density/sum(density), type='o', col='grey', cex=0.5))
with(hist11_14_25k_dp_chr, lines(mids, density/sum(density), type='o', col='grey', cex=0.5))
legend('topright', legend=c('all', '25kmer', '25k, depth', '25k, depth, chrom', 'fail 25k', 'fail depth', 'fail chrom', '2011-2014 pass'), lwd=1, pch=1, col=c(cols, 'grey'), cex=0.7, bty='n')

with(hist14_allk, plot(mids, density/sum(density), type='o', col=cols[1], cex=0.5, xlab='', ylab='', log='y', main='1907-2014', ylim=c(1e-6, 1e-1), xlim=c(0,1)))
with(hist14_25k, lines(mids, density/sum(density), type='o', col=cols[2], cex=0.5))
with(hist14_25k_dp, lines(mids, density/sum(density), type='o', col=cols[3], cex=0.5))
with(hist14_25k_dp_chr, lines(mids, density/sum(density), type='o', col=cols[4], cex=0.5))
with(hist14_n25k, lines(mids, density/sum(density), type='o', col=cols[5], cex=0.5))
with(hist14_ndp, lines(mids, density/sum(density), type='o', col=cols[6], cex=0.5))
with(hist14_nchr, lines(mids, density/sum(density), type='o', col=cols[7], cex=0.5))
with(hist11_14_25k_dp_chr, lines(mids, density/sum(density), type='o', col='grey', cex=0.5))

with(histCan_allk, plot(mids, density/sum(density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='Proportion of loci', log='y', main='Canada', ylim=c(1e-6, 1e-1), xlim=c(0,1)))
with(histCan_25k, lines(mids, density/sum(density), type='o', col=cols[2], cex=0.5))
with(histCan_25k_dp, lines(mids, density/sum(density), type='o', col=cols[3], cex=0.5))
with(histCan_25k_dp_chr, lines(mids, density/sum(density), type='o', col=cols[4], cex=0.5))
with(histCan_n25k, lines(mids, density/sum(density), type='o', col=cols[5], cex=0.5))
with(histCan_ndp, lines(mids, density/sum(density), type='o', col=cols[6], cex=0.5))
with(histCan_nchr, lines(mids, density/sum(density), type='o', col=cols[7], cex=0.5))

with(hist11_14_allk, plot(mids, density/sum(density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='', log='y', main='1911-2014', ylim=c(1e-6, 1e-1), xlim=c(0,1)))
with(hist11_14_25k, lines(mids, density/sum(density), type='o', col=cols[2], cex=0.5))
with(hist11_14_25k_dp, lines(mids, density/sum(density), type='o', col=cols[3], cex=0.5))
with(hist11_14_25k_dp_chr, lines(mids, density/sum(density), type='o', col=cols[4], cex=0.5))
with(hist11_14_n25k, lines(mids, density/sum(density), type='o', col=cols[5], cex=0.5))
with(hist11_14_ndp, lines(mids, density/sum(density), type='o', col=cols[6], cex=0.5))
with(hist11_14_nchr, lines(mids, density/sum(density), type='o', col=cols[7], cex=0.5))


dev.off()



#################################
# plot histograms of single-population p-values
#################################
require(RColorBrewer)
cols <- brewer.pal(8, 'RdYlBu')
bks <- seq(0,10,by=0.25)
pch=16

# histogram
histCan <- dat[, hist(-log10(pCan), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter
histCan_25k <- dat[kmer25==1, hist(-log10(pCan), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter
histCan_25k_dp <- dat[kmer25==1 & dpCanFlag==1, hist(-log10(pCan), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter
histCan_25k_dp_chr <- dat[kmer25==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),  hist(-log10(pCan), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter

histLof <- dat[, hist(-log10(pLof), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter
histLof_25k <- dat[kmer25==1, hist(-log10(pLof), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter
histLof_25k_dp <- dat[kmer25==1 & dpLofFlag==1, hist(-log10(pLof), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter
histLof_25k_dp_chr <- dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),  hist(-log10(pLof), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter


#hist_null <- hist(-log10(runif(1000000,0,1)), breaks=bks, plot=FALSE) # null model: uniform (resampling approach)
null <- data.frame(mids=seq(0.125,10,by=0.25), density=NA) # null model: exact calculation of uniform expectation with log10 bins set by bks
	for(i in 1:nrow(null)){
		null$density[i] <- (10^(-bks[i]) - 10^(-bks[i+1]))/(bks[i+1]-bks[i])
	}

quartz(width=5, height=5)
# pdf(width=5, height=5, file='figures/pLof_pCan_hist.pdf')
par(las=1, cex.axis=0.8, tcl=-0.2, mgp=c(2.5,0.5, 0), mai=c(0.75, 0.75, 0.2, 0.1))

with(histCan_25k_dp_chr, plot(mids, density/sum(density), type='o', col=cols[1], cex=0.5, xlab='-log10(p)', ylab='Proportion of loci', log='y', main='Null model p-values', ylim=c(1e-6, 1e-0), xlim=c(0,10), pch=pch))
with(histLof_25k_dp_chr, lines(mids, density/sum(density), type='o', col=cols[8], cex=0.5, pch=pch))

with(histCan_25k_dp, lines(mids, density/sum(density), type='o', col=cols[2], cex=0.5, pch=pch))
with(histLof_25k_dp, lines(mids, density/sum(density), type='o', col=cols[7], cex=0.5, pch=pch))

with(histCan_25k, lines(mids, density/sum(density), type='o', col=cols[3], cex=0.5, pch=pch))
with(histLof_25k, lines(mids, density/sum(density), type='o', col=cols[6], cex=0.5, pch=pch))

with(histCan, lines(mids, density/sum(density), type='o', col=cols[4], cex=0.5, pch=pch))
with(histLof, lines(mids, density/sum(density), type='o', col=cols[5], cex=0.5, pch=pch))

with(null, lines(mids, density/sum(density), type='o', col='light grey', cex=0.5, pch=4))

legend('topright', legend=c('Can pass k,dp,chr', 'Can pass k,dp', 'Can pass k', 'Can all', 'null', 'Lof pass k,dp,chr', 'Lof pass k,dp', 'Lof pass k', 'Lof all'), lwd=1, pch=c(rep(16,4),4, rep(16,4)), col=c(cols[1:4], 'light grey', cols[8:5]), cex=0.7, bty='n', ncol=2)


dev.off()


#################################
# plot histograms of combined p-values
#################################
require(RColorBrewer)
cols <- brewer.pal(10, 'Paired')
bks <- seq(0,10,by=0.25)

# histogram
hist_allk <- dat[, hist(-log10(p.comb), breaks=bks, plot=FALSE)] # all loci

hist_25k <- dat[kmer25==1, hist(-log10(p.comb), breaks=bks, plot=FALSE)] # 25kmer filter

hist_25k_dp <- dat[kmer25==1 & dpLofFlag==1, hist(-log10(p.comb), breaks=bks, plot=FALSE)] # 25kmer and depth filter

hist_25k_dp_chr <- dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),  hist(-log10(p.comb), breaks=bks, plot=FALSE)] # 25kmer, depth, and chromosome filter

hist_n25k <- dat[kmer25==0, hist(-log10(p.comb), breaks=bks, plot=FALSE)] # fail 25kmer filter

hist_ndp <- dat[dpLofFlag==0, hist(-log10(p.comb), breaks=bks, plot=FALSE)] # fail depth filter

hist_nchr <- dat[CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'), hist(-log10(p.comb), breaks=bks, plot=FALSE)] # fail chromosome filter
hist_unpl <- dat[CHROM %in% c('Unplaced'), hist(-log10(p.comb), breaks=bks, plot=FALSE)] # unplaced
hist_inv <- dat[CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12'), hist(-log10(p.comb), breaks=bks, plot=FALSE)] # inversions

hist_25k_dp_inv <- dat[kmer25==1 & dpLofFlag==1 &CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12'), hist(-log10(p.comb), breaks=bks, plot=FALSE)] # inversions but pass kmer and depth filters

#hist_null <- hist(-log10(runif(1000000,0,1)), breaks=bks, plot=FALSE) # null model: uniform (resampling approach)
null <- data.frame(mids=seq(0.125,10,by=0.25), density=NA) # null model: exact calculation of uniform expectation with log10 bins set by bks
	for(i in 1:nrow(null)){
		null$density[i] <- (10^(-bks[i]) - 10^(-bks[i+1]))/(bks[i+1]-bks[i])
	}

quartz(width=5, height=5)
# pdf(width=6, height=6, file='figures/p.comb_hist.pdf')
par(las=1, cex.axis=0.8, tcl=-0.2, mgp=c(2.5,0.5, 0), mai=c(0.75, 0.75, 0.2, 0.1))

with(hist_allk, plot(mids, density/sum(density), type='o', col=cols[1], cex=0.5, xlab='-log10(p)', ylab='Proportion of loci', log='y', main='Null model p-values', ylim=c(1e-6, 1e-0), xlim=c(0,10)))
with(hist_25k, lines(mids, density/sum(density), type='o', col=cols[2], cex=0.5))
with(hist_25k_dp, lines(mids, density/sum(density), type='o', col=cols[3], cex=0.5))
with(hist_25k_dp_chr, lines(mids, density/sum(density), type='o', col=cols[4], cex=0.5))
with(hist_n25k, lines(mids, density/sum(density), type='o', col=cols[5], cex=0.5))
with(hist_ndp, lines(mids, density/sum(density), type='o', col=cols[6], cex=0.5))
with(hist_nchr, lines(mids, density/sum(density), type='o', col=cols[7], cex=0.5))
with(hist_unpl, lines(mids, density/sum(density), type='o', col=cols[8], cex=0.5))
with(hist_inv, lines(mids, density/sum(density), type='o', col=cols[9], cex=0.5))
with(hist_25k_dp_inv, lines(mids, density/sum(density), type='o', col=cols[10], cex=0.5))
with(null, lines(mids, density/sum(density), type='o', col='light grey', cex=0.5, pch=4))
legend('topright', legend=c('all', 'pass kmer', 'pass k&depth', 'pass k,dp,chr', 'fail k', 'fail dp', 'fail chr', 'unplaced', 'inversions', 'inversions pass k,dp', 'null'), lwd=1, pch=c(rep(1,10),4), col=c(cols, 'light grey'), cex=0.7, bty='n')

dev.off()