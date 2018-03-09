## Calculate correlation in allele frequency and frequency change across datasets
## New comments!

# load functions
require(data.table)
require(RColorBrewer)

# read in data on all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
dat <- fread("gunzip -c analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz")
	dat


####################################################
# correlations in allele frequency among datasets
####################################################
# old to new within a site
dat[!is.na(Freq_07) & !is.na(Freq_11),cor.test(Freq_07, Freq_11, method='pearson')]
dat[!is.na(Freq_07) & !is.na(Freq_14),cor.test(Freq_07, Freq_14, method='pearson')]
dat[!is.na(Freq_11) & !is.na(Freq_14),cor.test(Freq_11, Freq_14, method='pearson')]
dat[!is.na(Freq_Can40) & !is.na(Freq_CanMod),cor.test(Freq_Can40, Freq_CanMod, method='pearson')]

# across sites
dat[!is.na(Freq_07) & !is.na(Freq_Can40),cor.test(Freq_07, Freq_Can40, method='pearson')]
dat[!is.na(Freq_11) & !is.na(Freq_CanMod),cor.test(Freq_11, Freq_CanMod, method='pearson')]
dat[!is.na(Freq_14) & !is.na(Freq_CanMod),cor.test(Freq_14, Freq_CanMod, method='pearson')]


###############################################################
# correlations in allele frequency change among populations
###############################################################
## All loci
# Lof 07-11 to Lof 07-14 (expect more correlation)
dat[!is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_14), cor.test(Freq_11-Freq_07, Freq_14 - Freq_07, method='pearson')]

# Lof 07-11 to Can
dat[!is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_Can40) & !is.na(Freq_CanMod), cor.test(Freq_11-Freq_07, Freq_CanMod - Freq_Can40, method='pearson')]

# Lof 07-14 to Can
dat[!is.na(Freq_07) & !is.na(Freq_14) & !is.na(Freq_Can40) & !is.na(Freq_CanMod), cor.test(Freq_14-Freq_07, Freq_CanMod - Freq_Can40, method='pearson')]


## Outliers (expect more correlation than all loci)
# Lof 07-11 to Lof 07-14 (expect more correlation)
dat[!is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_14) & outlierLof_q02==1, cor.test(Freq_11-Freq_07, Freq_14 - Freq_07, method='pearson')]

# Lof 07-11 to Can (pLof or pCan)
dat[!is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_Can40) & !is.na(Freq_CanMod) & (outlierLof_q02==1 | outlierCan_q02==1), cor.test(Freq_11-Freq_07, Freq_CanMod - Freq_Can40, method='pearson')]

# Lof 07-14 to Can (pLof or pCan)
dat[!is.na(Freq_07) & !is.na(Freq_14) & !is.na(Freq_Can40) & !is.na(Freq_CanMod) & (outlierLof_q02==1 | outlierCan_q02==1), cor.test(Freq_14-Freq_07, Freq_CanMod - Freq_Can40, method='pearson')]

# Lof 07-11 to Can (p.comb)
dat[!is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_Can40) & !is.na(Freq_CanMod) & outlierLof_Can_q02==1, cor.test(Freq_11-Freq_07, Freq_CanMod - Freq_Can40, method='pearson')]

# Lof 07-14 to Can (p.comb)
dat[!is.na(Freq_07) & !is.na(Freq_14) & !is.na(Freq_Can40) & !is.na(Freq_CanMod) & outlierLof_Can_q02==1, cor.test(Freq_14-Freq_07, Freq_CanMod - Freq_Can40, method='pearson')]


## Compare to noise (11-14)
# Lof 11-14 to Can (all loci)
dat[!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !is.na(Freq_CanMod), cor.test(Freq_14-Freq_11, Freq_CanMod - Freq_Can40, method='pearson')]

# Lof 11-14 to Can (outliers)
dat[!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !is.na(Freq_CanMod) & outlierLof_Can_q02==1, cor.test(Freq_14-Freq_11, Freq_CanMod - Freq_Can40, method='pearson')]


###############################################################
# plot correlations in allele frequency change among regions
###############################################################

# Plot correlations for all loci and outliers (latter in red)
cols <- brewer.pal(3, 'Set1')
colg <- rgb(0,0,0,0.1)
cexs <- c(0.5, 0.5, 0.6)

quartz(width=8, height=8)
# png(width=8, height=8, res=300, units='in', file='figures/freq_change.png')
par(mfrow=c(2,2))

# CAN and 07-14 
#dat[sample(.N,10000),plot(Freq_14-Freq_07, Freq_CanMod-Freq_Can40, col=colg, pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Can vs. Lof', xlab='Frequency change Lof', ylab='Frequency change Can')] # downsample for pdf or screen
dat[,plot(Freq_14-Freq_07, Freq_CanMod-Freq_Can40, col=colg, pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Can vs. Lof', xlab='Frequency change Lof', ylab='Frequency change Can')]

dat[outlierLof_Can_q02==1,points(Freq_14-Freq_07, Freq_CanMod-Freq_Can40, col=cols[1], pch=16, cex=0.5)]
dat[outlierLof_q02==1,points(Freq_14-Freq_07, Freq_CanMod-Freq_Can40, col=cols[2], pch=1, cex=0.5)]
dat[outlierCan_q02==1,points(Freq_14-Freq_07, Freq_CanMod-Freq_Can40, col=cols[3], pch=1, cex=0.6)]

legend('topleft', legend=c('all', 'combined outliers', 'Lof outliers', 'Can outliers'), pch=c(16,16,1,1), pt.cex=c(0.3,cexs), col=c(colg, cols), bty='n', cex=0.8)

# Lof 07-11 vs 07-14 all loci and outliers
#dat[sample(.N,10000),plot(Freq_11-Freq_07, Freq_14-Freq_07, col=colg, pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Lof 1907-2014 vs. Lof 1907-2011', xlab='Frequency change 1907-2011', ylab='Frequency change 1907-2014')]
dat[,plot(Freq_11-Freq_07, Freq_14-Freq_07, col=colg, pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Lof 1907-2014 vs. Lof 1907-2011', xlab='Frequency change 1907-2011', ylab='Frequency change 1907-2014')]

dat[outlierLof_Can_q02==1,points(Freq_11-Freq_07, Freq_14-Freq_07, col=cols[1], pch=16, cex=cexs[1])]
dat[outlierLof_q02==1,points(Freq_11-Freq_07, Freq_14-Freq_07, col=cols[2], pch=1, cex=cexs[2])]
dat[outlierCan_q02==1,points(Freq_11-Freq_07, Freq_14-Freq_07, col=cols[3], pch=1, cex=cexs[3])]

# Lof 11-14 (sampling noise) vs CAN all loci and outliers
#dat[sample(.N,10000),plot(Freq_14-Freq_11, Freq_CanMod-Freq_Can40, col=rgb(0,0,0,0.3), pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Can vs. sampling noise', xlab='Frequency change Lof 2011-2014', ylab='Frequency change Can')]
dat[,plot(Freq_14-Freq_11, Freq_CanMod-Freq_Can40, col=colg, pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Can vs. sampling noise', xlab='Frequency change Lof 2011-2014', ylab='Frequency change Can')]

dat[outlierLof_Can_q02==1,points(Freq_14-Freq_11, Freq_CanMod-Freq_Can40, col=cols[1], pch=16, cex=cexs[1])]
dat[outlierLof_q02==1,points(Freq_14-Freq_11, Freq_CanMod-Freq_Can40, col=cols[2], pch=1, cex=cexs[2])]
dat[outlierCan_q02==1,points(Freq_14-Freq_11, Freq_CanMod-Freq_Can40, col=cols[3], pch=1, cex=cexs[3])]

# Lof 07-14 vs. 11-14 (sampling noise) all loci and outliers
#dat[sample(.N,10000),plot(Freq_14-Freq_11, Freq_14-Freq_07, col=rgb(0,0,0,0.3), pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Lof vs. sampling noise', xlab='Frequency change Lof 2011-2014', ylab='Frequency change Lof 1907-2014')]
dat[,plot(Freq_14-Freq_11, Freq_14-Freq_07, col=colg, pch=16, cex=0.3, xlim=c(-1,1), ylim=c(-1,1), main='Lof vs. sampling noise', xlab='Frequency change Lof 2011-2014', ylab='Frequency change Lof 1907-2014')]

dat[outlierLof_Can_q02==1,points(Freq_14-Freq_11, Freq_14-Freq_07, col=cols[1], pch=16, cex=cexs[1])]
dat[outlierLof_q02==1,points(Freq_14-Freq_11, Freq_14-Freq_07, col=cols[2], pch=1, cex=cexs[2])]
dat[outlierCan_q02==1,points(Freq_14-Freq_11, Freq_14-Freq_07, col=cols[3], pch=1, cex=cexs[3])]

dev.off()
