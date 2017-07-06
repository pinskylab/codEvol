## Plot allele frequency change vs. genome position
## Compare across different data sets

# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){ # not cod node
	require(data.table)
	require(plyr)
	require(parallel)
	ncores=2
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	ncores=20
}

# read in data
dat14_25k <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(dat14_25k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
dat11_25k <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof11_25k.txt', header=TRUE); setnames(dat11_25k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
dat11_14_25k <- fread('data/data_29.06.17/Frequency_table_Lof11_Lof14_25k.txt', header=TRUE); setnames(dat11_14_25k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_1114')) # for 1907 vs. 2011
#datCAN <- fread('analysis/Frequency_table_Can_40_Can_TGA.txt', header=TRUE); setnames(datCAN, 3:7, c('N_CHR_40', 'Freq_40', 'N_CHR_TGA', 'Freq_TGA', 'ABS_DIFF_40TGA')) # for 1907 vs. 2011

dat14_150k <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE); setnames(dat14_150k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
dat11_150k <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof11_150k.txt', header=TRUE); setnames(dat11_150k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
dat11_14_150k <- fread('data/data_29.06.17/Frequency_table_Lof11_Lof14_150k.txt', header=TRUE); setnames(dat11_14_150k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_1114')) # for 1907 vs. 2011


#################################
# plot histograms of frequency change
#################################
require(RColorBrewer)
cols <- brewer.pal(3, 'Set1')
bks <- seq(0,1,by=0.01)
hist14_25k <- dat14_25k[,hist(ABS_DIFF_0714, breaks=bks, plot=FALSE)]
hist11_25k <- dat11_25k[,hist(ABS_DIFF_0711, breaks=bks, plot=FALSE)]
hist11_14_25k <- dat11_14_25k[,hist(ABS_DIFF_1114, breaks=bks, plot=FALSE)]
hist14_150k <- dat14_150k[,hist(ABS_DIFF_0714, breaks=bks, plot=FALSE)]
hist11_150k <- dat11_150k[,hist(ABS_DIFF_0711, breaks=bks, plot=FALSE)]
hist11_14_150k <- dat11_14_150k[,hist(ABS_DIFF_1114, breaks=bks, plot=FALSE)]


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

dev.off()