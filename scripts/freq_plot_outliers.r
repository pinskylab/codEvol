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
dat14 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat14, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
dat11 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat11, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
datCAN <- fread('data_2019_03_18/Frequency_table_CAN_40_TGA.txt', header=TRUE); setnames(datCAN, 3:7, c('N_CHR_Can40', 'Freq_Can40', 'N_CHR_CanMod', 'Freq_CanMod', 'ABS_DIFF_Can')) # for Canada

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- dat14[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, CHROM)

# merge nucleotide position into the frequency files
setkey(dat14, CHROM)
dat14 <- dat14[chrmax[,.(CHROM, start)], ]
dat14[,POSgen:=POS+start]

setkey(dat11, CHROM)
dat11 <- dat11[chrmax[,.(CHROM, start)], ]
dat11[,POSgen:=POS+start]

setkey(datCAN, CHROM)
datCAN <- datCAN[chrmax[,.(CHROM, start)], ]
datCAN[,POSgen:=POS+start]

# combine the datasets to look at correlations across them
setkey(dat11, CHROM, POS)
setkey(dat14, CHROM, POS)
setkey(datCAN, CHROM, POS)

#dat11_14 <- dat11[dat14, .(CHROM, POS, POSgen, N_CHR_07, N_CHR_11, N_CHR_14, Freq_07, Freq_11, Freq_14, ABS_DIFF_0711, ABS_DIFF_0714)]

datCAN_14 <- datCAN[dat14, .(CHROM, POS, POSgen, N_CHR_Can40, N_CHR_CanMod, N_CHR_07, N_CHR_14, Freq_Can40, Freq_CanMod, Freq_07, Freq_14, ABS_DIFF_Can, ABS_DIFF_0714)]

datCAN_11_14 <- datCAN_14[dat11, .(CHROM, POS, POSgen, N_CHR_Can40, N_CHR_CanMod, N_CHR_07, N_CHR_11, N_CHR_14, Freq_Can40, Freq_CanMod, Freq_07, Freq_11, Freq_14, ABS_DIFF_Can, ABS_DIFF_0711, ABS_DIFF_0714)]

##########################
# plot frequency change
# no smoothing
##########################
addchroms <- function(){
	col='black'
	lgs <- sort(unique(dat14[,CHROM]))
	for(j in 1:length(lgs)){
		rng <- range(dat14[CHROM==lgs[j], POSgen])
		if(j %% 2 == 0) lines(x=rng/1e6, y=c(0,0), col=col, lwd=2, lend=2)
		if(j %% 2 == 1) lines(x=rng/1e6, y=c(0.02,0.02), col=col, lwd=2, lend=2)
		text(x=mean(rng/1e6), y=0.03, labels=lgs[j], col=col, cex=0.5)
	
	}
}

cols = c('black', 'blue', 'red')
#dat14[,plot(POSgen, ABS_DIFF, pch=16, cex=0.3)]
quartz(height=12, width=8)
# png(height=12, width=8, units='in', res=300, file='figures/abs_diff_vs_pos_NEA_CAN_raw.png')
par(mfrow=c(4,1), mai=c(0.5, 1, 0.2, 0.5))
dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[1])]
dat11[,points(POSgen/1e6, ABS_DIFF_0711, cex=0.2, lwd=0.3, col=cols[2])]
datCAN[,points(POSgen/1e6, ABS_DIFF_Can, cex=0.2, lwd=0.3, col=cols[3])]

addchroms()
legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada'), lwd=1, col=cols, bty='n', cex=0.5)

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[1])]
addchroms()

dat11[,plot(POSgen/1e6, ABS_DIFF_0711, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[2])]
addchroms()

datCAN[,plot(POSgen/1e6, ABS_DIFF_Can, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[3])]
addchroms()


dev.off()


#############################
# running mean of allele frequency difference 
# (if saved, can skip this)
#############################
#stp = 1e4; windsz='1e4' 
stp = 1e6; windsz='1e6' 

mids <- seq(stp/4, max(dat14$POSgen), by=stp/2)
dat14mean <- data.frame(mids =mids, mean = rep(NA, length(mids))) # for mean of p-adj
nrow(dat14mean)
for(j in 1:nrow(dat14mean)){ # quick for stp=1e6, 10 min for stp=1e4
	if(j %% 100 == 0) cat(paste(j, ' ', sep=''))
	inds <- dat14$POSgen >= (dat14mean$mids[j]-stp/2) & dat14$POSgen < (dat14mean$mids[j]+stp/2)
	dat14mean$mean[j] <- mean(dat14$ABS_DIFF_0714[inds])
}
dat14mean <- as.data.table(dat14mean)

mids <- seq(stp/4, max(dat11$POSgen), by=stp/2)
dat11mean <- data.frame(mids=mids, mean=rep(NA, length(mids))) # for mean of p-adj
nrow(dat11mean)
for(j in 1:nrow(dat11mean)){ # quick for stp1e6. 10+ min for stp=1e4
	if(j %% 100 == 0) cat(paste(j, ' ', sep=''))
	inds <- dat11$POSgen >= (dat11mean$mids[j]-stp/2) & dat11$POSgen < (dat11mean$mids[j]+stp/2)
	dat11mean$mean[j] <- mean(dat11$ABS_DIFF_0711[inds])
}
dat11mean <- as.data.table(dat11mean)


mids <- seq(stp/4, max(datCAN$POSgen), by=stp/2)
datCANmean <- data.frame(mids=mids, mean = rep(NA, length(mids))) # for mean of p-adj
nrow(datCANmean)
for(j in 1:nrow(datCANmean)){ # takes a couple minutes
	if(j %% 100 == 0) cat(paste(j, ' ', sep=''))
	inds <- datCAN$POSgen >= (datCANmean$mids[j]-stp/2) & datCAN$POSgen < (datCANmean$mids[j]+stp/2)
	datCANmean$mean[j] <- mean(datCAN$ABS_DIFF_Can[inds])
}
datCANmean <- as.data.table(datCANmean)


# save
save(dat14mean, file=paste('analysis/Frequency_table_Lof07_Lof14_runmean', windsz, '.rdata', sep=''))
save(dat11mean, file=paste('analysis/Frequency_table_Lof07_Lof11_runmean', windsz, '.rdata', sep=''))
save(datCANmean, file=paste('analysis/Frequency_table_CAN_runmean', windsz, '.rdata', sep=''))

#################################
# load running means back in
#################################
# 1e4
# NOTE: loads dat14mean and dat11mean, without _25k or _150k suffixes
load('analysis/Frequency_table_Lof07_Lof14_runmean1e4.rdata'); load('analysis/Frequency_table_Lof07_Lof11_runmean1e4.rdata'); ylims=c(0,0.5); windsz='1e4'

# 1e6
# NOTE: loads dat14mean and dat11mean, without _25k or _150k suffixes
load('analysis/Frequency_table_Lof07_Lof14_runmean1e6.rdata'); load('analysis/Frequency_table_Lof07_Lof11_runmean1e6.rdata'); ylims=c(0,0.15); windsz='1e6'
load('analysis/Frequency_table_CAN_runmean1e6.rdata')


###########################################
# plot frequency change from running mean
###########################################
require(RColorBrewer)
cols = c('black', 'blue', 'red')
#cols = brewer.pal(4, 'BrBG') # for 25k and 150k, 2011 and 2014
quartz(height=4, width=8)
# png(height=4, width=8, units='in', res=300, file=paste('figures/abs_diff_vs_pos_NEA&CAN_runmean', windsz, '.png', sep=''))
dat14mean[,plot(mids/1e6, mean, type='l', lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=ylims, col=cols[1], main=paste('Running mean', windsz))] # for 1e6 window
dat11mean[,lines(mids/1e6, mean, lwd=0.3, col=cols[2])] # for 1e6 window
datCANmean[,lines(mids/1e6, mean, lwd=0.3, col=cols[3])] # for 1e6 window

col='red'
lgs <- sort(unique(dat14[,CHROM]))
for(j in 1:length(lgs)){
	rng <- range(dat14[CHROM==lgs[j], POSgen])
	if(j %% 2 == 0) lines(x=rng/1e6, y=c(0,0), col=col, lwd=2, lend=2)
	if(j %% 2 == 1) lines(x=rng/1e6, y=c(0.005,0.005), col=col, lwd=2, lend=2)
	text(x=mean(rng/1e6), y=0.01, labels=lgs[j], col=col, cex=0.5)
	
}

#legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011'), lwd=1, col=cols[1:2], bty='n', cex=0.5)
legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada'), lwd=1, col=cols, bty='n', cex=0.5)

dev.off()


#################################################################
# plot frequency difference for specific regions of the genome
#################################################################
#require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node
require(data.table) # not on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
dat <- fread("gzcat analysis/wfs_nullmodel_outliers_07-11-14_Can.csv.gz") # on mac
setkey(dat, CHROM, POS)

# find outliers
locinds <- dat[, which(q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2)]

# plot parameters
wndw <- 5e4
div <- 1e3; unit='kb'
cols = c('black', 'blue', 'red')
colstr = c('#00000055', '#0000FF55', '#FF000055')

# make plot
quartz(height=ceiling(2*length(locinds)/3), width=8)
# png(height=ceiling(2*length(locinds)/3), width=8, units='in', res=600, file='figures/abs_diff_vs_pos_outlierszoom_NEA&CAN.png')
par(mfrow=c(ceiling(length(locinds)/3),3), mai=c(0.3, 0.25, 0.1, 0.05), las=1, mgp=c(1, 0.3,0), tcl=-0.2)

for(i in 1:length(locinds)){
	# find region to plot
	indx <- dat[locinds[i], .(CHROM, start=POS-wndw, end=POS+wndw)]
	dat[CHROM==indx$CHROM & POS > indx$start & POS < indx$end, plot(POS/div, abs(Freq_11-Freq_07), type='p', cex=0.5, xlab=paste('Position (', unit, ')', sep=''), ylab='∆ frequency', ylim=c(0,1), col=cols[1], main=paste(dat[locinds[i], .(CHROM, POS)], collapse=' '), cex.axis=0.5)]
	dat[CHROM==indx$CHROM & POS > indx$start & POS < indx$end, points(POS/div, abs(Freq_14-Freq_07), cex=0.5, col=cols[2])]
	dat[CHROM==indx$CHROM & POS > indx$start & POS < indx$end, points(POS/div, abs(Freq_CanMod-Freq_Can40), cex=0.5, col=cols[3])]

	# loess lines
	lw1 <- dat[CHROM==indx$CHROM & POS > indx$start & POS < indx$end, loess(I(abs(Freq_11-Freq_07)) ~ I(POS/div))]
	lw2 <- dat[CHROM==indx$CHROM & POS > indx$start & POS < indx$end, loess(I(abs(Freq_14-Freq_07)) ~ I(POS/div))]
	lw3 <- dat[CHROM==indx$CHROM & POS > indx$start & POS < indx$end, loess(I(abs(Freq_CanMod-Freq_Can40)) ~ I(POS/div))]
	lines(lw1$x,lw1$fitted,col=cols[1],lwd=1)
	lines(lw2$x,lw2$fitted,col=cols[2],lwd=1)
	lines(lw3$x,lw3$fitted,col=cols[3],lwd=1)

	# legend
	if(i==1) legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada'), pch=1, cex=0.5, col=cols, bty='o')

	# vertical lines for outlier loci
	dat[(q3.Lof0711<0.2 | q3.comb0711Can<0.2) & CHROM==indx$CHROM, abline(v=(POS-wndw/150)/div, lty=2, col=colstr[1])]
	dat[(q3.Lof0714<0.2 | q3.comb0714Can<0.2) & CHROM==indx$CHROM, abline(v=POS/div, lty=2, col=colstr[2])]
	dat[(q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2) & CHROM==indx$CHROM, abline(v=(POS+wndw/150)/div, lty=2, col=colstr[3])]

}



dev.off()
