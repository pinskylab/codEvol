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
dat14 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat14, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
dat11 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat11, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
datCAN <- fread('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', header=TRUE); setnames(datCAN, 3:7, c('N_CHR_Can40', 'Freq_Can40', 'N_CHR_CanMod', 'Freq_CanMod', 'ABS_DIFF_Can')) # for Canada

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
cols = c('black', 'blue', 'red')
#dat14[,plot(POSgen, ABS_DIFF, pch=16, cex=0.3)]
quartz(height=4, width=8)
# png(height=4, width=8, units='in', res=300, file='figures/abs_diff_vs_pos_NEA_CAN_raw.png')
dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='l', lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[1])]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=0.3, col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_Can, lwd=0.3, col=cols[3])]

col='black'
lgs <- sort(unique(dat14[,CHROM]))
for(j in 1:length(lgs)){
	rng <- range(dat14[CHROM==lgs[j], POSgen])
	if(j %% 2 == 0) lines(x=rng/1e6, y=c(0,0), col=col, lwd=2, lend=2)
	if(j %% 2 == 1) lines(x=rng/1e6, y=c(0.02,0.02), col=col, lwd=2, lend=2)
	text(x=mean(rng/1e6), y=0.03, labels=lgs[j], col=col, cex=0.5)
	
}

#legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011'), lwd=1, col=cols[1:2], bty='n', cex=0.5)
legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada'), lwd=1, col=cols, bty='n', cex=0.5)


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
wndw <- 5e4
lg1 <- 'LG04'; xlims1 <- c(11755415,11755733); nm1<-'LG04 11755415-733'

cols = c('black', 'blue', 'red')

quartz(height=4, width=8)
# png(height=4, width=8, units='in', res=300, file='figures/abs_diff_vs_pos_outlierszoom_NEA&CAN.png')
par(mfrow=c(1,1), mai=c(0.7, 0.7, 0.2, 0.1), las=1, mgp=c(2.5, 1,0))

lg <- lg1; xlims <- xlims1
dat14_25k[CHROM==lg,plot(POS/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), cex=0.5, col=cols[1], main=nm1, xlim=(xlims+c(-wndw, wndw))/1e6)]
dat11_25k[CHROM==lg,lines(POS/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN_25k[CHROM==lg,lines(POS/1e6, ABS_DIFF_Can, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14_150k[CHROM==lg,lines(POS/1e6, ABS_DIFF_0714, lwd=1,cex=0.5,  type='o', col=cols[1], lty=2)]
dat11_150k[CHROM==lg,lines(POS/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2], lty=2)]
datCAN_150k[CHROM==lg,lines(POS/1e6, ABS_DIFF_Can, lwd=1,cex=0.5,  type='o', col=cols[3], lty=2)]


legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada', '25k', '150k'), lwd=1, col=c(cols, 'black', 'black'), lty=c(1,1,1,1,2), bty='n', cex=0.6)


dev.off()
