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
dat14 <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(dat14, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
dat11 <- fread('analysis/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat11, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
datCAN <- fread('analysis/Frequency_table_Can_40_Can_TGA.txt', header=TRUE); setnames(datCAN, 3:7, c('N_CHR_40', 'Freq_40', 'N_CHR_TGA', 'Freq_TGA', 'ABS_DIFF_40TGA')) # for 1907 vs. 2011

	# 2016-10 data (higher SNP density, but lower quality)
dat14orig <- fread('analysis/Frequency_table_Lof07_Lof14_2016-10.txt', header=TRUE); setnames(dat14orig, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714'))

# make a nucleotide position for the whole genome
#chrmax <- dat14[,.(len=max(POS)), by=CHROM]
chrmax <- dat14orig[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, CHROM)

setkey(dat14, CHROM)
dat14 <- dat14[chrmax[,.(CHROM, start)], ]
dat14[,POSgen:=POS+start]

setkey(dat11, CHROM)
dat11 <- dat11[chrmax[,.(CHROM, start)], ]
dat11[,POSgen:=POS+start]

setkey(datCAN, CHROM)
datCAN <- datCAN[chrmax[,.(CHROM, start)], ]
datCAN[,POSgen:=POS+start]

setkey(dat14orig, CHROM)
dat14orig <- dat14orig[chrmax[,.(CHROM, start)], ]
dat14orig[,POSgen:=POS+start]


# combine the datasets
setkey(dat11, CHROM, POS)
setkey(dat14, CHROM, POS)
setkey(datCAN, CHROM, POS)

dat11_14 <- dat11[dat14, .(CHROM, POS, POSgen, N_CHR_07, N_CHR_11, N_CHR_14, Freq_07, Freq_11, Freq_14, ABS_DIFF_0711, ABS_DIFF_0714)]

datCAN_14 <- datCAN[dat14, .(CHROM, POS, POSgen, N_CHR_40, N_CHR_TGA, N_CHR_07, N_CHR_14, Freq_40, 
Freq_TGA, Freq_07, Freq_14, ABS_DIFF_40TGA, ABS_DIFF_0714)]

datCAN_11_14 <- datCAN_14[dat11, .(CHROM, POS, POSgen, N_CHR_40, N_CHR_TGA, N_CHR_07, N_CHR_11, N_CHR_14, Freq_40, Freq_TGA, Freq_07, Freq_11, Freq_14, ABS_DIFF_40TGA, ABS_DIFF_0711, ABS_DIFF_0714)]

# running mean of allele frequency difference (saved, so can skip this)
#stp = 1e6
#mids <- seq(stp/4, max(dat14$POSgen), by=stp/2)
#dat14mean <- data.frame(mids =mids, mean = rep(NA, length(mids))) # for mean of p-adj
#nrow(dat14mean)
#for(j in 1:nrow(dat14mean)){ # takes a couple minutes
#if(j %% 100 == 0) cat(j)
#inds <- dat14$POSgen >= (dat14mean$mids[j]-stp/2) & dat14$POSgen < (dat14mean$mids[j]+stp/2)
#dat14mean$mean[j] <- mean(dat14$ABS_DIFF_0714[inds])
#}
#dat14mean <- as.data.table(dat14mean)
#
#mids <- seq(stp/4, max(dat11$POSgen), by=stp/2)
#dat11mean <- data.frame(mids=mids, mean=rep(NA, length(mids))) # for mean of p-adj
#nrow(dat11mean)
#for(j in 1:nrow(dat11mean)){ # takes a couple minutes
#if(j %% 100 == 0) cat(j)
#inds <- dat11$POSgen >= (dat11mean$mids[j]-stp/2) & dat11$POSgen < (dat11mean$mids[j]+stp/2)
#dat11mean$mean[j] <- mean(dat11$ABS_DIFF_0711[inds])
#}
#dat11mean <- as.data.table(dat11mean)
#
#mids <- seq(stp/4, max(datCAN$POSgen), by=stp/2)
#datCANmean <- data.frame(mids=mids, mean = rep(NA, length(mids))) # for mean of p-adj
#nrow(datCANmean)
#for(j in 1:nrow(datCANmean)){ # takes a couple minutes
#if(j %% 100 == 0) cat(j)
#inds <- datCAN$POSgen >= (datCANmean$mids[j]-stp/2) & datCAN$POSgen < (datCANmean$mids[j]+stp/2)
#datCANmean$mean[j] <- mean(datCAN$ABS_DIFF_40TGA[inds])
#}
#datCANmean <- as.data.table(datCANmean)
#
#
## save
#save(dat14mean, file='analysis/Frequency_table_Lof07_Lof14_runmean.rdata')
#save(dat11mean, file='analysis/Frequency_table_Lof07_Lof11_runmean.rdata')
#save(datCANmean, file='analysis/Frequency_table_CAN_runmean.rdata')

# load running means back in
load('analysis/Frequency_table_Lof07_Lof14_runmean.rdata')
load('analysis/Frequency_table_Lof07_Lof11_runmean.rdata')
load('analysis/Frequency_table_CAN_runmean.rdata')


# plot frequency difference
cols = c('black', 'blue', 'red')
#dat14[,plot(POSgen, ABS_DIFF, pch=16, cex=0.3)]
quartz(height=4, width=8)
# png(height=4, width=8, units='in', res=300, file='analysis/figures/abs_diff_vs_pos.png')
# png(height=4, width=8, units='in', res=300, file='analysis/figures/abs_diff_vs_pos_NEA&CAN.png')
dat14mean[,plot(mids/1e6, mean, type='l', lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,0.15), col=cols[1])] # for 1e6 window
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

legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada'), lwd=1, col=cols, bty='n', cex=0.5)

dev.off()


# plot frequency difference for specific regions
wndw <- 5e4
xlims1 <- dat14[CHROM=='LG04' & POS>=30790553 & POS <=30790644, range(POSgen)] + c(-wndw, wndw); nm1<-'LG04 30790553'
xlims2 <- dat14[CHROM=='LG08' & POS>=16347278 & POS <=16347801, range(POSgen)] + c(-wndw, wndw); nm2<-'LG08 16347278'
xlims3 <- dat14[CHROM=='LG09' & POS>=8665949 & POS <=8667347, range(POSgen)] + c(-wndw, wndw); nm3<-'LG09 8665949'
xlims4 <- dat14[CHROM=='LG10' & POS>=10275862 & POS <=10275897, range(POSgen)] + c(-wndw, wndw); nm4<-'LG10 10275862'
xlims5 <- dat14[CHROM=='LG11' & POS>=8187239 & POS <=8187315, range(POSgen)] + c(-wndw, wndw); nm5<-'LG11 8187239'
xlims6 <- dat14[CHROM=='LG18' & POS>=18099619 & POS <=18099632, range(POSgen)] + c(-wndw, wndw); nm6<-'LG18 18099619'
xlims7 <- dat14[CHROM=='LG23' & POS>=6727974 & POS <=6728060, range(POSgen)] + c(-wndw, wndw); nm7<-'LG23 6727974'
xlims8 <- dat14[CHROM=='LG23' & POS>=7614593 & POS <=7614595, range(POSgen)] + c(-wndw, wndw); nm8<-'LG23 8665949'
cols = c('black', 'blue', 'red', 'grey')

quartz(height=8, width=8)
# png(height=16, width=8, units='in', res=300, file='analysis/figures/abs_diff_vs_pos_outlierszoom_NEA&CAN.png')
par(mfrow=c(8,1), mai=c(0.5, 0.5, 0.2, 0.1), las=1, mgp=c(2.5, 1,0))

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims1/1e6, cex=0.5, col=cols[1], main=nm1)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]

legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada', 'NEA 1907-2014 lowqual'), lwd=1, col=cols, bty='n', cex=0.6)

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims2/1e6, cex=0.5, col=cols[1], main=nm2)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims3/1e6, cex=0.5, col=cols[1], main=nm3)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims4/1e6, cex=0.5, col=cols[1], main=nm4)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims5/1e6, cex=0.5, col=cols[1], main=nm5)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims6/1e6, cex=0.5, col=cols[1], main=nm6)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims7/1e6, cex=0.5, col=cols[1], main=nm7)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='o', lwd=1, xlab='Position (Mb)', ylab='∆ frequency', ylim=c(0,0.5), xlim=xlims8/1e6, cex=0.5, col=cols[1], main=nm8)]
dat11[,lines(POSgen/1e6, ABS_DIFF_0711, lwd=1, cex=0.5, type='o', col=cols[2])]
datCAN[,lines(POSgen/1e6, ABS_DIFF_40TGA, lwd=1,cex=0.5,  type='o', col=cols[3])]
dat14orig[,lines(POSgen/1e6, ABS_DIFF_0714, lwd=1, cex=0.5, type='o', col=cols[4])]


dev.off()


# plot correlations among regions
	# 07-11 and 07-14
dat11_14[,plot(ABS_DIFF_0711, ABS_DIFF_0714, col=rgb(0,0,0,0.5), pch=16, cex=0.3)]
	dat11_14[,cor.test(ABS_DIFF_0711, ABS_DIFF_0714)]

	# CAN and 07-14
datCAN_14[,plot(ABS_DIFF_40TGA, ABS_DIFF_0714, col=rgb(0,0,0,0.5), pch=16, cex=0.3)]
	datCAN_14[,cor.test(ABS_DIFF_40TGA, ABS_DIFF_0714)]

	datCAN_14[ABS_DIFF_40TGA>0.3 & ABS_DIFF_0714>0.3,]
	
	# all 3 comparisons
	datCAN_11_14[ABS_DIFF_40TGA>0.3 & ABS_DIFF_0711>0.3 & ABS_DIFF_0714>0.3,]