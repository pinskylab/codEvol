## Plot allele frequency change vs. genome position
## Compare across different data sets

###########################
# load functions
###########################
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

addchroms <- function(dat){
	col='black'
	lgs <- sort(unique(dat[,CHROM]))
	for(j in 1:length(lgs)){
		rng <- range(dat[CHROM==lgs[j], POSgen])
		if(j %% 2 == 0) lines(x=rng/1e6, y=c(0,0), col=col, lwd=2, lend=2)
		if(j %% 2 == 1) lines(x=rng/1e6, y=c(0.02,0.02), col=col, lwd=2, lend=2)
		text(x=mean(rng/1e6), y=0.03, labels=lgs[j], col=col, cex=0.5)
	
	}
}

#####################
# read and prep in data
#####################

dat14 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat14, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
dat11 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat11, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
datCAN <- fread('data_2019_03_18/Frequency_table_CAN_40_TGA.txt', header=TRUE); setnames(datCAN, 3:7, c('N_CHR_Can40', 'Freq_Can40', 'N_CHR_CanMod', 'Freq_CanMod', 'ABS_DIFF_Can')) # for Canada
	nrow(dat11)
	nrow(dat14)
	nrow(datCAN)


# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- dat14[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, CHROM)

# merge nucleotide position into the frequency files
setkey(dat14, CHROM)
dat14 <- dat14[chrmax[,.(CHROM, start)], ]
dat14[,POSgen:=POS+start]
dat14[,start:=NULL]

setkey(dat11, CHROM)
dat11 <- dat11[chrmax[,.(CHROM, start)], ]
dat11[,POSgen:=POS+start]
dat11[,start:=NULL]

setkey(datCAN, CHROM)
datCAN <- datCAN[chrmax[,.(CHROM, start)], ]
datCAN[,POSgen:=POS+start]
datCAN[,start:=NULL]

# combine the datasets to look at correlations across them
setkey(dat11, CHROM, POS)
setkey(dat14, CHROM, POS)
setkey(datCAN, CHROM, POS)

dat <- merge(datCAN, dat14[,.(CHROM, POS, N_CHR_07, N_CHR_14, Freq_07, Freq_14, ABS_DIFF_0714)], all=TRUE)
dat <- merge(dat, dat11[,.(CHROM, POS, N_CHR_11, Freq_11, ABS_DIFF_0711)], all=TRUE)
	nrow(dat)
	dat


#############################
# running mean of allele frequency difference 
# (if saved, can skip this)
#############################
width <- 1e4; stp <- 1e4; windsz='1e4' 

# midpoints of windows
datmean <- data.table(CHROM=character(), POSmid=numeric())
chroms <- dat[,unique(CHROM)]
for(i in 1:length(chroms)){
	temp <- data.table(CHROM=chroms[i], POSmid=seq(0, dat[CHROM==chroms[i], max(POS)], by=stp))
	datmean <- rbind(datmean, temp)
}

nrow(datmean)
setkey(dat, CHROM, POS)
setkey(datmean, CHROM, POSmid)

# round POS to nearest POSmid, once for each step in width
for(j in 1:(width/stp)){
	nm <- paste('POSmid', j, sep='') # create column name
	for(i in 1:length(chroms)){
		dat[CHROM==chroms[i], eval(nm):=floor((POS+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
	}
}
posnms <- grep('POSmid', colnames(dat), value=TRUE) # get column names created
	posnms

# calculate moving window means
datmean <- data.table(CHROM=character(), POSmid=numeric(), ABS_DIFF_0711=numeric(), ABS_DIFF_0714=numeric(), ABS_DIFF_Can=numeric())
for(j in 1:length(posnms)){
	for(i in 1:length(chroms)){
		temp <- dat[CHROM==chroms[i], .(CHROM=unique(CHROM), ABS_DIFF_0711=mean(ABS_DIFF_0711, na.rm=TRUE), ABS_DIFF_0714=mean(ABS_DIFF_0714, na.rm=TRUE), ABS_DIFF_Can=mean(ABS_DIFF_Can, na.rm=TRUE)), by=eval(posnms[j])]
		setnames(temp, eval(posnms[j]), 'POSmid')
		datmean <- rbind(datmean, temp)
	}
}
dim(dat)
dim(datmean)

# make POSgen column
setkey(datmean, CHROM)
setkey(chrmax, CHROM)
datmean <- datmean[chrmax[,.(CHROM, start)], ]
datmean[,POSgen:=POSmid+start]
datmean[,start:=NULL]


# sort
setkey(datmean, CHROM, POSmid)

# calculate percentiles (like Oziolor et al. 2019 Science)
ecdf0711 <- ecdf(datmean$ABS_DIFF_0711)
ecdf0714 <- ecdf(datmean$ABS_DIFF_0714)
ecdfCan <- ecdf(datmean$ABS_DIFF_Can)
datmean[,perc0711 := ecdf0711(ABS_DIFF_0711)]
datmean[,perc0714 := ecdf0714(ABS_DIFF_0714)]
datmean[,percCAN := ecdfCan(ABS_DIFF_Can)]

	# examine
	datmean[perc0711>0.999 & perc0714>0.999 & percCAN>0.999,]
	
	datmean[plot(perc0711, perc0714, cex=0.5, col='#00000005')]
	datmean[plot(perc0711, percCAN, cex=0.5, col='#00000005')]
	
# group the 1% outliers
	# find distance among outlier regions
	setkey(datmean, CHROM, POSmid)
	datmean[,dist0711 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	datmean[,dist0714 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	datmean[,distCAN := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	for(i in 1:length(chroms)){
		temp0711 <- c(NA, datmean[CHROM==chroms[i] & perc0711>0.99, POSmid])
		temp0714 <- c(NA, datmean[CHROM==chroms[i] & perc0714>0.99, POSmid])
		tempCAN <- c(NA, datmean[CHROM==chroms[i] & percCAN>0.99, POSmid])

		temp0711 <- temp0711[-length(temp0711)]
		temp0714 <- temp0714[-length(temp0714)]
		tempCAN <- tempCAN[-length(tempCAN)]

		datmean[CHROM==chroms[i] & perc0711>0.99, dist0711:= POSmid - temp0711]
		datmean[CHROM==chroms[i] & perc0714>0.99, dist0714:= POSmid - temp0714]
		datmean[CHROM==chroms[i] & percCAN>0.99, distCAN:= POSmid - tempCAN]
	} 

	# label clusters (a little slow)
	indx <- 1
	datmean[,cluster0711 := as.numeric(NA)]
	outls <- which(datmean$perc0711>0.99)
	for(i in 1:length(outls)){
		datmean$cluster0711[outls[i]] <- indx
		if(i < length(outls)){
			if(datmean$dist0711[outls[i+1]]>=width | is.na(datmean$dist0711[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	datmean[,cluster0714 := as.numeric(NA)]
	outls <- which(datmean$perc0714>0.99)
	for(i in 1:length(outls)){
		datmean$cluster0714[outls[i]] <- indx
		if(i < length(outls)){
			if(datmean$dist0714[outls[i+1]]>=width | is.na(datmean$dist0714[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	datmean[,clusterCAN := as.numeric(NA)]
	outls <- which(datmean$percCAN>0.99)
	for(i in 1:length(outls)){
		datmean$clusterCAN[outls[i]] <- indx
		if(i < length(outls)){
			if(datmean$distCAN[outls[i+1]]>=width | is.na(datmean$distCAN[outls[i+1]])) indx <- indx+1
		}
	}

	

# save region means data
filenm <- paste('analysis/Frequency_table_ABS_DIFF_runmean', windsz, '.rds', sep='')
filenm
saveRDS(datmean, file=filenm)


# add region outlier stats to locus-by-locus frequency change file
dat2 <- merge(dat, datmean[, .(CHROM, POSmid, perc0711, perc0714, percCAN, cluster0711, cluster0714, clusterCAN)], by.x=c('CHROM', 'POSmid1'), by.y=c('CHROM', 'POSmid'))
	dim(dat)
	dim(dat2)
	dim(datmean)
setnames(dat2, c('perc0711', 'perc0714', 'percCAN', 'cluster0711', 'cluster0714', 'clusterCAN'), c('region10kb_perc0711', 'region10kb_perc0714', 'region10kb_percCAN', 'region10kb_cluster0711', 'region10kb_cluster0714', 'region10kb_clusterCAN'))

# save locus data with region data appended
filenm2 <- paste('analysis/Frequency_table_ABS_DIFF_w_region_outliers_', windsz, '.rds', sep='')
filenm2
saveRDS(dat2, file=filenm2)


###########################
## rank the outlier regions (see Oziolor et al. 2019 Science)
##########################
outlreg0711 = datmean[,.(CHROM=unique(CHROM), POSmin=min(POSmid), POSmax=max(POSmid), score=sum(ABS_DIFF_0711)), by=cluster0711]
outlreg0714 = datmean[,.(CHROM=unique(CHROM), POSmin=min(POSmid), POSmax=max(POSmid), score=sum(ABS_DIFF_0714)), by=cluster0714]
outlregCAN = datmean[,.(CHROM=unique(CHROM), POSmin=min(POSmid), POSmax=max(POSmid), score=sum(ABS_DIFF_Can)), by=clusterCAN]

setkey(outlreg0711, score)
setkey(outlreg0714, score)
setkey(outlregCAN, score)

tail(outlreg0711[CHROM != 'LG01'], 20)
tail(outlreg0714[CHROM != 'LG01'], 20)
tail(outlregCAN[CHROM != 'LG01'], 20)

##########################
# plot frequency change
# no smoothing
##########################

cols = c('black', 'blue', 'red')
#dat14[,plot(POSgen, ABS_DIFF, pch=16, cex=0.3)]
quartz(height=12, width=8)
# png(height=12, width=8, units='in', res=300, file='figures/abs_diff_vs_pos_NEA_CAN_raw.png')
par(mfrow=c(4,1), mai=c(0.5, 1, 0.2, 0.5))
dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[1])]
dat11[,points(POSgen/1e6, ABS_DIFF_0711, cex=0.2, lwd=0.3, col=cols[2])]
datCAN[,points(POSgen/1e6, ABS_DIFF_Can, cex=0.2, lwd=0.3, col=cols[3])]

addchroms(dat14)
legend('topright', legend=c('NEA 1907-2014', 'NEA 1907-2011', 'Canada'), lwd=1, col=cols, bty='n', cex=0.5)

dat14[,plot(POSgen/1e6, ABS_DIFF_0714, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[1])]
addchroms(dat14)

dat11[,plot(POSgen/1e6, ABS_DIFF_0711, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[2])]
addchroms(dat14)

datCAN[,plot(POSgen/1e6, ABS_DIFF_Can, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,1), col=cols[3])]
addchroms(dat14)


dev.off()


##########################
# plot frequency change
# no smoothing
# color if in a 10kb outlier region shared across populations
##########################

# 1e4
dat2 <- readRDS('analysis/Frequency_table_ABS_DIFF_w_region_outliers_1e4.rds'); width='1e4'

# 1e6
# NOTE: loads dat14mean and dat11mean, without _25k or _150k suffixes
# load('analysis/Frequency_table_Lof07_Lof14_runmean1e6.rdata'); load('analysis/Frequency_table_Lof07_Lof11_runmean1e6.rdata'); ylims=c(0,0.15); width='1e6'
# load('analysis/Frequency_table_CAN_runmean1e6.rdata')


cols = c('grey', '#e0ecf4', '#9ebcda', '#8856a7') # not outlier, single pop outlier, 2-pop outlier, all pop outlier
quartz(height=6, width=8)
# png(height=6, width=8, units='in', res=300, file=paste('figures/abs_diff_vs_pos_NEA_CAN_bylocus_regionoutliers_', width, '.png', sep=''))
par(mfrow=c(3,1), mai=c(0.5, 1, 0.2, 0.5))
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN<=0.99, plot(POSgen/1e6, ABS_DIFF_0711, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Freq change Lof0711', ylim=c(0,1), col=cols[1])]
dat2[region10kb_perc0711>0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_0711, cex=0.2, col=cols[2])] # 1 pop
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714>0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_0711, cex=0.2, col=cols[2])]
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0711, cex=0.2, col=cols[2])]
dat2[region10kb_perc0711>0.99 & region10kb_perc0714>0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_0711, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711>0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0711, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714>0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0711, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711>0.99 & region10kb_perc0714>0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0711, cex=0.5, col=cols[4])]
addchroms(dat2)

legend('topright', legend=c(paste('Not in', width, 'bp 1% outlier region'), 'In single pop outlier region', 'In 2-pop outlier region', 'In 3-pop outlier region'), pch=1, col=cols, bty='n', cex=0.6)

dat2[region10kb_perc0711<=0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN<=0.99,plot(POSgen/1e6, ABS_DIFF_0714, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Freq change Lof0714', ylim=c(0,1), col=cols[1])]
dat2[region10kb_perc0711>0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_0714, cex=0.2, col=cols[2])] # 1 pop
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714>0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_0714, cex=0.2, col=cols[2])]
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0714, cex=0.2, col=cols[2])]
dat2[region10kb_perc0711>0.99 & region10kb_perc0714>0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_0714, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711>0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0714, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714>0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0714, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711>0.99 & region10kb_perc0714>0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_0714, cex=0.5, col=cols[4])]
addchroms(dat2)

dat2[region10kb_perc0711<=0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN<=0.99,plot(POSgen/1e6, ABS_DIFF_Can, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Freq change CAN', ylim=c(0,1), col=cols[1])]
dat2[region10kb_perc0711>0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_Can, cex=0.2, col=cols[2])] # 1 pop
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714>0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_Can, cex=0.2, col=cols[2])]
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_Can, cex=0.2, col=cols[2])]
dat2[region10kb_perc0711>0.99 & region10kb_perc0714>0.99 & region10kb_percCAN<=0.99, points(POSgen/1e6, ABS_DIFF_Can, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711>0.99 & region10kb_perc0714<=0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_Can, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711<=0.99 & region10kb_perc0714>0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_Can, cex=0.5, col=cols[3])]  # 2 pop
dat2[region10kb_perc0711>0.99 & region10kb_perc0714>0.99 & region10kb_percCAN>0.99, points(POSgen/1e6, ABS_DIFF_Can, cex=0.5, col=cols[4])]
addchroms(dat2)

dev.off()


###########################################
# plot frequency change from running mean
###########################################
# 1e4
datmean <- readRDS('analysis/Frequency_table_ABS_DIFF_runmean1e4.rds'); width='1e4'

# 1e6
# NOTE: loads dat14mean and dat11mean, without _25k or _150k suffixes
# load('analysis/Frequency_table_Lof07_Lof14_runmean1e6.rdata'); load('analysis/Frequency_table_Lof07_Lof11_runmean1e6.rdata'); ylims=c(0,0.15); width='1e6'
# load('analysis/Frequency_table_CAN_runmean1e6.rdata')


cols = c('grey', 'purple')
quartz(height=6, width=8)
# png(height=6, width=8, units='in', res=300, file=paste('figures/abs_diff_vs_pos_NEA_CAN_', width, '.png', sep=''))
par(mfrow=c(3,1), mai=c(0.5, 1, 0.2, 0.5))
datmean[,plot(POSgen/1e6, ABS_DIFF_0711, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Freq change Lof0711', ylim=c(0,1), col=cols[1])]
datmean[perc0711>0.99 & perc0714>0.99 & percCAN>0.99,points(POSgen/1e6, ABS_DIFF_0711, cex=0.5, col=cols[2])]
addchroms(datmean)

legend('topright', legend=c(paste(width, 'bp moving window average'), 'Outlier shared among all 3 pops'), pch=1, col=cols, bty='n')

datmean[,plot(POSgen/1e6, ABS_DIFF_0714, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Freq change Lof0714', ylim=c(0,1), col=cols[1])]
datmean[perc0711>0.99 & perc0714>0.99 & percCAN>0.99,points(POSgen/1e6, ABS_DIFF_0714, cex=0.5, col=cols[2])]
addchroms(datmean)

datmean[,plot(POSgen/1e6, ABS_DIFF_Can, type='p', cex=0.2, lwd=0.3, xlab='Position (Mb)', ylab='Freq change CAN', ylim=c(0,1), col=cols[1])]
datmean[perc0711>0.99 & perc0714>0.99 & percCAN>0.99,points(POSgen/1e6, ABS_DIFF_Can, cex=0.5, col=cols[2])]
addchroms(datmean)

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







###############################################
# plot histogram of allele frequency changes
###############################################
datmean <- readRDS('analysis/Frequency_table_ABS_DIFF_runmean1e4.rds'); width='1e4'

bks <- seq(0,1,by=0.1)
cols = c('black', 'blue', 'red')


hist0711 <- hist(datmean$ABS_DIFF_0711, breaks=bks, plot=FALSE)
hist0714 <- hist(datmean$ABS_DIFF_0714, breaks=bks, plot=FALSE)
histCAN <- hist(datmean$ABS_DIFF_Can, breaks=bks, plot=FALSE)

quartz(height=4, width=4)
par(las=1)
plot(hist0711$mids, log10(hist0711$counts), type='l', col=cols[1], xlab='Allele frequency change', ylab='log10(counts)')
lines(hist0714$mids, log10(hist0714$counts), col=cols[2])
lines(histCAN$mids, log10(histCAN$counts), col=cols[3])

