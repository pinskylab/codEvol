# calculate Tajima's D change through time

# require(data.table)
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

width <- 1e4; stp <- 1e4; windsz='1e4'# window parameters



######################
# calculate Tajima's D change
######################

# read in locus data for one population
datloc <- fread('data_2019_03_18/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(datloc, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014

# read in LD data on cod node
dat14 <- fread('analysis/LOF_S_14_10kb.Tajima.D')
dat11 <- fread('analysis/LOF_S_11_10kb.Tajima.D')
dat07 <- fread('analysis/LOF_07_10kb.Tajima.D')
dat40 <- fread('analysis/CAN40_10kb.Tajima.D')
datMod <- fread('analysis/CANMod_10kb.Tajima.D')

# set col names
setnames(dat14, 3:4, c('N_SNPS_14', 'TajimaD_14'))
setnames(dat11, 3:4, c('N_SNPS_11', 'TajimaD_11'))
setnames(dat07, 3:4, c('N_SNPS_07', 'TajimaD_07'))
setnames(dat40, 3:4, c('N_SNPS_40', 'TajimaD_40'))
setnames(datMod, 3:4, c('N_SNPS_Mod', 'TajimaD_Mod'))

# merge
bins <- merge(dat14, dat11, by=c('CHROM', 'BIN_START'))
bins <- merge(bins, dat07, by=c('CHROM', 'BIN_START'))
bins <- merge(bins, dat40, by=c('CHROM', 'BIN_START'))
bins <- merge(bins, datMod, by=c('CHROM', 'BIN_START'))

# calculate change
bins[,D_diff_0711:=TajimaD_11-TajimaD_07]
bins[,D_diff_0714:=TajimaD_14-TajimaD_07]
bins[,D_diff_CAN:=TajimaD_Mod-TajimaD_40]

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- datloc[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, CHROM)
setkey(bins, CHROM)
bins <- merge(bins, chrmax[,.(CHROM, start)], by='CHROM')
bins[,POSgen:=BIN_START+start]
bins[,start:=NULL]

# calculate percentiles (like Oziolor et al. 2019 Science)
ecdf0711 <- ecdf(bins$D_diff_0711)
ecdf0714 <- ecdf(bins$D_diff_0714)
ecdfCan <- ecdf(bins$D_diff_CAN)
bins[,D_region10kb_perc0711 := ecdf0711(D_diff_0711)]
bins[,D_region10kb_perc0714 := ecdf0714(D_diff_0714)]
bins[,D_region10kb_percCAN := ecdfCan(D_diff_CAN)]

	# examine
	bins[D_region10kb_perc0711>0.99 & D_region10kb_perc0714>0.99 & D_region10kb_percCAN>0.99,]
	bins[D_region10kb_perc0711<0.01 & D_region10kb_perc0714<0.01 & D_region10kb_percCAN<0.01,]
		
# group the 1% outliers (increasing or decreasing)
	# find distance among increasing outlier regions
	chroms <- datloc[,unique(CHROM)]
	setkey(bins, CHROM, BIN_START)
	bins[,dist0711 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,dist0714 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,distCAN := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	for(i in 1:length(chroms)){
		# offset the locations by one outlier locus
		temp0711 <- c(NA, bins[CHROM==chroms[i] & D_region10kb_perc0711>0.995, BIN_START])
		temp0714 <- c(NA, bins[CHROM==chroms[i] & D_region10kb_perc0714>0.995, BIN_START])
		tempCAN <- c(NA, bins[CHROM==chroms[i] & D_region10kb_percCAN>0.995, BIN_START])

		temp0711 <- temp0711[-length(temp0711)]
		temp0714 <- temp0714[-length(temp0714)]
		tempCAN <- tempCAN[-length(tempCAN)]

		# subtract to get the distance
		bins[CHROM==chroms[i] & D_region10kb_perc0711>0.995, dist0711:= as.numeric(BIN_START - temp0711)]
		bins[CHROM==chroms[i] & D_region10kb_perc0714>0.995, dist0714:= as.numeric(BIN_START - temp0714)]
		bins[CHROM==chroms[i] & D_region10kb_percCAN>0.995, distCAN:= as.numeric(BIN_START - tempCAN)]
	} 

	# label increasing D clusters (a little slow)
	indx <- 1
	bins[,D_region10kb_cluster0711 := as.numeric(NA)]
	outls <- which(bins$D_region10kb_perc0711>0.995)
	for(i in 1:length(outls)){
		bins$D_region10kb_cluster0711[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0711[outls[i+1]]>=width | is.na(bins$dist0711[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,D_region10kb_cluster0714 := as.numeric(NA)]
	outls <- which(bins$D_region10kb_perc0714>0.995)
	for(i in 1:length(outls)){
		bins$D_region10kb_cluster0714[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0714[outls[i+1]]>=width | is.na(bins$dist0714[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,D_region10kb_clusterCAN := as.numeric(NA)]
	outls <- which(bins$D_region10kb_percCAN>0.995)
	for(i in 1:length(outls)){
		bins$D_region10kb_clusterCAN[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$distCAN[outls[i+1]]>=width | is.na(bins$distCAN[outls[i+1]])) indx <- indx+1
		}
	}

	# find distance among decreasing outlier regions
	chroms <- datloc[,unique(CHROM)]
	setkey(bins, CHROM, BIN_START)
	bins[,dist0711dec := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,dist0714dec := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,distCANdec := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	for(i in 1:length(chroms)){
		# offset the locations by one outlier locus
		temp0711 <- c(NA, bins[CHROM==chroms[i] & D_region10kb_perc0711<0.005, BIN_START])
		temp0714 <- c(NA, bins[CHROM==chroms[i] & D_region10kb_perc0714<0.005, BIN_START])
		tempCAN <- c(NA, bins[CHROM==chroms[i] & D_region10kb_percCAN<0.005, BIN_START])

		temp0711 <- temp0711[-length(temp0711)]
		temp0714 <- temp0714[-length(temp0714)]
		tempCAN <- tempCAN[-length(tempCAN)]

		# subtract to get the distance
		bins[CHROM==chroms[i] & D_region10kb_perc0711<0.005, dist0711dec:= as.numeric(BIN_START - temp0711)]
		bins[CHROM==chroms[i] & D_region10kb_perc0714<0.005, dist0714dec:= as.numeric(BIN_START - temp0714)]
		bins[CHROM==chroms[i] & D_region10kb_percCAN<0.005, distCANdec:= as.numeric(BIN_START - tempCAN)]
	} 

	# label increasing D clusters (a little slow)
	indx <- 1
	bins[,D_region10kb_cluster0711dec := as.numeric(NA)]
	outls <- which(bins$D_region10kb_perc0711<0.005)
	for(i in 1:length(outls)){
		bins$D_region10kb_cluster0711dec[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0711dec[outls[i+1]]>=width | is.na(bins$dist0711dec[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,D_region10kb_cluster0714dec := as.numeric(NA)]
	outls <- which(bins$D_region10kb_perc0714<0.005)
	for(i in 1:length(outls)){
		bins$D_region10kb_cluster0714dec[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0714dec[outls[i+1]]>=width | is.na(bins$dist0714dec[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,D_region10kb_clusterCANdec := as.numeric(NA)]
	outls <- which(bins$D_region10kb_percCAN<0.005)
	for(i in 1:length(outls)){
		bins$D_region10kb_clusterCANdec[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$distCANdec[outls[i+1]]>=width | is.na(bins$distCANdec[outls[i+1]])) indx <- indx+1
		}
	}


# save region means data
filenm <- paste('analysis/tajimasD_change_region_', windsz, '.rds', sep='')
filenm
saveRDS(bins, file=filenm)


###########################
## rank the outlier regions
## NOT YET UPDATED FOR D
##########################
outlreg0711 = bins[,.(CHROM=unique(CHR), POSmin=min(POS1mid1), POSmax=max(POS1mid1), score=sum(ld_diff_0711)), by=ld_region10kb_cluster0711]
outlreg0714 = bins[,.(CHROM=unique(CHR), POSmin=min(POS1mid1), POSmax=max(POS1mid1), score=sum(ld_diff_0714)), by=ld_region10kb_cluster0714]
outlregCAN = bins[,.(CHROM=unique(CHR), POSmin=min(POS1mid1), POSmax=max(POS1mid1), score=sum(ld_diff_CAN)), by=ld_region10kb_clusterCAN]

setkey(outlreg0711, score)
setkey(outlreg0714, score)
setkey(outlregCAN, score)

tail(outlreg0711, 5)
tail(outlreg0714, 5)
tail(outlregCAN, 5)



###########################################
# plot D change from regions
###########################################
require(data.table)

# read in data
bins <- readRDS('analysis/tajimasD_change_region_1e4.rds'); width='1e4'

# create chromosome labels from an approximate midpoint
chrs <- bins[, .(pos=mean(POSgen)), by=CHROM]

# add a vector for color by LG
bins[,lgcol := 'grey40']
bins[CHROM %in% chrs$CHROM[seq(2, nrow(chrs),by=2)], lgcol := 'grey60']


cols = c('grey', 'purple', 'red')
quartz(height=6, width=8)
# png(height=6, width=8, units='in', res=300, file=paste('figures/D_change_vs_pos_NEA_CAN_runmean', width, '.png', sep=''))
par(mfrow=c(3,1), mai=c(0.3, 0.7, 0.1, 0.5), mgp=c(3.2, 0.6, 0), las=1, tcl=-0.3)
bins[(D_region10kb_perc0711<=0.995 | D_region10kb_perc0714<=0.995 | D_region10kb_percCAN<=0.995) & (D_region10kb_perc0711>=0.005 | D_region10kb_perc0714>=0.005 | D_region10kb_percCAN>=0.005), plot(POSgen, D_diff_0711, type='p', cex=0.2, lwd=0.3, xaxt='n', xlab='', ylab='D change Lof0711', col=lgcol, las=1)]
bins[D_region10kb_perc0711>0.995 & D_region10kb_perc0714>0.995 & D_region10kb_percCAN>0.995, points(POSgen, D_diff_0711, cex=0.5, col=cols[2])]
bins[D_region10kb_perc0711<0.005 & D_region10kb_perc0714<0.005 & D_region10kb_percCAN<0.005, points(POSgen, D_diff_0711, cex=0.5, col=cols[3])]

axis(side=1, at=chrs$pos, labels=gsub('LG','', chrs$CHROM), tcl=-0.3, cex.axis=0.7)

legend('topright', legend=c(paste(width, 'bp moving window average'), '+Outlier shared among all 3 pops', '-Outlier shared among all 3 pops'), pch=1, col=cols, bty='n')

bins[(D_region10kb_perc0711<=0.995 | D_region10kb_perc0714<=0.995 | D_region10kb_percCAN<=0.995) & (D_region10kb_perc0711>=0.005 | D_region10kb_perc0714>=0.005 | D_region10kb_percCAN>=0.005), plot(POSgen, D_diff_0714, type='p', cex=0.2, lwd=0.3, xaxt='n', xlab='', ylab='D change Lof0714', col=lgcol, las=1)]
bins[D_region10kb_perc0711>0.995 & D_region10kb_perc0714>0.995 & D_region10kb_percCAN>0.995,points(POSgen, D_diff_0714, cex=0.5, col=cols[2])]
bins[D_region10kb_perc0711<0.005 & D_region10kb_perc0714<0.005 & D_region10kb_percCAN<0.005, points(POSgen, D_diff_0714, cex=0.5, col=cols[3])]

axis(side=1, at=chrs$pos, labels=gsub('LG','', chrs$CHROM), tcl=-0.3, cex.axis=0.7)

bins[D_region10kb_perc0711<=0.99 | D_region10kb_perc0714<=0.99 | D_region10kb_percCAN<=0.99, plot(POSgen, D_diff_CAN, type='p', cex=0.2, lwd=0.3, xaxt='n', xlab='', ylab='D change CAN', col=lgcol, las=1)]
bins[D_region10kb_perc0711>0.995 & D_region10kb_perc0714>0.995 & D_region10kb_percCAN>0.995, points(POSgen, D_diff_CAN, cex=0.5, col=cols[2])]
bins[D_region10kb_perc0711<0.005 & D_region10kb_perc0714<0.005 & D_region10kb_percCAN<0.005, points(POSgen, D_diff_CAN, cex=0.5, col=cols[3])]

axis(side=1, at=chrs$pos, labels=gsub('LG','', chrs$CHROM), tcl=-0.3, cex.axis=0.7)

dev.off()
