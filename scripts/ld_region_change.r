# average LD within windows and calculate change through time


######################
# calculate LD change
# run on cod node
######################
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

# width <- 1e4; stp <- 1e4; windsz='1e4'; windnm='10kb' # window parameters
width <- 3e4; stp <- 3e4; windsz='3e4'; windnm='30kb' # window parameters

# read in locus data for one population
datloc <- fread('data_2019_06_06/DP7_Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(datloc, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014

# read in LD data on cod node
dat14 <- fread(paste('zcat analysis/LOF_S_14_', windnm, '.geno.ld', sep=''))
dat11 <- fread(paste('zcat analysis/LOF_S_11_', windnm, '.geno.ld', sep=''))
dat07 <- fread(paste('zcat analysis/LOF_07_', windnm, '.geno.ld', sep=''))
dat40 <- fread(paste('zcat analysis/CAN40_', windnm, '.geno.ld', sep=''))
datMod <- fread(paste('zcat analysis/CANMod_', windnm, '.geno.ld', sep=''))

# remove ^ from a column name
setnames(dat14, 5, 'r2')
setnames(dat11, 5, 'r2')
setnames(dat07, 5, 'r2')
setnames(dat40, 5, 'r2')
setnames(datMod, 5, 'r2')

# midpoints of windows for averaging
datmean <- data.table(CHROM=character(), POSmid=numeric())
chroms <- datloc[,unique(CHROM)]
for(i in 1:length(chroms)){
	temp <- data.table(CHROM=chroms[i], POSmid=seq(width/2, datloc[CHROM==chroms[i], max(POS)], by=stp))
	datmean <- rbind(datmean, temp)
}
nrow(datmean) # 64418 (1e4), 21440 (3e4)

# round POS1 and POS2 to nearest window midpoint, once for each step in width
for(j in 1:(width/stp)){
	nm <- paste('POS1mid', j, sep='') # create column name
	for(i in 1:length(chroms)){
		dat14[CHR==chroms[i], eval(nm):=floor((POS1+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		dat11[CHR==chroms[i], eval(nm):=floor((POS1+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		dat07[CHR==chroms[i], eval(nm):=floor((POS1+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		dat40[CHR==chroms[i], eval(nm):=floor((POS1+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		datMod[CHR==chroms[i], eval(nm):=floor((POS1+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
	}
	nm <- paste('POS2mid', j, sep='') # create column name
	for(i in 1:length(chroms)){
		dat14[CHR==chroms[i], eval(nm):=floor((POS2+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		dat11[CHR==chroms[i], eval(nm):=floor((POS2+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		dat07[CHR==chroms[i], eval(nm):=floor((POS2+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		dat40[CHR==chroms[i], eval(nm):=floor((POS2+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
		datMod[CHR==chroms[i], eval(nm):=floor((POS2+(j-1)*stp)/width)*width+(j-1)*stp+width/2]
	}
}
posnms <- grep('POS[[:digit:]]mid', colnames(dat14), value=TRUE) # get column names created
	posnms


# calculate average LD within regions
bins14 <- dat14[!is.na(r2) & POS1mid1==POS2mid1,.(r2ave_14=mean(r2), r2n_14=.N), by=.(CHR, POS1mid1)]
bins11 <- dat11[!is.na(r2) & POS1mid1==POS2mid1,.(r2ave_11=mean(r2), r2n_11=.N), by=.(CHR, POS1mid1)]
bins07 <- dat07[!is.na(r2) & POS1mid1==POS2mid1,.(r2ave_07=mean(r2), r2n_07=.N), by=.(CHR, POS1mid1)] 
bins40 <- dat40[!is.na(r2) & POS1mid1==POS2mid1,.(r2ave_40=mean(r2), r2n_40=.N), by=.(CHR, POS1mid1)] 
binsMod <- datMod[!is.na(r2) & POS1mid1==POS2mid1,.(r2ave_Mod=mean(r2), r2n_Mod=.N), by=.(CHR, POS1mid1)] 

setkey(bins14, CHR, POS1mid1)
setkey(bins11, CHR, POS1mid1)
setkey(bins07, CHR, POS1mid1)
setkey(bins40, CHR, POS1mid1)
setkey(binsMod, CHR, POS1mid1)

# merge
bins <- merge(bins14, bins11, by=c('CHR', 'POS1mid1'), all=TRUE)
bins <- merge(bins, bins07, by=c('CHR', 'POS1mid1'), all=TRUE)
bins <- merge(bins, bins40, by=c('CHR', 'POS1mid1'), all=TRUE)
bins <- merge(bins, binsMod, by=c('CHR', 'POS1mid1'), all=TRUE)

nrow(bins)

# add BIN_START
bins[,BIN_START:=POS1mid1 - width/2]

# calculate change
bins[,ld_diff_0711:=r2ave_11-r2ave_07]
bins[,ld_diff_0714:=r2ave_14-r2ave_07]
bins[,ld_diff_CAN:=r2ave_Mod-r2ave_40]

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- datloc[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, CHROM)
setkey(bins, CHR)
bins <- merge(bins, chrmax[,.(CHROM, start)], by.x='CHR', by.y='CHROM')
bins[,POSgen:=POS1mid1+start]
bins[,start:=NULL]


# calculate percentiles (like Oziolor et al. 2019 Science)
ecdf0711 <- ecdf(bins$ld_diff_0711)
ecdf0714 <- ecdf(bins$ld_diff_0714)
ecdfCan <- ecdf(bins$ld_diff_CAN)
bins[,ld_region_perc0711 := ecdf0711(ld_diff_0711)]
bins[,ld_region_perc0714 := ecdf0714(ld_diff_0714)]
bins[,ld_region_percCAN := ecdfCan(ld_diff_CAN)]

	# examine
	bins[ld_region_perc0711>0.99 & ld_region_perc0714>0.99 & ld_region_percCAN>0.99,]
		
# group the 1% outliers
	# find distance among outlier regions
	setkey(bins, CHR, POS1mid1)
	bins[,dist0711 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,dist0714 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,distCAN := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	for(i in 1:length(chroms)){
		# offset the locations by one outlier locus
		temp0711 <- c(NA, bins[CHR==chroms[i] & ld_region_perc0711>0.99, POS1mid1])
		temp0714 <- c(NA, bins[CHR==chroms[i] & ld_region_perc0714>0.99, POS1mid1])
		tempCAN <- c(NA, bins[CHR==chroms[i] & ld_region_percCAN>0.99, POS1mid1])

		temp0711 <- temp0711[-length(temp0711)]
		temp0714 <- temp0714[-length(temp0714)]
		tempCAN <- tempCAN[-length(tempCAN)]

		# subtract to get the distance
		bins[CHR==chroms[i] & ld_region_perc0711>0.99, dist0711:= POS1mid1 - temp0711]
		bins[CHR==chroms[i] & ld_region_perc0714>0.99, dist0714:= POS1mid1 - temp0714]
		bins[CHR==chroms[i] & ld_region_percCAN>0.99, distCAN:= POS1mid1 - tempCAN]
	} 

	# label clusters (a little slow)
	indx <- 1
	bins[,ld_region_cluster0711 := as.numeric(NA)]
	outls <- which(bins$ld_region_perc0711>0.99)
	for(i in 1:length(outls)){
		bins$ld_region_cluster0711[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0711[outls[i+1]]>=width | is.na(bins$dist0711[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,ld_region_cluster0714 := as.numeric(NA)]
	outls <- which(bins$ld_region_perc0714>0.99)
	for(i in 1:length(outls)){
		bins$ld_region_cluster0714[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0714[outls[i+1]]>=width | is.na(bins$dist0714[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,ld_region_clusterCAN := as.numeric(NA)]
	outls <- which(bins$ld_region_percCAN>0.99)
	for(i in 1:length(outls)){
		bins$ld_region_clusterCAN[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$distCAN[outls[i+1]]>=width | is.na(bins$distCAN[outls[i+1]])) indx <- indx+1
		}
	}

	# examine
	bins[,range(table(ld_region_cluster0711))] # all clusters only 1 window wide
	bins[,range(table(ld_region_cluster0714))] # all clusters only 1 window wide
	bins[,range(table(ld_region_clusterCAN))] # all clusters only 1 window wide


# save region means data
filenm <- paste('analysis/ld_change_region_', windsz, '.rds', sep='')
filenm
saveRDS(bins, file=filenm)


###########################
## rank the outlier regions
##########################
outlreg0711 = bins[,.(CHROM=unique(CHR), POSmin=min(POS1mid1), POSmax=max(POS1mid1), score=sum(ld_diff_0711)), by=ld_region_cluster0711]
outlreg0714 = bins[,.(CHROM=unique(CHR), POSmin=min(POS1mid1), POSmax=max(POS1mid1), score=sum(ld_diff_0714)), by=ld_region_cluster0714]
outlregCAN = bins[,.(CHROM=unique(CHR), POSmin=min(POS1mid1), POSmax=max(POS1mid1), score=sum(ld_diff_CAN)), by=ld_region_clusterCAN]

setkey(outlreg0711, score)
setkey(outlreg0714, score)
setkey(outlregCAN, score)

tail(outlreg0711, 5)
tail(outlreg0714, 5)
tail(outlregCAN, 5)



###########################################
# plot LD change from regions
###########################################
require(data.table)

# 1e4
bins <- readRDS('analysis/ld_change_region_1e4.rds'); width='1e4'

# create chromosome labels from an approximate midpoint
chrs <- bins[, .(pos=mean(POSgen)), by=CHR]

# add a vector for color by LG
bins[,lgcol := 'grey40']
bins[CHR %in% chrs$CHR[seq(2, nrow(chrs),by=2)], lgcol := 'grey60']


cols = c('grey', 'purple')
quartz(height=6, width=8)
# png(height=6, width=8, units='in', res=300, file=paste('figures/ld_change_vs_pos_NEA_CAN_runmean', width, '.png', sep=''))
par(mfrow=c(3,1), mai=c(0.3, 0.7, 0.1, 0.5), mgp=c(3.2, 0.6, 0), las=1, tcl=-0.3)
bins[ld_region_perc0711<=0.99 | ld_region_perc0714<=0.99 | ld_region_percCAN<=0.99, plot(POSgen, ld_diff_0711, type='p', cex=0.2, lwd=0.3, xaxt='n', xlab='', ylab='LD change Lof0711', ylim=c(-1,1), col=lgcol, las=1)]
bins[ld_region_perc0711>0.99 & ld_region_perc0714>0.99 & ld_region_percCAN>0.99, points(POSgen, ld_diff_0711, cex=0.5, col=cols[2])]

axis(side=1, at=chrs$pos, labels=gsub('LG','', chrs$CHR), tcl=-0.3, cex.axis=0.7)

legend('topright', legend=c(paste(width, 'bp moving window average'), 'Outlier shared among all 3 pops'), pch=1, col=cols, bty='n')

bins[ld_region_perc0711<=0.99 | ld_region_perc0714<=0.99 | ld_region_percCAN<=0.99, plot(POSgen, ld_diff_0714, type='p', cex=0.2, lwd=0.3, xaxt='n', xlab='', ylab='LD change Lof0714', ylim=c(-1,1), col=lgcol, las=1)]
bins[ld_region_perc0711>0.99 & ld_region_perc0714>0.99 & ld_region_percCAN>0.99,points(POSgen, ld_diff_0714, cex=0.5, col=cols[2])]

axis(side=1, at=chrs$pos, labels=gsub('LG','', chrs$CHR), tcl=-0.3, cex.axis=0.7)

bins[ld_region_perc0711<=0.99 | ld_region_perc0714<=0.99 | ld_region_percCAN<=0.99, plot(POSgen, ld_diff_CAN, type='p', cex=0.2, lwd=0.3, xaxt='n', xlab='', ylab='LD change CAN', ylim=c(-1,1), col=lgcol, las=1)]
bins[ld_region_perc0711>0.99 & ld_region_perc0714>0.99 & ld_region_percCAN>0.99, points(POSgen, ld_diff_CAN, cex=0.5, col=cols[2])]

axis(side=1, at=chrs$pos, labels=gsub('LG','', chrs$CHR), tcl=-0.3, cex.axis=0.7)

dev.off()
