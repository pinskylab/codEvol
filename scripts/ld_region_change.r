# average LD within windows and calculate change through time

#require(data.table)
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

width <- 1e4; stp <- 1e4; windsz='1e4'# window parameters


######################
# calculate LD decay
######################

# read in locus data for one population
datloc <- fread('data_2019_03_18/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(datloc, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014

# read in LD data on cod node
dat14 <- fread('zcat analysis/LOF_S_14_10kb.geno.ld')
dat11 <- fread('zcat analysis/LOF_S_11_10kb.geno.ld')
dat07 <- fread('zcat analysis/LOF_07_10kb.geno.ld')
dat40 <- fread('zcat analysis/CAN40_10kb.geno.ld')
datMod <- fread('zcat analysis/CANMod_10kb.geno.ld')

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
	temp <- data.table(CHROM=chroms[i], POSmid=seq(0, datloc[CHROM==chroms[i], max(POS)], by=stp))
	datmean <- rbind(datmean, temp)
}
nrow(datmean) # 64429

# round POS1 and POS2 to nearest window midpoint, once for each step in width
for(j in 1:(width/stp)){
	nm <- paste('POS1mid', j, sep='') # create column name
	for(i in 1:length(chroms)){
		dat14[CHR==chroms[i], eval(nm):=floor((POS1+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		dat11[CHR==chroms[i], eval(nm):=floor((POS1+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		dat07[CHR==chroms[i], eval(nm):=floor((POS1+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		dat40[CHR==chroms[i], eval(nm):=floor((POS1+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		datMod[CHR==chroms[i], eval(nm):=floor((POS1+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
	}
	nm <- paste('POS2mid', j, sep='') # create column name
	for(i in 1:length(chroms)){
		dat14[CHR==chroms[i], eval(nm):=floor((POS2+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		dat11[CHR==chroms[i], eval(nm):=floor((POS2+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		dat07[CHR==chroms[i], eval(nm):=floor((POS2+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		dat40[CHR==chroms[i], eval(nm):=floor((POS2+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
		datMod[CHR==chroms[i], eval(nm):=floor((POS2+width/2-(j-1)*stp)/width)*width+(j-1)*stp]
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
bins <- merge(bins14, bins11, by=c('CHR', 'POS1mid1'))
bins <- merge(bins, bins07, by=c('CHR', 'POS1mid1'))
bins <- merge(bins, bins40, by=c('CHR', 'POS1mid1'))
bins <- merge(bins, binsMod, by=c('CHR', 'POS1mid1'))

# calculate change
bins[,ld_diff_0711:=r2ave_11-r2ave_07]
bins[,ld_diff_0714:=r2ave_14-r2ave_07]
bins[,ld_diff_CAN:=r2ave_Mod-r2ave_40]

# calculate percentiles (like Oziolor et al. 2019 Science)
ecdf0711 <- ecdf(bins$ld_diff_0711)
ecdf0714 <- ecdf(bins$ld_diff_0714)
ecdfCan <- ecdf(bins$ld_diff_CAN)
bins[,ld_region10kb_perc0711 := ecdf0711(ld_diff_0711)]
bins[,ld_region10kb_perc0714 := ecdf0714(ld_diff_0714)]
bins[,ld_region10kb_percCAN := ecdfCan(ld_diff_CAN)]

	# examine
	bins[ld_region10kb_perc0711>0.99 & ld_region10kb_perc0714>0.99 & ld_region10kb_percCAN>0.99,]
		
# group the 1% outliers
	# find distance among outlier regions
	setkey(bins, CHR, POS1mid1)
	bins[,dist0711 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,dist0714 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	bins[,distCAN := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	for(i in 1:length(chroms)){
		# offset the locations by one outlier locus
		temp0711 <- c(NA, bins[CHR==chroms[i] & ld_region10kb_perc0711>0.99, POS1mid1])
		temp0714 <- c(NA, bins[CHR==chroms[i] & ld_region10kb_perc0714>0.99, POS1mid1])
		tempCAN <- c(NA, bins[CHR==chroms[i] & ld_region10kb_percCAN>0.99, POS1mid1])

		temp0711 <- temp0711[-length(temp0711)]
		temp0714 <- temp0714[-length(temp0714)]
		tempCAN <- tempCAN[-length(tempCAN)]

		# subtract to get the distance
		bins[CHR==chroms[i] & ld_region10kb_perc0711>0.99, dist0711:= POS1mid1 - temp0711]
		bins[CHR==chroms[i] & ld_region10kb_perc0714>0.99, dist0714:= POS1mid1 - temp0714]
		bins[CHR==chroms[i] & ld_region10kb_percCAN>0.99, distCAN:= POS1mid1 - tempCAN]
	} 

	# label clusters (a little slow)
	indx <- 1
	bins[,ld_region10kb_cluster0711 := as.numeric(NA)]
	outls <- which(bins$ld_region10kb_perc0711>0.99)
	for(i in 1:length(outls)){
		bins$ld_region10kb_cluster0711[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0711[outls[i+1]]>=width | is.na(bins$dist0711[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,ld_region10kb_cluster0714 := as.numeric(NA)]
	outls <- which(bins$ld_region10kb_perc0714>0.99)
	for(i in 1:length(outls)){
		bins$ld_region10kb_cluster0714[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$dist0714[outls[i+1]]>=width | is.na(bins$dist0714[outls[i+1]])) indx <- indx+1
		}
	}

	indx <- 1
	bins[,ld_region10kb_clusterCAN := as.numeric(NA)]
	outls <- which(bins$ld_region10kb_percCAN>0.99)
	for(i in 1:length(outls)){
		bins$ld_region10kb_clusterCAN[outls[i]] <- indx
		if(i < length(outls)){
			if(bins$distCAN[outls[i+1]]>=width | is.na(bins$distCAN[outls[i+1]])) indx <- indx+1
		}
	}


# save region means data
filenm <- paste('analysis/ld_change_region_', windsz, '.rds', sep='')
filenm
saveRDS(bins, file=filenm)

###########################################
# plot LD change from regions
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
