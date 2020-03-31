## Plot FST vs. genome position
## Compare across different data sets

###########################
# load functions
###########################
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){ # not cod node
	require(data.table)
	require(plyr)
  require(ggplot2)
  require(RColorBrewer)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
}


#####################
# read and prep in data
#####################

# sliding windows from ANGSD (all sites)
# header is missing the fst column, so have to skip and make our own
can <- fread('analysis/Can_40.Can_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof11 <- fread('analysis/Lof_07.Lof_11.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof14 <- fread('analysis/Lof_07.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof1114 <- fread('analysis/Lof_11.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 

# sliding windows from ANGSD (GATK sites)
cangatk <- fread('analysis/Can_40.Can_14.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof11gatk <- fread('analysis/Lof_07.Lof_11.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof14gatk <- fread('analysis/Lof_07.Lof_14.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof1114gatk <- fread('analysis/Lof_11.Lof_14.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 


# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- can[,.(len=max(midPos)), by=chr]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, chr)

# merge nucleotide position into the frequency files
setkey(can, chr)
can <- can[chrmax[, .(chr, start)], ]
can[, posgen := midPos + start]
can[,start := NULL]

setkey(lof11, chr)
lof11 <- lof11[chrmax[, .(chr, start)], ]
lof11[, posgen := midPos + start]
lof11[,start := NULL]

setkey(lof14, chr)
lof14 <- lof14[chrmax[, .(chr, start)], ]
lof14[, posgen := midPos + start]
lof14[,start := NULL]

setkey(lof1114, chr)
lof1114 <- lof1114[chrmax[, .(chr, start)], ]
lof1114[, posgen := midPos + start]
lof1114[,start := NULL]

setkey(cangatk, chr)
cangatk <- cangatk[chrmax[, .(chr, start)], ]
cangatk[, posgen := midPos + start]
cangatk[,start := NULL]

setkey(lof11gatk, chr)
lof11gatk <- lof11gatk[chrmax[, .(chr, start)], ]
lof11gatk[, posgen := midPos + start]
lof11gatk[,start := NULL]

setkey(lof14gatk, chr)
lof14gatk <- lof14gatk[chrmax[, .(chr, start)], ]
lof14gatk[, posgen := midPos + start]
lof14gatk[,start := NULL]

setkey(lof1114gatk, chr)
lof1114gatk <- lof1114gatk[chrmax[, .(chr, start)], ]
lof1114gatk[, posgen := midPos + start]
lof1114gatk[,start := NULL]


# combine the datasets to look at correlations across them
can[, pop := 'can']
lof11[, pop := 'lof11']
lof14[, pop := 'lof14']
lof1114[, pop := 'lof1114']

cangatk[, pop := 'can']
lof11gatk[, pop := 'lof11']
lof14gatk[, pop := 'lof14']
lof1114gatk[, pop := 'lof1114']

dat <- rbind(lof11, lof14, lof1114, can)
datgatk <- rbind(lof11gatk, lof14gatk, lof1114gatk, cangatk)

nrow(dat)
nrow(datgatk)

dat
datgatk

dat[, pop := factor(pop, levels = c('can', 'lof11', 'lof14', 'lof1114'))]
datgatk[, pop := factor(pop, levels = c('can', 'lof11', 'lof14', 'lof1114'))]


##############
# plots
##############

# plot fst vs. #snps
ggplot(dat, aes(Nsites, fst, color = pop)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  scale_x_log10()

ggplot(datgatk, aes(Nsites, fst, color = pop)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  scale_x_log10()


# plot fst vs. position
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]

p1 <- ggplot(dat, aes(posgen, fst, color = chr)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols)
p1
ggsave(plot = p1, device = 'png', filename = 'figures/fst_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)


p2 <- ggplot(datgatk, aes(posgen, fst, color = chr)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols)
p2
ggsave(plot = p2, device = 'png', filename = 'figures/fst_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)


############
# to be updated below here

# calculate percentiles (like Oziolor et al. 2019 Science)
ecdf0711 <- ecdf(datmean$ABS_DIFF_0711)
ecdf0714 <- ecdf(datmean$ABS_DIFF_0714)
ecdf1114 <- ecdf(datmean$ABS_DIFF_1114)
ecdfCan <- ecdf(datmean$ABS_DIFF_Can)
datmean[,perc0711 := ecdf0711(ABS_DIFF_0711)]
datmean[,perc0714 := ecdf0714(ABS_DIFF_0714)]
datmean[,perc1114 := ecdf1114(ABS_DIFF_1114)]
datmean[,percCAN := ecdfCan(ABS_DIFF_Can)]

	# examine
	datmean[perc0711>0.99 & perc0714>0.99 & percCAN>0.99,]
	datmean[perc0711>0.999 & perc0714>0.999 & percCAN>0.999,]
	
	datmean[plot(perc0711, perc0714, cex=0.5, col='#00000005')]
	datmean[plot(perc0711, percCAN, cex=0.5, col='#00000005')]
	
# group the 1% outliers
	# find distance among outlier regions
	setkey(datmean, CHROM, POSmid)
	datmean[,dist0711 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	datmean[,dist0714 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	datmean[,dist1114 := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	datmean[,distCAN := as.numeric(NA)] # nearest neighbor at an earlier outlier (measure only to left)
	for(i in 1:length(chroms)){
		temp0711 <- c(NA, datmean[CHROM==chroms[i] & perc0711>0.99, POSmid])
		temp0714 <- c(NA, datmean[CHROM==chroms[i] & perc0714>0.99, POSmid])
		temp1114 <- c(NA, datmean[CHROM==chroms[i] & perc1114>0.99, POSmid])
		tempCAN <- c(NA, datmean[CHROM==chroms[i] & percCAN>0.99, POSmid])

		temp0711 <- temp0711[-length(temp0711)]
		temp0714 <- temp0714[-length(temp0714)]
		temp1114 <- temp1114[-length(temp1114)]
		tempCAN <- tempCAN[-length(tempCAN)]

		datmean[CHROM==chroms[i] & perc0711>0.99, dist0711:= POSmid - temp0711]
		datmean[CHROM==chroms[i] & perc0714>0.99, dist0714:= POSmid - temp0714]
		datmean[CHROM==chroms[i] & perc1114>0.99, dist1114:= POSmid - temp1114]
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
	datmean[,cluster1114 := as.numeric(NA)]
	outls <- which(datmean$perc1114>0.99)
	for(i in 1:length(outls)){
		datmean$cluster1114[outls[i]] <- indx
		if(i < length(outls)){
			if(datmean$dist1114[outls[i+1]]>=width | is.na(datmean$dist1114[outls[i+1]])) indx <- indx+1
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

	# examine
	datmean[,range(table(cluster0711))] # all clusters only 1 window wide
	datmean[,range(table(cluster0714))] # all clusters only 1 window wide
	datmean[,range(table(cluster1114))] # all clusters only 1 window wide
	datmean[,range(table(clusterCAN))] # all clusters only 1 window wide
	

# save region means data
filenm <- paste('analysis/Frequency_table_ABS_DIFF_runmean', windsz, '.rds', sep='')
filenm
saveRDS(datmean, file=filenm)


# add region outlier stats to locus-by-locus frequency change file
dat2 <- merge(dat, datmean[, .(CHROM, POSmid, BIN_START, perc0711, perc0714, perc1114, percCAN, cluster0711, cluster0714, cluster1114, clusterCAN)], by.x=c('CHROM', 'POSmid1'), by.y=c('CHROM', 'POSmid'))
	dim(dat)
	dim(dat2)
	dim(datmean)
setnames(dat2, c('perc0711', 'perc0714', 'perc1114', 'percCAN', 'cluster0711', 'cluster0714', 'cluster1114', 'clusterCAN'), c('region_perc0711', 'region_perc0714', 'region_perc1114', 'region_percCAN', 'region_cluster0711', 'region_cluster0714', 'region_cluster1114', 'region_clusterCAN'))

	# examine outliers
	dat2[region_perc0711>0.99 & region_perc0714>0.99 & region_percCAN>0.99,]

# save locus data with region data appended
filenm2 <- paste('analysis/Frequency_table_ABS_DIFF_w_region_outliers_', windsz, '.rds', sep='')
filenm2
saveRDS(dat2, file=filenm2)


###########################
## rank the outlier regions (see Oziolor et al. 2019 Science)
##########################
datmean <- readRDS('analysis/Frequency_table_ABS_DIFF_runmean1e4.rds'); width='1e4'

outlreg0711 = datmean[,.(CHROM=unique(CHROM), POSmin=min(POSmid), POSmax=max(POSmid), score=sum(ABS_DIFF_0711)), by=cluster0711]
outlreg0714 = datmean[,.(CHROM=unique(CHROM), POSmin=min(POSmid), POSmax=max(POSmid), score=sum(ABS_DIFF_0714)), by=cluster0714]
outlregCAN = datmean[,.(CHROM=unique(CHROM), POSmin=min(POSmid), POSmax=max(POSmid), score=sum(ABS_DIFF_Can)), by=clusterCAN]

setkey(outlreg0711, score)
setkey(outlreg0714, score)
setkey(outlregCAN, score)

tail(outlreg0711[CHROM != 'LG01'], 5)
tail(outlreg0714[CHROM != 'LG01'], 5)
tail(outlregCAN[CHROM != 'LG01'], 5)

