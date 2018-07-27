# compare Tajima's D near outlier loci and near non-outlier loci


######################
# calculate Tajima's D in windows
# to run on a cod node!
######################
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
outl <- fread("zcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz") # 270MB file
#	outl

# read in Tajima's D data on cod node (note: includes inversions. Need to keep track of bp windows with any nan windows trimmed out)
	# for 50 bp windows
#bp <- 50
#dat14 <- fread('analysis/LOF_S_14.Tajima.D')
#dat11 <- fread('analysis/LOF_S_11.Tajima.D')
#dat07 <- fread('analysis/LOF_07.Tajima.D')
#dat40 <- fread('analysis/CAN40.Tajima.D')
#datMod <- fread('analysis/CANMod.Tajima.D')

	# for 100 bp windows
#bp <- 100
#dat14 <- fread('analysis/LOF_S_14_100bp.Tajima.D')
#dat11 <- fread('analysis/LOF_S_11_100bp.Tajima.D')
#dat07 <- fread('analysis/LOF_07_100bp.Tajima.D')
#dat40 <- fread('analysis/CAN40_100bp.Tajima.D')
#datMod <- fread('analysis/CANMod_100bp.Tajima.D')

	# for 1000 bp windows (also masked)
bp <- 1000
dat14 <- fread('analysis/LOF_S_14_1000bp_mask.Tajima.D')
dat11 <- fread('analysis/LOF_S_11_1000bp_mask.Tajima.D')
dat07 <- fread('analysis/LOF_07_1000bp_mask.Tajima.D')
dat40 <- fread('analysis/CAN40_1000bp_mask.Tajima.D')
datMod <- fread('analysis/CANMod_1000bp_mask.Tajima.D')

# add nearest XX bp bin start
# dat14[,sort(unique(substr(BIN_START, nchar(BIN_START)-1, nchar(BIN_START))))]
if(bp == 50) {# (either 00, 20, 40, 60, or 80)
	outl[,BIN_START := floor(POS/20)*20] 
}
if(bp == 100) {# (all 00)
	outl[,BIN_START := floor(POS/100)*100] 
}
if(bp == 1000) {# (all 000)
	outl[,BIN_START := floor(POS/1000)*1000] 
}

# remove inversions and unplaced (also remove unmappable?)
dat14 <- dat14[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]
dat11 <- dat11[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]
dat07 <- dat07[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]
dat40 <- dat40[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]
datMod <- datMod[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]
outl <- outl[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]


# merge in outlier information to dat
# use q<0.3 for more loci in calculations
nrow(dat14) # 920397 (50bp) 721949 (100bp) 136809 (1000bp)
nrow(dat11) # 825344
nrow(dat07) # 819717
nrow(dat40) # 625016
nrow(datMod) #634704

dat14 <- merge(dat14, outl[,.(CHROM, POS, BIN_START, kmer25, dpLofFlag, outlier=p.Lof.adj3<0.3)], by.x=c('CHROM', 'BIN_START'), by.y=c('CHROM', 'BIN_START'), all=TRUE)

dat11 <- merge(dat11, outl[,.(CHROM, POS, BIN_START, kmer25, dpLofFlag, outlier=p.Lof.adj3<0.3)], by.x=c('CHROM', 'BIN_START'), by.y=c('CHROM', 'BIN_START'), all=TRUE)
	
dat07 <- merge(dat07, outl[,.(CHROM, POS, BIN_START, kmer25, dpLofFlag, outlier=p.Lof.adj3<0.3)], by.x=c('CHROM', 'BIN_START'), by.y=c('CHROM', 'BIN_START'), all=TRUE)

dat40 <- merge(dat40, outl[,.(CHROM, POS, BIN_START, kmer25, dpCanFlag, outlier=p.Can.adj3<0.3)], by.x=c('CHROM', 'BIN_START'), by.y=c('CHROM', 'BIN_START'), all=TRUE)

datMod <- merge(datMod, outl[,.(CHROM, POS, BIN_START, kmer25, dpCanFlag, outlier=p.Can.adj3<0.3)], by.x=c('CHROM', 'BIN_START'), by.y=c('CHROM', 'BIN_START'), all=TRUE)

nrow(dat14) # 1370545
nrow(dat11) # 2006940
nrow(dat07) # 2002440
nrow(dat40) # 1851140
nrow(datMod) # 1858579

# remove low-quality rows (SNPs that fail a filter)
# note: could also remove all bins with SNPs that fail a filter
dat14 <- dat14[kmer25==1 & dpLofFlag==TRUE,]
dat11 <- dat11[kmer25==1 & dpLofFlag==TRUE,]
dat07 <- dat07[kmer25==1 & dpLofFlag==TRUE,]
dat40 <- dat40[kmer25==1 & dpCanFlag==TRUE,]
datMod <- datMod[kmer25==1 & dpCanFlag==TRUE,]

nrow(dat14) # 330524
nrow(dat11) 
nrow(dat07) 
nrow(dat40) # 254047
nrow(datMod) 

# how many outlier comparisons?
dat14[,summary(outlier)] # 60
dat11[,summary(outlier)]
dat07[,summary(outlier)]
dat40[,summary(outlier)] # 5
datMod[,summary(outlier)]

# set up 1000 not-outlier loci
	dat14[,nearoutlier:=FALSE]
	dat11[,nearoutlier:=FALSE]
	dat07[,nearoutlier:=FALSE]
	dat40[,nearoutlier:=FALSE]
	datMod[,nearoutlier:=FALSE]

	# for Norway
	dat14[outlier==TRUE, nearoutlier:=TRUE] # list of loci that are outliers
	dat11[outlier==TRUE, nearoutlier:=TRUE] # list of loci that are outliers
	dat07[outlier==TRUE, nearoutlier:=TRUE] # list of loci that are outliers
	for(i in dat14[,which(outlier)]) { # add loci that are <= 5000 bp from outliers
		dat14[CHROM==dat14$CHROM[i] & abs(POS - dat14$POS[i]) <= 5000, nearoutlier:=TRUE]
		dat11[CHROM==dat14$CHROM[i] & abs(POS - dat14$POS[i]) <= 5000, nearoutlier:=TRUE]
		dat07[CHROM==dat14$CHROM[i] & abs(POS - dat14$POS[i]) <= 5000, nearoutlier:=TRUE]
	}
	allloci <- dat14[nearoutlier==FALSE, .(CHROM, POS)] # list of all loci that aren't near outliers
	allloci$dist = NA
	allloci$dist[2:nrow(allloci)] <- allloci$POS[2:nrow(allloci)] - allloci$POS[1:(nrow(allloci)-1)] # calculate distance to next locus to the left
		dim(allloci)
	allloci <- allloci[dist > 5000 & POS > 5000,] # trim out loci <5000 bp from start of CHROM or from another locus (doesn't check end of CHROM...)
		dim(allloci)
	notoutl <- allloci[sample(.N, 1000)] # sample 1000 loci to be the "not outliers"
	notoutl[,notoutlier:=TRUE] # label these as the 1000 not outliers
	notoutl[,dist:=NULL] # remove dist column
	dat14 <- merge(dat14, notoutl, by=c('CHROM', 'POS'), all.x=TRUE) # label the not-outliers
	dat11 <- merge(dat11, notoutl, by=c('CHROM', 'POS'), all.x=TRUE)
	dat07 <- merge(dat07, notoutl, by=c('CHROM', 'POS'), all.x=TRUE)
	dat14[is.na(notoutlier), notoutlier:=FALSE]	
	dat11[is.na(notoutlier), notoutlier:=FALSE]	
	dat07[is.na(notoutlier), notoutlier:=FALSE]	


	# for Canada
	dat40[outlier==TRUE, nearoutlier:=TRUE] # list of loci that are outliers
	datMod[outlier==TRUE, nearoutlier:=TRUE] # list of loci that are outliers
	for(i in dat40[,which(outlier)]) { # add loci that are <= 5000 bp from outliers
		dat40[CHROM==dat40$CHROM[i] & abs(POS - dat40$POS[i]) <= 5000, nearoutlier:=TRUE]
		datMod[CHROM==dat40$CHROM[i] & abs(POS - dat40$POS[i]) <= 5000, nearoutlier:=TRUE]
	}
	allloci <- dat40[nearoutlier==FALSE, .(CHROM, POS)] # list of all loci that aren't near outliers
	allloci$dist = NA
	allloci$dist[2:nrow(allloci)] <- allloci$POS[2:nrow(allloci)] - allloci$POS[1:(nrow(allloci)-1)] # calculate distance to next locus to the left
		dim(allloci)
	allloci <- allloci[dist > 5000 & POS > 5000,] # trim out loci <5000 bp from start of CHROM or from another locus (doesn't check end of CHROM...)
		dim(allloci)
	notoutl <- allloci[sample(.N, 1000)] # sample 1000 loci to be the "not outliers"
	notoutl[,notoutlier:=TRUE] # label these as the 1000 not outliers
	notoutl[,dist:=NULL] # remove dist column
	dat40 <- merge(dat40, notoutl, by=c('CHROM', 'POS'), all.x=TRUE) # label the not-outliers
	datMod <- merge(datMod, notoutl, by=c('CHROM', 'POS'), all.x=TRUE)
	dat40[is.na(notoutlier), notoutlier:=FALSE]	
	datMod[is.na(notoutlier), notoutlier:=FALSE]	


	dat14[,summary(notoutlier)] # 1000
	dat11[,summary(notoutlier)]
	dat07[,summary(notoutlier)]
	dat40[,summary(notoutlier)]
	datMod[,summary(notoutlier)]

# calculate distance in units of Tajima's D bin sizes

	if('distclass' %in% colnames(dat14)) dat14[,':=' (distclass=NULL, disttype=NULL)]
	dat14[, distclass:=as.numeric(NA)]
	dat14[, disttype:=as.factor(NA)]
	levels(dat14$disttype) <- c('outlier', 'notoutlier')
	if('distclass' %in% colnames(dat11)) dat11[,':=' (distclass=NULL, disttype=NULL)]
	dat11[, distclass:=as.numeric(NA)]
	dat11[, disttype:=as.factor(NA)]
	levels(dat11$disttype) <- c('outlier', 'notoutlier')
	if('distclass' %in% colnames(dat07)) dat07[,':=' (distclass=NULL, disttype=NULL)]
	dat07[, distclass:=as.numeric(NA)]
	dat07[, disttype:=as.factor(NA)]
	levels(dat07$disttype) <- c('outlier', 'notoutlier')

	if('distclass' %in% colnames(dat40)) dat40[,':=' (distclass=NULL, disttype=NULL)]
	dat40[, distclass:=as.numeric(NA)]
	dat40[, disttype:=as.factor(NA)]
	levels(dat40$disttype) <- c('outlier', 'notoutlier')
	if('distclass' %in% colnames(datMod)) datMod[,':=' (distclass=NULL, disttype=NULL)]
	datMod[, distclass:=as.numeric(NA)]
	datMod[, disttype:=as.factor(NA)]
	levels(datMod$disttype) <- c('outlier', 'notoutlier')

	# from each outlier
	# note this will write over any outliers whose 5000 bp overlap any earlier values
		# Norway
	dim(dat14) # 330524
	dim(dat11)
	dim(dat07)
	for(i in which(dat14$outlier)){
		dat14[CHROM==dat14$CHROM[i] & abs(BIN_START-dat14$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat14$BIN_START[i]), disttype='outlier')]
		dat11[CHROM==dat14$CHROM[i] & abs(BIN_START-dat14$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat14$BIN_START[i]), disttype='outlier')]
		dat07[CHROM==dat14$CHROM[i] & abs(BIN_START-dat14$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat14$BIN_START[i]), disttype='outlier')]
	}

		# Canada
	dim(dat40) # 254047
	dim(datMod)
	for(i in which(dat40$outlier)){
		dat40[CHROM==dat40$CHROM[i] & abs(BIN_START-dat40$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat40$BIN_START[i]), disttype='outlier')]
		datMod[CHROM==dat40$CHROM[i] & abs(BIN_START-dat40$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat40$BIN_START[i]), disttype='outlier')]
	}
	

	# and from each nonoutlier
		# Norway
	for(i in which(dat14$notoutlier)){
		dat14[CHROM==dat14$CHROM[i] & abs(BIN_START-dat14$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat14$BIN_START[i]), disttype='notoutlier')]
		dat11[CHROM==dat14$CHROM[i] & abs(BIN_START-dat14$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat14$BIN_START[i]), disttype='notoutlier')]
		dat07[CHROM==dat14$CHROM[i] & abs(BIN_START-dat14$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat14$BIN_START[i]), disttype='notoutlier')]
	}

		# Canada
	for(i in which(dat40$notoutlier)){
		dat40[CHROM==dat40$CHROM[i] & abs(BIN_START-dat40$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat40$BIN_START[i]), disttype='notoutlier')]
		datMod[CHROM==dat40$CHROM[i] & abs(BIN_START-dat40$BIN_START[i]) <= 5000, ':=' (distclass=abs(BIN_START-dat40$BIN_START[i]), disttype='notoutlier')]
	}


# trim to unique rows (distance bins near outliers or not outliers)
nrow(dat14) # 330524
nrow(dat11)
nrow(dat07)
nrow(dat40) # 254047
nrow(datMod)

dat14t <- unique(dat14[!is.na(TajimaD)], by=c('CHROM', 'BIN_START', 'TajimaD', 'distclass', 'disttype'))
dat11t <- unique(dat11[!is.na(TajimaD)], by=c('CHROM', 'BIN_START', 'TajimaD', 'distclass', 'disttype'))
dat07t <- unique(dat07[!is.na(TajimaD)], by=c('CHROM', 'BIN_START', 'TajimaD', 'distclass', 'disttype'))
dat40t <- unique(dat40[!is.na(TajimaD)], by=c('CHROM', 'BIN_START', 'TajimaD', 'distclass', 'disttype'))
datModt <- unique(datMod[!is.na(TajimaD)], by=c('CHROM', 'BIN_START', 'TajimaD', 'distclass', 'disttype'))

nrow(dat14t) # 294275 (50bp) 221475 (100bp) 133329 (1000bp)
nrow(dat11t)
nrow(dat07t)
nrow(dat40t) # 229878 (50bp) 182012 (100bp) 120795 (1000bp)
nrow(datModt)


# calculate averages within bins (like Bario et al. 2016 eLife)
bins14 <- dat14t[!is.na(TajimaD),.(Dave=mean(TajimaD), Dsd=sd(TajimaD), Dl95=quantile(TajimaD,probs=0.025), Du95=quantile(TajimaD,probs=0.975), Dl75=quantile(TajimaD,probs=0.125), Du75=quantile(TajimaD,probs=0.875), Dn=.N), by=list(disttype, distclass)] # average within distance classes for outlier or not outlier (but not neither)

bins11 <- dat11t[!is.na(TajimaD),.(Dave=mean(TajimaD), Dsd=sd(TajimaD), Dl95=quantile(TajimaD,probs=0.025), Du95=quantile(TajimaD,probs=0.975), Dl75=quantile(TajimaD,probs=0.125), Du75=quantile(TajimaD,probs=0.875), Dn=.N), by=list(disttype, distclass)]

bins07 <- dat07t[!is.na(TajimaD),.(Dave=mean(TajimaD), Dsd=sd(TajimaD), Dl95=quantile(TajimaD,probs=0.025), Du95=quantile(TajimaD,probs=0.975), Dl75=quantile(TajimaD,probs=0.125), Du75=quantile(TajimaD,probs=0.875), Dn=.N), by=list(disttype, distclass)]

bins40 <- dat40t[!is.na(TajimaD),.(Dave=mean(TajimaD), Dsd=sd(TajimaD), Dl95=quantile(TajimaD,probs=0.025), Du95=quantile(TajimaD,probs=0.975), Dl75=quantile(TajimaD,probs=0.125), Du75=quantile(TajimaD,probs=0.875), Dn=.N), by=list(disttype, distclass)]

binsMod <- datModt[!is.na(TajimaD),.(Dave=mean(TajimaD), Dsd=sd(TajimaD), Dl95=quantile(TajimaD,probs=0.025), Du95=quantile(TajimaD,probs=0.975), Dl75=quantile(TajimaD,probs=0.125), Du75=quantile(TajimaD,probs=0.875), Dn=.N), by=list(disttype, distclass)]


bins14[,Dse:=Dsd/Dn]
bins11[,Dse:=Dsd/Dn]
bins07[,Dse:=Dsd/Dn]
bins40[,Dse:=Dsd/Dn]
binsMod[,Dse:=Dsd/Dn]

setkey(bins14, disttype, distclass)
setkey(bins11, disttype, distclass)
setkey(bins07, disttype, distclass)
setkey(bins40, disttype, distclass)
setkey(binsMod, disttype, distclass)

# examine
bins14[,.(disttype,distclass,Dave, Dsd, Dn, Dse)]
bins40[,.(disttype,distclass,Dave, Dsd, Dn, Dse)]
binsMod[,.(disttype,distclass,Dave, Dsd, Dn, Dse)]

# label by population
bins14[,pop:='LOF_S_14']
bins11[,pop:='LOF_S_11']
bins07[,pop:='LOF_07']
bins40[,pop:='CAN40']
binsMod[,pop:='CANMod']

# merge
bins <- rbind(bins07, bins11, bins14, bins40, binsMod)

# write out
write.csv(bins, file=paste('analysis/TajimaD_outliers_', bp, 'bp.csv', sep=''), row.names=FALSE)

###############
# plot
# to run on macbook
###############
require(RColorBrewer)
require(data.table)
require(mgcv)

# read in data
#bins <- fread('analysis/TajimaD_outliers.csv', ); bp=50
bins <- fread('analysis/TajimaD_outliers_100bp.csv', ); bp=100
bins <- fread('analysis/TajimaD_outliers_1000bp.csv', ); bp=1000

# calculate smooth fits
smooth14 <- bins[pop=='LOF_S_14' & disttype=='outlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smooth11 <- bins[pop=='LOF_S_11' & disttype=='outlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smooth07 <- bins[pop=='LOF_07' & disttype=='outlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smooth14n <- bins[pop=='LOF_S_14' & disttype=='notoutlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smooth11n <- bins[pop=='LOF_S_11' & disttype=='notoutlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smooth07n <- bins[pop=='LOF_07' & disttype=='notoutlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smooth40 <- bins[pop=='CAN40' & disttype=='outlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smoothMod <- bins[pop=='CANMod' & disttype=='outlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smooth40n <- bins[pop=='CAN40' & disttype=='notoutlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]
smoothModn <- bins[pop=='CANMod' & disttype=='notoutlier', predict(loess(Dave ~ I(distclass)), se=TRUE)]

# plotting parameters
cols <- brewer.pal(5, 'Set1')
colstrans <- paste(cols, '55', sep='')
cex=0.5
lwd=2
ylims <- bins[!is.na(disttype),range(c(Dave-1.96*Dse, Dave+1.96*Dse), na.rm=TRUE)]



### plot data with loess smoothers
quartz(width=6, height=3)
# pdf(width=6, height=3, file=paste('figures/TajimasD_decay_outliers_', bp, 'bp_loess.pdf', sep=''))
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Lof
		# points outlier
bins[pop=='LOF_S_14' & disttype=='outlier',plot(distclass, Dave, ylim=ylims, type='p', xlab='Distance (bp)', ylab='Pi per site', cex=cex, main='Lof', log='', col=colstrans[1], pch=16)]
bins[pop=='LOF_S_11' & disttype=='outlier',points(distclass, Dave, cex=cex, col=colstrans[2], pch=16)]
bins[pop=='LOF_07' & disttype=='outlier',points(distclass, Dave, cex=cex, col=colstrans[3], pch=16)]
		# points notoutlier
bins[pop=='LOF_S_14' & disttype=='notoutlier',points(distclass, Dave, cex=cex, col=colstrans[1], pch=1)]
bins[pop=='LOF_S_11' & disttype=='notoutlier',points(distclass, Dave, cex=cex, col=colstrans[2], pch=1)]
bins[pop=='LOF_07' & disttype=='notoutlier',points(distclass, Dave, cex=cex, col=colstrans[3], pch=1)]

		# polygon areas outlier
bins[pop=='LOF_S_14' & disttype=='outlier',polygon(c(distclass, rev(distclass)), c(smooth14$fit - smooth14$se.fit, rev(smooth14$fit + smooth14$se.fit)), cex=cex, col=colstrans[1], border=NA)]
bins[pop=='LOF_S_11' & disttype=='outlier',polygon(c(distclass, rev(distclass)), c(smooth11$fit - smooth11$se.fit, rev(smooth11$fit + smooth11$se.fit)), cex=cex, col=colstrans[2], border=NA)]
bins[pop=='LOF_07' & disttype=='outlier',polygon(c(distclass, rev(distclass)), c(smooth07$fit - smooth07$se.fit, rev(smooth07$fit + smooth07$se.fit)), cex=cex, col=colstrans[3], border=NA)]
		# polygon areas notoutlier
bins[pop=='LOF_S_14' & disttype=='notoutlier',polygon(c(distclass, rev(distclass)), c(smooth14n$fit - smooth14n$se.fit, rev(smooth14n$fit + smooth14n$se.fit)), cex=cex, col=colstrans[1], border=NA)]
bins[pop=='LOF_S_11' & disttype=='notoutlier',polygon(c(distclass, rev(distclass)), c(smooth11n$fit - smooth11n$se.fit, rev(smooth11n$fit + smooth11n$se.fit)), cex=cex, col=colstrans[2], border=NA)]
bins[pop=='LOF_07' & disttype=='notoutlier',polygon(c(distclass, rev(distclass)), c(smooth07n$fit - smooth07n$se.fit, rev(smooth07n$fit + smooth07n$se.fit)), cex=cex, col=colstrans[3], border=NA)]

		# lines outlier
bins[pop=='LOF_S_14' & disttype=='outlier',lines(distclass, smooth14$fit, cex=cex, col=cols[1], lwd=lwd)]
bins[pop=='LOF_S_11' & disttype=='outlier',lines(distclass, smooth11$fit, cex=cex, col=cols[2], lwd=lwd)]
bins[pop=='LOF_07' & disttype=='outlier',lines(distclass, smooth07$fit, cex=cex, col=cols[3], lwd=lwd)]
		# lines notoutlier
bins[pop=='LOF_S_14' & disttype=='notoutlier',lines(distclass, smooth14n$fit, cex=cex, col=cols[1], lwd=lwd, lty=3)]
bins[pop=='LOF_S_11' & disttype=='notoutlier',lines(distclass, smooth11n$fit, cex=cex, col=cols[2], lwd=lwd, lty=3)]
bins[pop=='LOF_07' & disttype=='notoutlier',lines(distclass, smooth07n$fit, cex=cex, col=cols[3], lwd=lwd, lty=3)]



legend('bottomleft', legend=c('LOF_S_14', 'LOF_S_11', 'LOF_07', 'Outlier', 'Not outlier'), col=c(cols[1:3], 'black', 'black'), lty=c(rep(1,4), 3), bty='n', cex=0.4)

	# Canada
		# points outlier
bins[pop=='CANMod' & disttype=='outlier',plot(distclass, Dave, ylim=ylims, xlim=c(1,5000), type='p', xlab='Distance (bp)', ylab='Pi per site', cex=cex, main='Can', log='', col=colstrans[4], pch=16)]
bins[pop=='CAN40' & disttype=='outlier',points(distclass, Dave, cex=cex, col=colstrans[5], pch=16)]
		# points notoutlier
bins[pop=='CANMod' & disttype=='notoutlier',points(distclass, Dave, cex=cex, col=colstrans[4], pch=1)]
bins[pop=='CAN40' & disttype=='notoutlier',points(distclass, Dave, cex=cex, col=colstrans[5], pch=1)]

		# polygon areas outlier
bins[pop=='CANMod' & disttype=='outlier',polygon(c(distclass, rev(distclass)), c(smoothMod$fit - smoothMod$se.fit, rev(smoothMod$fit + smoothMod$se.fit)), cex=cex, col=colstrans[4], border=NA)]
bins[pop=='CAN40' & disttype=='outlier',polygon(c(distclass, rev(distclass)), c(smooth40$fit - smooth40$se.fit, rev(smooth40$fit + smooth40$se.fit)), cex=cex, col=colstrans[5], border=NA)]
		# polygon areas notoutlier
bins[pop=='CANMod' & disttype=='notoutlier',polygon(c(distclass, rev(distclass)), c(smoothModn$fit - smoothModn$se.fit, rev(smoothModn$fit + smoothModn$se.fit)), cex=cex, col=colstrans[4], border=NA)]
bins[pop=='CAN40' & disttype=='notoutlier',polygon(c(distclass, rev(distclass)), c(smooth40n$fit - smooth40n$se.fit, rev(smooth40n$fit + smooth40n$se.fit)), cex=cex, col=colstrans[5], border=NA)]

		# lines outlier
bins[pop=='CANMod' & disttype=='outlier',lines(distclass, smoothMod$fit, cex=cex, col=cols[4], lwd=lwd)]
bins[pop=='CAN40' & disttype=='outlier',lines(distclass, smooth40$fit, cex=cex, col=cols[5], lwd=lwd)]
		# lines notoutlier
bins[pop=='CANMod' & disttype=='notoutlier',lines(distclass, smoothModn$fit, cex=cex, col=cols[4], lwd=lwd, lty=3)]
bins[pop=='CAN40' & disttype=='notoutlier',lines(distclass, smooth40n$fit, cex=cex, col=cols[5], lwd=lwd, lty=3)]

legend('topright', legend=c('CANMod', 'CAN40', 'Outlier', 'Not outlier'), col=c(cols[4:5], 'black', 'black'), lty=c(1,1,1,3), bty='n', cex=0.4)

dev.off()





### plot binned data
quartz(width=6, height=3)
# pdf(width=6, height=3, file=paste('figures/TajimasD_decay_outliers_', bp, 'bp.pdf', sep=''))
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Lof
bins[pop=='LOF_S_14' & disttype=='outlier',plot(distclass, Dave, ylim=ylims, type='l', xlab='Distance (bp)', ylab="Tajima's D", cex=cex, main='Lof', log='', col=cols[1], lwd=lwd)]
bins[pop=='LOF_S_11' & disttype=='outlier',lines(distclass+100, Dave, type='l', cex=cex, col=cols[2], lwd=lwd)]
bins[pop=='LOF_07' & disttype=='outlier',lines(distclass+200, Dave, type='l', cex=cex, col=cols[3], lwd=lwd)]

bins[pop=='LOF_S_14' & disttype=='nooutlier',lines(distclass, Dave, type='l', cex=cex, col=cols[1], lty=3, lwd=lwd)]
bins[pop=='LOF_S_11' & disttype=='nooutlier',lines(distclass+100, Dave, type='l', cex=cex, col=cols[2], lty=3, lwd=lwd)]
bins[pop=='LOF_07' & disttype=='notoutlier',lines(distclass+200, Dave, type='l', cex=cex, col=cols[3], lty=3, lwd=lwd)]

bins[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[1]), by=distclass] # 95% CI on the mean.
bins[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass+100, distclass+100), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & disttype=='outlier',lines(c(distclass+200, distclass+200), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[3]), by=distclass]

bins[pop=='LOF_S_14' & disttype=='notoutlier',lines(c(distclass, distclass), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[1]), by=distclass]
bins[pop=='LOF_S_11' & disttype=='notoutlier',lines(c(distclass+100, distclass+100), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & disttype=='notoutlier',lines(c(distclass+200, distclass+200), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[3]), by=distclass]

legend('topright', legend=c('LOF_S_14', 'LOF_S_11', 'LOF_07', 'Outlier', 'Not outlier'), col=c(cols[1:3], 'black', 'black'), lty=c(rep(1,4), 3), bty='n', cex=0.4)

	# Canada
bins[pop=='CANMod' & disttype=='outlier',plot(distclass, Dave, ylim=ylims, type='l', xlab='Distance (bp)', ylab="Tajima's D", cex=cex, main='Can', log='', col=cols[4], lwd=lwd)]
bins[pop=='CAN40' & disttype=='outlier',lines(distclass+100, Dave, type='l', cex=cex, col=cols[5], lwd=lwd)]

bins[pop=='CANMod' & disttype=='notoutlier',lines(distclass, Dave, type='l', cex=cex, col=cols[4], lty=3, lwd=lwd)]
bins[pop=='CAN40' & disttype=='notoutlier',lines(distclass+100, Dave, type='l', cex=cex, col=cols[5], lty=3, lwd=lwd)]

bins[pop=='CANMod' & disttype=='outlier',lines(c(distclass, distclass), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[4]), by=distclass]
bins[pop=='CAN40' & disttype=='outlier',lines(c(distclass+100, distclass+100), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[5]), by=distclass]
bins[pop=='CANMod' & disttype=='notoutlier',lines(c(distclass, distclass), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[4]), by=distclass]
bins[pop=='CAN40' & disttype=='notoutlier',lines(c(distclass+100, distclass+100), c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols[5]), by=distclass]

legend('topright', legend=c('CANMod', 'CAN40', 'Outlier', 'Not outlier'), col=c(cols[4:5], 'black', 'black'), lty=c(1,1,1,3), bty='n', cex=0.4)

dev.off()