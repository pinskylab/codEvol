# compare LD near outlier loci and near non-outlier loci

# flags
#outliertype <- 'bypop' # use q3.Lof071114 and q3.Can < 0.3 to define outlier loci
outliertype <- 'combinedpop' # use q3.comb071114Can < 0.3

######################
# calculate LD decay
# to run on a cod node!
######################
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
outl <- fread("zcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz")
#	outl

# read in LD data on cod node
dat14 <- fread('zcat analysis/LOF_S_14.geno.ld')
dat11 <- fread('zcat analysis/LOF_S_11.geno.ld')
dat07 <- fread('zcat analysis/LOF_07.geno.ld')
dat40 <- fread('zcat analysis/CAN40.geno.ld')
datMod <- fread('zcat analysis/CANMod.geno.ld')

# remove ^ from a column name
setnames(dat14, 5, 'r2')
setnames(dat11, 5, 'r2')
setnames(dat07, 5, 'r2')
setnames(dat40, 5, 'r2')
setnames(datMod, 5, 'r2')

# calculate distance between loci
dat14[,dist:=abs(POS2-POS1)]
dat11[,dist:=abs(POS2-POS1)]
dat07[,dist:=abs(POS2-POS1)]
dat40[,dist:=abs(POS2-POS1)]
datMod[,dist:=abs(POS2-POS1)]

# merge in outlier information
# use q<0.3 for more loci in LD calculations
nrow(dat14) # 19706634
nrow(dat11)
nrow(dat07)
nrow(dat40) # 9918703
nrow(datMod)

if(outliertype == 'bypop'){
	dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.Lof071114<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE) # on POS1
	dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.Lof071114<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE) # on POS2

	dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.Lof071114<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.Lof071114<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	
	dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.Lof071114<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.Lof071114<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)

	dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.Can<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.Can<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)

	datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.Can<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.Can<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)
}

if(outliertype == 'combinedpop'){
	dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE) # on POS1
	dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE) # on POS2

	dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	
	dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)

	dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)

	datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=q3.comb071114Can<0.3)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)
}

nrow(dat14) # 19706634
nrow(dat11)
nrow(dat07)
nrow(dat40) # 9918703
nrow(datMod)

# combine information from POS1 and POS2 quality filters and outliers
# has to pass kmer25 filter or dp filter for both loci, but only one locus needs to be an outlier to be flagged
# use filters within populations
dat14[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpLofFlag.x&dpLofFlag.y), outlier=(outlier.x==1)|(outlier.y==1))]
dat11[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpLofFlag.x&dpLofFlag.y), outlier=(outlier.x==1)|(outlier.y==1))]
dat07[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpLofFlag.x&dpLofFlag.y), outlier=(outlier.x==1)|(outlier.y==1))]
dat40[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpCanFlag.x&dpCanFlag.y), outlier=(outlier.x==1)|(outlier.y==1))]
datMod[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpCanFlag.x&dpCanFlag.y), outlier=(outlier.x==1)|(outlier.y==1))]

# remove unneeded columns
dat14[,':='(kmer25.x=NULL, kmer25.y=NULL, dpLofFlag.x=NULL, dpLofFlag.y=NULL, outlier.x=NULL, outlier.y=NULL)]
dat11[,':='(kmer25.x=NULL, kmer25.y=NULL, dpLofFlag.x=NULL, dpLofFlag.y=NULL, outlier.x=NULL, outlier.y=NULL)]
dat07[,':='(kmer25.x=NULL, kmer25.y=NULL, dpLofFlag.x=NULL, dpLofFlag.y=NULL, outlier.x=NULL, outlier.y=NULL)]
dat40[,':='(kmer25.x=NULL, kmer25.y=NULL, dpCanFlag.x=NULL, dpCanFlag.y=NULL, outlier.x=NULL, outlier.y=NULL)]
datMod[,':='(kmer25.x=NULL, kmer25.y=NULL, dpCanFlag.x=NULL, dpCanFlag.y=NULL, outlier.x=NULL, outlier.y=NULL)]

# remove low-quality rows (comparisons that involve at least one locus that fails a filter)
nrow(dat14)
nrow(dat11)
nrow(dat07)
nrow(dat40)
nrow(datMod)

dat14 <- dat14[kmer25==TRUE & dpFlag==TRUE,]
dat11 <- dat11[kmer25==TRUE & dpFlag==TRUE,]
dat07 <- dat07[kmer25==TRUE & dpFlag==TRUE,]
dat40 <- dat40[kmer25==TRUE & dpFlag==TRUE,]
datMod <- datMod[kmer25==TRUE & dpFlag==TRUE,]

nrow(dat14) # 2137914 (p.comb)
nrow(dat11)
nrow(dat07)
nrow(dat40) # 1180832 (p.comb)
nrow(datMod)

# how many outlier comparisons?
dat14[,sum(outlier, na.rm=TRUE)] # 330 (p.comb<0.2), 410 (p.comb<0.3)
dat11[,sum(outlier, na.rm=TRUE)]
dat07[,sum(outlier, na.rm=TRUE)]
dat40[,sum(outlier, na.rm=TRUE)] # 157 (p.comb<0.2), 218 (p.comb<0.3)
datMod[,sum(outlier, na.rm=TRUE)]

# set up 1000 not-outlier loci
dat14[,notoutlier:=FALSE]
dat11[,notoutlier:=FALSE]
dat07[,notoutlier:=FALSE]
dat40[,notoutlier:=FALSE]
datMod[,notoutlier:=FALSE]

exclude <- dat14[outlier==TRUE, unique(c(paste(CHR, POS1), paste(CHR, POS2)))] # list of loci that are outliers or are <5000 bp from outliers (latter because we calculated r2 between loci up to 5000bp apart)
allloci <- dat14[, unique(c(paste(CHR, POS1), paste(CHR, POS2)))] # list of all loci
notoutl <- sample(setdiff(allloci, exclude), 1000) # sample 1000 loci to be the "not outliers"
dat14[(paste(CHR, POS1) %in% notoutl) | (paste(CHR, POS2) %in% notoutl), notoutlier:=TRUE] # label the not-outliers

exclude <- dat11[outlier==TRUE, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
allloci <- dat11[, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
notoutl <- sample(setdiff(allloci, exclude), 1000)
dat11[(paste(CHR, POS1) %in% notoutl) | (paste(CHR, POS2) %in% notoutl), notoutlier:=TRUE] # label the not-outliers

exclude <- dat07[outlier==TRUE, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
allloci <- dat07[, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
notoutl <- sample(setdiff(allloci, exclude), 1000)
dat07[(paste(CHR, POS1) %in% notoutl) | (paste(CHR, POS2) %in% notoutl), notoutlier:=TRUE] # label the not-outliers

exclude <- dat40[outlier==TRUE, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
allloci <- dat40[, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
notoutl <- sample(setdiff(allloci, exclude), 1000)
dat40[(paste(CHR, POS1) %in% notoutl) | (paste(CHR, POS2) %in% notoutl), notoutlier:=TRUE] # label the not-outliers

exclude <- datMod[outlier==TRUE, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
allloci <- datMod[, unique(c(paste(CHR, POS1), paste(CHR, POS2)))]
notoutl <- sample(setdiff(allloci, exclude), 1000)
datMod[(paste(CHR, POS1) %in% notoutl) | (paste(CHR, POS2) %in% notoutl), notoutlier:=TRUE] # label the not-outliers

dat14[,summary(notoutlier)]
dat11[,summary(notoutlier)]
dat07[,summary(notoutlier)]
dat40[,summary(notoutlier)]
datMod[,summary(notoutlier)]

# remove comparisons that aren't near outliers or chosen not-outliers
# also remove NA values for r2
dat14 <- dat14[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
dat11 <- dat11[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
dat07 <- dat07[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
dat40 <- dat40[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
datMod <- datMod[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]

nrow(dat14) # 11306 (p.comb<0.2) 11248 (p.comb<0.3): these numbers seem off
nrow(dat11) # 10661 (p.comb<0.2) 11211 (p.comb<0.3)
nrow(dat07) # 10257 (p.comb<0.2) 11228 (p.comb<0.3)
nrow(dat40) # 7491 (p.comb<0.2) 7710 (p.comb<0.3)
nrow(datMod) # 8120 (p.comb<0.2) 8547 (p.comb<0.3)

# examine outliers
dat07[outlier==TRUE, .(CHR, POS1, POS2, dist, N_INDV, r2, kmer25, dpFlag)]
dat40[outlier==TRUE, .(CHR, POS1, POS2, dist, N_INDV, r2, kmer25, dpFlag)]

# set up distance bins
if(outliertype=='bypop'){
	stp1 <- 10 # step size for small distances
	thresh <- 50 # use stp2 above this distance
	stp2 <- 1000
}
if(outliertype=='combinedpop'){
	stp1 <- 50
	thresh <- 100
	stp2 <- 1000
}

dat14[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat14[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

dat11[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat11[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

dat07[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat07[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

dat40[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat40[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

datMod[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
datMod[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class


# calculate averages within bins (like Bario et al. 2016 eLife)
bins14 <- dat14[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlier, distclass)] # average within distance classes for outlier or not outlier (but not neither)

bins11 <- dat11[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlier, distclass)]

bins07 <- dat07[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlier, distclass)]

bins40 <- dat40[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlier, distclass)]

binsMod <- datMod[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlier, distclass)]

bins14[,r2se:=r2sd/r2n]
bins11[,r2se:=r2sd/r2n]
bins07[,r2se:=r2sd/r2n]
bins40[,r2se:=r2sd/r2n]
binsMod[,r2se:=r2sd/r2n]

setkey(bins14, outlier, distclass)
setkey(bins11, outlier, distclass)
setkey(bins07, outlier, distclass)
setkey(bins40, outlier, distclass)
setkey(binsMod, outlier, distclass)

# examine
bins14[,.(outlier,distclass,r2ave, r2se)]
bins11[,.(outlier,distclass,r2ave, r2se)]
bins07[,.(outlier,distclass,r2ave, r2se)]
binsMod[,.(outlier,distclass,r2ave, r2se)]
bins40[,.(outlier,distclass,r2ave, r2se)]

# label by population
bins14[,pop:='LOF_S_14']
bins11[,pop:='LOF_S_11']
bins07[,pop:='LOF_07']
bins40[,pop:='CAN40']
binsMod[,pop:='CANMod']

# merge
bins <- rbind(bins07, bins11, bins14, bins40, binsMod)

# write out
write.csv(bins, file=paste('analysis/ld_outliers_', outliertype, '.csv', sep=''))

###############
# plot
# to run on macbook
###############
require(RColorBrewer)
require(data.table)

#bins <- fread('analysis/ld_outliers.csv'); outliertype<-'bypop'
bins <- fread('analysis/ld_outliers_combinedpop.csv'); outliertype <- 'combinedpop'
cols <- brewer.pal(5, 'Set1')
cex=0.5
lwd=2

# plot binned data
quartz(width=6, height=3)
# pdf(width=6, height=3, file=paste('figures/ld_decay_outliers_', outliertype, '.pdf', sep=''))
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Lof
bins[pop=='LOF_S_14' & outlier==TRUE,plot(distclass, r2ave, ylim=c(0,1), type='l', xlab='Distance (bp)', ylab='Average correlation (r2)', cex=cex, main='Lof', log='x', col=cols[1], lwd=lwd)]
bins[pop=='LOF_S_11' & outlier==TRUE,lines(distclass*1.1, r2ave, type='l', cex=cex, col=cols[2], lwd=lwd)]
bins[pop=='LOF_07' & outlier==TRUE,lines(distclass*1.2, r2ave, type='l', cex=cex, col=cols[3], lwd=lwd)]

bins[pop=='LOF_S_14' & outlier==FALSE,lines(distclass*1.3, r2ave, type='l', cex=cex, col=cols[1], lty=3, lwd=lwd)]
bins[pop=='LOF_S_11' & outlier==FALSE,lines(distclass*1.4, r2ave, type='l', cex=cex, col=cols[2], lty=3, lwd=lwd)]
bins[pop=='LOF_07' & outlier==FALSE,lines(distclass*1.5, r2ave, type='l', cex=cex, col=cols[3], lty=3, lwd=lwd)]

bins[pop=='LOF_S_14' & outlier==TRUE,lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[1]), by=distclass] # 95% CI on the mean.
bins[pop=='LOF_S_11' & outlier==TRUE,lines(c(distclass*1.1, distclass*1.1), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & outlier==TRUE,lines(c(distclass*1.2, distclass*1.2), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[3]), by=distclass]

bins[pop=='LOF_S_14' & outlier==FALSE,lines(c(distclass*1.3, distclass*1.3), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[1]), by=distclass]
bins[pop=='LOF_S_11' & outlier==FALSE,lines(c(distclass*1.4, distclass*1.4), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & outlier==FALSE,lines(c(distclass*1.5, distclass*1.5), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[3]), by=distclass]

legend('topright', legend=c('LOF_S_14', 'LOF_S_11', 'LOF_07', 'Outlier', 'Not outlier'), col=c(cols[1:3], 'black', 'black'), lty=c(rep(1,4), 3), bty='n', cex=0.5)

	# Canada
bins[pop=='CANMod' & outlier==TRUE,plot(distclass, r2ave, ylim=c(0,1), type='l', xlab='Distance (bp)', ylab='Average correlation (r2)', cex=cex, main='Can', log='x', col=cols[4], lwd=lwd)]
bins[pop=='CAN40' & outlier==TRUE,lines(distclass*1.1, r2ave, type='l', cex=cex, col=cols[5], lwd=lwd)]

bins[pop=='CANMod' & outlier==FALSE,lines(distclass*1.2, r2ave, type='l', cex=cex, col=cols[4], lty=3, lwd=lwd)]
bins[pop=='CAN40' & outlier==FALSE,lines(distclass*1.3, r2ave, type='l', cex=cex, col=cols[5], lty=3, lwd=lwd)]

bins[pop=='CANMod' & outlier==TRUE,lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[4]), by=distclass]
bins[pop=='CAN40' & outlier==TRUE,lines(c(distclass*1.1, distclass*1.1), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[5]), by=distclass]
bins[pop=='CANMod' & outlier==FALSE,lines(c(distclass*1.2, distclass*1.2), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[4]), by=distclass]
bins[pop=='CAN40' & outlier==FALSE,lines(c(distclass*1.3, distclass*1.3), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[5]), by=distclass]

legend('topright', legend=c('CANMod', 'CAN40', 'Outlier', 'Not outlier'), col=c(cols[4:5], 'black', 'black'), lty=c(1,1,1,3), bty='n', cex=0.5)

dev.off()