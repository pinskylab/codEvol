# compare LD near outlier loci and near non-outlier loci

# flags
#outliertype <- 'bypop' # use q3.Lof071114 and q3.Can < 0.3 to define outlier loci
#outliertype <- 'combinedpop' # use q3.comb071114Can < 0.3
outliertype <- 'union' # use q3.comb071114Can < 0.2 | q3.Lof071114 < 0.2 | q3.Can < 0.2

dpfilter <- TRUE
#dpfilter <- FALSE

mapfilter <- TRUE
#mapfilter <- FALSE

if((mapfilter==FALSE & dpfilter==TRUE) | ((dpfilter==FALSE | mapfilter==FALSE) & outliertype != 'union')) warning('NO CODE YET for this flag combination. THINK CAREFULLY')

######################
# calculate LD decay
# to run on a cod node!
######################
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
outl <- fread("zcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz")
#	outl

# read in LD data on cod node (each file about 600MB)
dat14 <- fread('zcat analysis/LOF_S_14.geno.ld.gz')
dat11 <- fread('zcat analysis/LOF_S_11.geno.ld.gz')
dat07 <- fread('zcat analysis/LOF_07.geno.ld.gz')
dat40 <- fread('zcat analysis/CAN40.geno.ld.gz')
datMod <- fread('zcat analysis/CANMod.geno.ld.gz')

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

# define outlier information
if(outliertype == 'bypop' & dpfilter==TRUE & mapfilter==TRUE){
	outl[, outlierLof:=q3.Lof071114<0.2]
	outl[, outlierCan:=q3.Can<0.2]
}

if(outliertype == 'combinedpop' & dpfilter==TRUE & mapfilter==TRUE){
	outl[, outlierLof:=q3.comb071114Can<0.2]
	outl[, outlierCan:=q3.comb071114Can<0.2]
}

if(outliertype == 'union' & dpfilter==TRUE & mapfilter==TRUE){
	outl[, outlierLof:=(q3.comb071114Can<0.2 | q3.Lof071114<0.2)]
	outl[, outlierCan:=(q3.comb071114Can<0.2 | q3.Can<0.2)]
}

if(outliertype == 'union' & dpfilter==FALSE & mapfilter==TRUE){
	print('No depth filter!')

	# have to do new FDR adjustments if we relax these filters
	outl[kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q4.comb071114Can := p.adjust(p.comb071114Can, method='fdr')]
	outl[kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q4.Lof071114 := p.adjust(pLof071114, method='fdr')]
	outl[kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q4.Can := p.adjust(pCan, method='fdr')]

	outl[, outlierLof:=(q4.comb071114Can<0.2 | q4.Lof071114<0.2)]
	outl[, outlierCan:=(q4.comb071114Can<0.2 | q4.Can<0.2)]
}

if(outliertype == 'union' & dpfilter==FALSE & mapfilter==FALSE){
	print('No depth or map filter!')

	# have to do new FDR adjustments if we relax these filters
	outl[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q4.comb071114Can := p.adjust(p.comb071114Can, method='fdr')]
	outl[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q4.Lof071114 := p.adjust(pLof071114, method='fdr')]
	outl[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q4.Can := p.adjust(pCan, method='fdr')]

	outl[, outlierLof:=(q4.comb071114Can<0.2 | q4.Lof071114<0.2)]
	outl[, outlierCan:=(q4.comb071114Can<0.2 | q4.Can<0.2)]
}

outl[,sum(outlierLof, na.rm=TRUE)] # 50 (union<0.2) XX (union<0.2 no dpfilter)
outl[,sum(outlierCan, na.rm=TRUE)] # 45 (union<0.2) 112 (union<0.2 no dpfilter)


# merge in outlier information
nrow(dat14) # 6148211
nrow(dat11)
nrow(dat07)
nrow(dat40) # 6148211
nrow(datMod)

dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlierPOS1=outlierLof)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE) # on POS1
dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlierPOS2=outlierLof)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE) # on POS2

dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlierPOS1=outlierLof)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlierPOS2=outlierLof)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	
dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlierPOS1=outlierLof)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlierPOS2=outlierLof)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)

dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlierPOS1=outlierCan)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlierPOS2=outlierCan)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)

datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlierPOS1=outlierCan)], by.x=c('CHR', 'POS1'), by.y=c('CHROM', 'POS'), all.x=TRUE)
datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlierPOS2=outlierCan)], by.x=c('CHR', 'POS2'), by.y=c('CHROM', 'POS'), all.x=TRUE)


nrow(dat14) # 6148211
nrow(dat11)
nrow(dat07)
nrow(dat40) # 6148211
nrow(datMod)

# combine information from POS1 and POS2 quality filters and outliers
# has to pass kmer25 filter or dp filter for both loci, but only one locus needs to be an outlier to be flagged
# use filters within populations
if(dpfilter==TRUE & mapfilter==TRUE){
	dat14[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpLofFlag.x&dpLofFlag.y), outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat11[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpLofFlag.x&dpLofFlag.y), outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat07[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpLofFlag.x&dpLofFlag.y), outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat40[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpCanFlag.x&dpCanFlag.y), outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	datMod[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=(dpCanFlag.x&dpCanFlag.y), outlier=(outlierPOS1==1)|(outlierPOS2==1))]
}
if(dpfilter==FALSE & mapfilter==TRUE){
	print('No depth filter!')
	dat14[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat11[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat07[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat40[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	datMod[,':='(kmer25=(kmer25.x==1)&(kmer25.y==1), dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
}
if(dpfilter==FALSE & mapfilter==FALSE){
	print('No depth or map filter!')
	dat14[,':='(kmer25=TRUE, dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat11[,':='(kmer25=TRUE, dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat07[,':='(kmer25=TRUE, dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	dat40[,':='(kmer25=TRUE, dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
	datMod[,':='(kmer25=TRUE, dpFlag=TRUE, outlier=(outlierPOS1==1)|(outlierPOS2==1))]
}

# remove unneeded columns
dat14[,':='(kmer25.x=NULL, kmer25.y=NULL, dpLofFlag.x=NULL, dpLofFlag.y=NULL)]
dat11[,':='(kmer25.x=NULL, kmer25.y=NULL, dpLofFlag.x=NULL, dpLofFlag.y=NULL)]
dat07[,':='(kmer25.x=NULL, kmer25.y=NULL, dpLofFlag.x=NULL, dpLofFlag.y=NULL)]
dat40[,':='(kmer25.x=NULL, kmer25.y=NULL, dpCanFlag.x=NULL, dpCanFlag.y=NULL)]
datMod[,':='(kmer25.x=NULL, kmer25.y=NULL, dpCanFlag.x=NULL, dpCanFlag.y=NULL)]

# remove low-quality rows (comparisons that involve at least one locus that fails a filter)
nrow(dat14) # 6148211
nrow(dat11)
nrow(dat07)
nrow(dat40)
nrow(datMod)

dat14 <- dat14[kmer25==TRUE & dpFlag==TRUE,]
dat11 <- dat11[kmer25==TRUE & dpFlag==TRUE,]
dat07 <- dat07[kmer25==TRUE & dpFlag==TRUE,]
dat40 <- dat40[kmer25==TRUE & dpFlag==TRUE,]
datMod <- datMod[kmer25==TRUE & dpFlag==TRUE,]

nrow(dat14) # 2673672 (p.comb, union) xx (union no dpfilter)  (union no dpfilter no mapfilter)
nrow(dat11)
nrow(dat07)
nrow(dat40) # 1763983 (p.comb, union) xx (union no dpfilter)  (union no dpfilter no mapfilter)
nrow(datMod)

# how many outlier comparisons?
dat14[,sum(outlier, na.rm=TRUE)] # xx (p.comb<0.2), x (p.comb<0.3), 665 (union<0.2) xx (union<0.2 no dpfilter)  (union<0.2 no dpfilter no mapfilter)
dat11[,sum(outlier, na.rm=TRUE)]
dat07[,sum(outlier, na.rm=TRUE)]
dat40[,sum(outlier, na.rm=TRUE)] # xx (p.comb<0.2), xx (p.comb<0.3), 394 (union<0.3) xx (union<0.2 no dpfilter)  (union<0.2 no dpfilter no mapfilter)
datMod[,sum(outlier, na.rm=TRUE)]

# set up 1000 not-outlier loci
dat14[,notoutlier:=FALSE]
dat11[,notoutlier:=FALSE]
dat07[,notoutlier:=FALSE]
dat40[,notoutlier:=FALSE]
datMod[,notoutlier:=FALSE]

exclude <- dat14[outlier==TRUE, unique(c(paste(CHR, POS1), paste(CHR, POS2)))] # list of loci that are outliers or are <5000 bp from outliers (latter happens because we calculated r2 between loci up to 5000bp apart)
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

dat14[,summary(notoutlier)] # 14302 (union<0.2), but varies per run and by table
dat11[,summary(notoutlier)]
dat07[,summary(notoutlier)]
dat40[,summary(notoutlier)] # 11782 (union<0.2), but varies per run and by table
datMod[,summary(notoutlier)]

# remove comparisons that aren't near outliers or chosen not-outliers
# also remove NA values for r2
dat14 <- dat14[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
dat11 <- dat11[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
dat07 <- dat07[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
dat40 <- dat40[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]
datMod <- datMod[(outlier==TRUE | notoutlier==TRUE) & !is.na(r2),]

nrow(dat14) #  (p.comb<0.2)  (p.comb<0.3) 12233 (union<0.2): vary by run
nrow(dat11) #  
nrow(dat07) #  
nrow(dat40) #  (p.comb<0.2)  (p.comb<0.3) 9029 (union<0.2)
nrow(datMod) # 

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
if(outliertype=='union'){
	stp1 <- 30 # step size for small distances
	thresh <- 60 # use stp2 above this distance
	stp2 <- 1000
}

dat14[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat14[dist>=thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

dat11[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat11[dist>=thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

dat07[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat07[dist>=thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

dat40[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat40[dist>=thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

datMod[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
datMod[dist>=thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

# mark outlier or not
dat14[,outlierclass := '']
dat14[outlier==TRUE,outlierclass := 'outlier']
dat14[notoutlier==TRUE,outlierclass := 'notoutlier']

dat11[,outlierclass := '']
dat11[outlier==TRUE,outlierclass := 'outlier']
dat11[notoutlier==TRUE,outlierclass := 'notoutlier']

dat07[,outlierclass := '']
dat07[outlier==TRUE,outlierclass := 'outlier']
dat07[notoutlier==TRUE,outlierclass := 'notoutlier']

dat40[,outlierclass := '']
dat40[outlier==TRUE,outlierclass := 'outlier']
dat40[notoutlier==TRUE,outlierclass := 'notoutlier']

datMod[,outlierclass := '']
datMod[outlier==TRUE,outlierclass := 'outlier']
datMod[notoutlier==TRUE,outlierclass := 'notoutlier']

# calculate averages within bins (like Bario et al. 2016 eLife)
bins14 <- dat14[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlierclass, distclass)] # average within distance classes for outlier or not outlier (but not neither)

bins11 <- dat11[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlierclass, distclass)]

bins07 <- dat07[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlierclass, distclass)]

bins40 <- dat40[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlierclass, distclass)]

binsMod <- datMod[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=list(outlierclass, distclass)]

bins14[,r2se:=r2sd/r2n]
bins11[,r2se:=r2sd/r2n]
bins07[,r2se:=r2sd/r2n]
bins40[,r2se:=r2sd/r2n]
binsMod[,r2se:=r2sd/r2n]

setkey(bins14, outlierclass, distclass)
setkey(bins11, outlierclass, distclass)
setkey(bins07, outlierclass, distclass)
setkey(bins40, outlierclass, distclass)
setkey(binsMod, outlierclass, distclass)

# examine
bins14[,.(outlierclass,distclass,r2ave, r2se)]
bins11[,.(outlierclass,distclass,r2ave, r2se)]
bins07[,.(outlierclass,distclass,r2ave, r2se)]
binsMod[,.(outlierclass,distclass,r2ave, r2se, r2n)]
bins40[,.(outlierclass,distclass,r2ave, r2se, r2n)] # check for low sample size and small bins

# label by population
bins14[,pop:='LOF_S_14']
bins11[,pop:='LOF_S_11']
bins07[,pop:='LOF_07']
bins40[,pop:='CAN40']
binsMod[,pop:='CANMod']

# merge
bins <- rbind(bins07, bins11, bins14, bins40, binsMod)

# write out average LD around outliers and not outliers
if(dpfilter==TRUE & mapfilter==TRUE) outfile <- paste('analysis/ld_outliers_', outliertype, '.csv', sep='')
if(dpfilter==FALSE & mapfilter==TRUE) outfile <- paste('analysis/ld_outliers_', outliertype, '_nodpfilter.csv', sep='')
if(dpfilter==FALSE & mapfilter==FALSE) outfile <- paste('analysis/ld_outliers_', outliertype, '_nodpfilter_nomapfilter.csv', sep='')
outfile
write.csv(bins, file=outfile)

# write out all LD calculations around outliers
dat14[,pop:='LOF_S_14']
dat11[,pop:='LOF_S_11']
dat07[,pop:='LOF_07']
dat40[,pop:='CAN40']
datMod[,pop:='CANMod']

outlLD <- rbind(dat14[outlierclass=='outlier',], dat11[outlierclass=='outlier',], dat07[outlierclass=='outlier',], dat40[outlierclass=='outlier',], datMod[outlierclass=='outlier',]) # merge into one file
	nrow(outlLD) # 1472
	
outlLD2 <- outlLD[outlierPOS1==TRUE, .(CHR, POSoutl=POS1, POScomp=POS2, N_INDV, r2, dist, kmer25, dpFlag, pop)] # reformat so outlier position is always in the same column
outlLD2 <- rbind(outlLD2, outlLD[outlierPOS2==TRUE, .(CHR, POSoutl=POS2, POScomp=POS1, N_INDV, r2, dist, kmer25, dpFlag, pop)])
	nrow(outlLD2) # 1576, since some loci were close to each other

outlLD2[, dist:=POScomp - POSoutl] # pos or neg

setkey(outlLD2, pop, CHR, POSoutl, dist) # sort
	
if(dpfilter==TRUE & mapfilter==TRUE) outfile2 <- paste('analysis/ld_outliers_allcalcs_', outliertype, '.csv', sep='')
if(dpfilter==FALSE & mapfilter==TRUE) outfile2 <- paste('analysis/ld_outliers_allcalcs_', outliertype, '_nodpfilter.csv', sep='')
if(dpfilter==FALSE & mapfilter==FALSE) outfile2 <- paste('analysis/ld_outliers_allcalcs_', outliertype, '_nodpfilter_nomapfilter.csv', sep='')
outfile2
write.csv(outlLD2, file=outfile2, row.names=FALSE)



#################################################
# plot average LD around outliers/not outliers
# to run on macbook
#################################################
require(RColorBrewer)
require(data.table)

#bins <- fread('analysis/ld_outliers.csv'); outliertype<-'bypop'
#bins <- fread('analysis/ld_outliers_combinedpop.csv'); outliertype <- 'combinedpop'
#bins <- fread('analysis/ld_outliers_union.csv'); outliertype <- 'union'
bins <- fread('analysis/ld_outliers_union_nodpfilter.csv'); outliertype <- 'union_nodpfilter'
cols <- brewer.pal(5, 'Set1')
cex=0.5
lwd=2

# plot binned data
quartz(width=6, height=3)
# pdf(width=6, height=3, file=paste('figures/ld_decay_outliers_', outliertype, '.pdf', sep=''))
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Lof
bins[pop=='LOF_S_14' & outlierclass=='outlier',plot(distclass, r2ave, ylim=c(0,1), type='l', xlab='Distance (bp)', ylab='Average correlation (r2)', cex=cex, main='Lof', log='x', col=cols[1], lwd=lwd)]
bins[pop=='LOF_S_11' & outlierclass=='outlier',lines(distclass*1.1, r2ave, type='l', cex=cex, col=cols[2], lwd=lwd)]
bins[pop=='LOF_07' & outlierclass=='outlier',lines(distclass*1.2, r2ave, type='l', cex=cex, col=cols[3], lwd=lwd)]

bins[pop=='LOF_S_14' & outlierclass=='notoutlier',lines(distclass*1.3, r2ave, type='l', cex=cex, col=cols[1], lty=3, lwd=lwd)]
bins[pop=='LOF_S_11' & outlierclass=='notoutlier',lines(distclass*1.4, r2ave, type='l', cex=cex, col=cols[2], lty=3, lwd=lwd)]
bins[pop=='LOF_07' & outlierclass=='notoutlier',lines(distclass*1.5, r2ave, type='l', cex=cex, col=cols[3], lty=3, lwd=lwd)]

bins[pop=='LOF_S_14' & outlierclass=='outlier',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[1]), by=distclass] # 95% CI on the mean.
bins[pop=='LOF_S_11' & outlierclass=='outlier',lines(c(distclass*1.1, distclass*1.1), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & outlierclass=='outlier',lines(c(distclass*1.2, distclass*1.2), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[3]), by=distclass]

bins[pop=='LOF_S_14' & outlierclass=='notoutlier',lines(c(distclass*1.3, distclass*1.3), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[1]), by=distclass]
bins[pop=='LOF_S_11' & outlierclass=='notoutlier',lines(c(distclass*1.4, distclass*1.4), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & outlierclass=='notoutlier',lines(c(distclass*1.5, distclass*1.5), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[3]), by=distclass]

legend('topright', legend=c('LOF_S_14', 'LOF_S_11', 'LOF_07', 'Outlier', 'Not outlier'), col=c(cols[1:3], 'black', 'black'), lty=c(rep(1,4), 3), bty='n', cex=0.5)

	# Canada
bins[pop=='CANMod' & outlierclass=='outlier',plot(distclass, r2ave, ylim=c(0,1), type='l', xlab='Distance (bp)', ylab='Average correlation (r2)', cex=cex, main='Can', log='x', col=cols[4], lwd=lwd)]
bins[pop=='CAN40' & outlierclass=='outlier',lines(distclass*1.1, r2ave, type='l', cex=cex, col=cols[5], lwd=lwd)]

bins[pop=='CANMod' & outlierclass=='notoutlier',lines(distclass*1.2, r2ave, type='l', cex=cex, col=cols[4], lty=3, lwd=lwd)]
bins[pop=='CAN40' & outlierclass=='notoutlier',lines(distclass*1.3, r2ave, type='l', cex=cex, col=cols[5], lty=3, lwd=lwd)]

bins[pop=='CANMod' & outlierclass=='outlier',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[4]), by=distclass]
bins[pop=='CAN40' & outlierclass=='outlier',lines(c(distclass*1.1, distclass*1.1), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[5]), by=distclass]
bins[pop=='CANMod' & outlierclass=='notoutlier',lines(c(distclass*1.2, distclass*1.2), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[4]), by=distclass]
bins[pop=='CAN40' & outlierclass=='notoutlier',lines(c(distclass*1.3, distclass*1.3), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[5]), by=distclass]

legend('topright', legend=c('CANMod', 'CAN40', 'Outlier', 'Not outlier'), col=c(cols[4:5], 'black', 'black'), lty=c(1,1,1,3), bty='n', cex=0.5)

dev.off()



#################################################
# plot LD calcs around outliers/not outliers
# to run on macbook
#################################################
require(RColorBrewer)
require(data.table)

outlLD <- fread('analysis/ld_outliers_allcalcs_union.csv'); outliertype <- 'union'

cols <- brewer.pal(5, 'Set1')
names(cols) <- c('LOF_07', 'LOF_S_11', 'LOF_S_14', 'CAN40', 'CANMod')
cex=0.5
lwd=2

# make list of outlier loci
locs <- unique(outlLD[,.(CHR, POSoutl)])
setkey(locs, CHR, POSoutl)

# plot binned data
quartz(width=8.5, height=11)
# pdf(width=8.5, height=11, file=paste('figures/ld_decay_byoutlier_', outliertype, '.pdf', sep=''))
par(mfrow=c(6,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

for(i in 1:nrow(locs)){
#for(i in 1:12){
	pops <- outlLD[CHR==locs$CHR[i] & POSoutl == locs$POSoutl[i], sort(unique(pop))]
	outlLD[CHR==locs$CHR[i] & POSoutl == locs$POSoutl[i] & pop==pops[1],plot(dist, r2, type='o', xlim=c(-5000,5000), ylim=c(0,1), main=paste(locs$CHR[i], locs$POSoutl[i]), col=cols[[pops[1]]])]
	if(length(pops)>1){
		for(j in 2:length(pops)){
			outlLD[CHR==locs$CHR[i] & POSoutl == locs$POSoutl[i] & pop==pops[j],points(dist, r2, col=cols[pops[j]], type='o')]
		}
	}

	if(i%%12 == 1) legend('topright', legend=names(cols), col=cols, pch=1)

}


dev.off()