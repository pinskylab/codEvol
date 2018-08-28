# compare snp density near outlier loci and near non-outlier loci
# account for which loci can be called as SNPs/not SNPs
# use unfiltered data (no kmer or depth filter)

#outliertype <- 'bypop' # use q3.Lof071114 and q3.Can < 0.3 to define outlier loci
#outliertype <- 'combinedpop' # use q3.comb071114Can < 0.3
outliertype <- 'union' # use q3.comb071114Can < 0.2 | q3.Lof071114 < 0.2 | q3.Can < 0.2

dpfilter <- FALSE
mapfilter <- FALSE

######################
# calculate pi in windows
# to run on a cod node!
######################
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
outl <- fread("zcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz") # 270MB file
#	outl

# read in pi data on cod node as way to find polymorphic SNPs
# note: files include inversions. have not been trimmed for kmer25 or depth
dat14 <- fread('analysis/LOF_S_14.sites.pi')
dat11 <- fread('analysis/LOF_S_11.sites.pi')
dat07 <- fread('analysis/LOF_07.sites.pi')
dat40 <- fread('analysis/CAN40.sites.pi')
datMod <- fread('analysis/CANMod.sites.pi')

# read in all callable sites information
if(mapfilter) allsites <- fread('zcat all_sites_data_29_06_18/AllSites.kept.sites.kmer25.gz') # gzipped file
if(!mapfilter) allsites <- fread('zcat all_sites_data_29_06_18/AllSites.kept.sites.gz') # gzipped file 6.5GB
setkey(allsites, CHROM, POS)

# define outlier information
# use q<0.3 for more loci in LD calculations (except union Norway)
if(outliertype == 'bypop' & dpfilter==TRUE & mapfilter==TRUE){
	outl[, outlierLof:=q3.Lof071114<0.3]
	outl[, outlierCan:=q3.Can<0.3]
}

if(outliertype == 'combinedpop' & dpfilter==TRUE & mapfilter==TRUE){
	outl[, outlierLof:=q3.comb071114Can<0.3]
	outl[, outlierCan:=q3.comb071114Can<0.3]
}

if(outliertype == 'union'){
	# use the outliers identified with the filters on, whether or not filters turned on now
	outl[, outlierLof:=(q3.comb071114Can<0.2 | q3.Lof071114<0.2)]
	outl[, outlierCan:=(q3.comb071114Can<0.2 | q3.Can<0.2)]
}


outl[,sum(outlierLof, na.rm=TRUE)] # 59 (union<0.2) 182 (union<0.2 no dpfilter)
outl[,sum(outlierCan, na.rm=TRUE)] # 27 (union<0.2) 94 (union<0.2 no dpfilter)


# merge in outlier information to dat
nrow(dat14) # 1402283
nrow(dat11)
nrow(dat07)
nrow(dat40) # 1018417
nrow(datMod)

dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=outlierLof)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)
dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=outlierLof)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)
dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=outlierLof)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)
dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=outlierCan)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)
datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=outlierCan)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)

nrow(dat14)
nrow(dat11)
nrow(dat07)
nrow(dat40)
nrow(datMod)

# remove low-quality rows (SNPs that fail a filter)
nrow(dat14)
nrow(dat11)
nrow(dat07)
nrow(dat40)
nrow(datMod)

if(dpfilter==TRUE & mapfilter==TRUE){
	dat14 <- dat14[kmer25==TRUE & dpLofFlag==TRUE,]
	dat11 <- dat11[kmer25==TRUE & dpLofFlag==TRUE,]
	dat07 <- dat07[kmer25==TRUE & dpLofFlag==TRUE,]
	dat40 <- dat40[kmer25==TRUE & dpCanFlag==TRUE,]
	datMod <- datMod[kmer25==TRUE & dpCanFlag==TRUE,]
}
if(dpfilter==FALSE & mapfilter==TRUE){
	print('No depth filter!')
	dat14 <- dat14[kmer25==TRUE,]
	dat11 <- dat11[kmer25==TRUE,]
	dat07 <- dat07[kmer25==TRUE,]
	dat40 <- dat40[kmer25==TRUE,]
	datMod <- datMod[kmer25==TRUE,]
}
if(dpfilter==FALSE & mapfilter==FALSE){
	print('No depth or map filter!')
}

nrow(dat14) # 408771 (both filters) 429195 (no dpfilter) 1402283 (no depth or map filter)
nrow(dat11)
nrow(dat07)
nrow(dat40) # 313787 (both filters) 329301 (no dpfilter) 1018417 (no depth or map filter)
nrow(datMod)

# how many outlier comparisons?
# WHY SO MANY NAs IN THE OUTLIER COLUMN? APPARENTLY WEREN'T ANY FOR THE LD CALCULATIONS
dat14[,summary(outlier)] # should be same as in outl
dat11[,summary(outlier)]
dat07[,summary(outlier)]
dat40[,summary(outlier)]
datMod[,summary(outlier)]

# set up not-outlier loci
	dat14[,notoutlier:=FALSE]
	dat11[,notoutlier:=FALSE]
	dat07[,notoutlier:=FALSE]
	dat40[,notoutlier:=FALSE]
	datMod[,notoutlier:=FALSE]

	nnotoutl <- 100 # number of not-outliers to choose

	# for Norway
	exclude <- dat14[outlier==TRUE, .(CHROM, POS)] # list of loci that are outliers
	exclude2 <- exclude
	for(i in 1:nrow(exclude)){ # add loci that are <= 5000 bp from outliers
		newlocs <- dat14[CHROM==exclude$CHROM[i] & abs(POS - exclude$POS[i]) <= 5000, .(CHROM, POS)]
		exclude2 <- rbind(exclude2, newlocs)
	}
	allloci <- dat14[outlier==FALSE, .(CHROM, POS)] # list of all loci
	allloci$dist = NA
	allloci$dist[2:nrow(allloci)] <- allloci$POS[2:nrow(allloci)] - allloci$POS[1:(nrow(allloci)-1)]
		dim(allloci)
	allloci <- allloci[dist > 5000 & POS > 5000,] # trim out loci <5000 bp from start of CHROM or from another locus (doesn't check end of CHROM...)
		dim(allloci)
	notoutl <- sample(setdiff(paste(allloci$CHROM, allloci$POS), paste(exclude2$CHROM, exclude2$POS)), nnotoutl) # sample nnotoutl loci to be the "not outliers"
	dat14[paste(CHROM, POS) %in% notoutl, notoutlier:=TRUE] # label the not-outliers
	dat11[paste(CHROM, POS) %in% notoutl, notoutlier:=TRUE] # label the not-outliers
	dat07[paste(CHROM, POS) %in% notoutl, notoutlier:=TRUE] # label the not-outliers


	# for Canada
	exclude <- dat40[outlier==TRUE, .(CHROM, POS)]
	exclude2 <- exclude
	for(i in 1:nrow(exclude)){
		newlocs <- dat40[CHROM==exclude$CHROM[i] & abs(POS - exclude$POS[i]) <= 5000, .(CHROM, POS)]
		exclude2 <- rbind(exclude2, newlocs)
	}
	allloci <- dat40[outlier==FALSE, .(CHROM, POS)] # list of all loci
	allloci$dist = NA
	allloci$dist[2:nrow(allloci)] <- allloci$POS[2:nrow(allloci)] - allloci$POS[1:(nrow(allloci)-1)]
		dim(allloci)
	allloci <- allloci[dist > 5000 & POS > 5000,] # trim out loci <5000 bp from start of CHROM or from another locus (doesn't check end of CHROM...)
		dim(allloci)
	notoutl <- sample(setdiff(paste(allloci$CHROM, allloci$POS), paste(exclude2$CHROM, exclude2$POS)), nnotoutl) # sample nnotoutl loci to be the "not outliers"
	dat40[paste(CHROM, POS) %in% notoutl, notoutlier:=TRUE] # label the not-outliers
	datMod[paste(CHROM, POS) %in% notoutl, notoutlier:=TRUE] # label the not-outliers


	dat14[,summary(notoutlier)] # should == nnotoutl
	dat11[,summary(notoutlier)]
	dat07[,summary(notoutlier)]
	dat40[,summary(notoutlier)]
	datMod[,summary(notoutlier)]

# calculate distance, add monomorphic positions, and set up distance bins from each outlier and not-outlier
# also label which outlier or not-outlier the locus is near (nearPOS column)s
	stp1 <- 100 # step size
	thresh <- 5000 # how far out to go

	if('distclass' %in% colnames(dat14)) dat14[,':=' (distclass=NULL, disttype=NULL, nearPOS=NULL)]
	dat14[, distclass:=as.numeric(NA)]
	dat14[, disttype:=as.factor(NA)]
	dat14[, nearPOS:=as.numeric(NA)]
	levels(dat14$disttype) <- c('outlier', 'notoutlier')
	if('distclass' %in% colnames(dat11)) dat11[,':=' (distclass=NULL, disttype=NULL, nearPOS=NULL)]
	dat11[, distclass:=as.numeric(NA)]
	dat11[, disttype:=as.factor(NA)]
	dat11[, nearPOS:=as.numeric(NA)]
	levels(dat11$disttype) <- c('outlier', 'notoutlier')
	if('distclass' %in% colnames(dat07)) dat07[,':=' (distclass=NULL, disttype=NULL, nearPOS=NULL)]
	dat07[, distclass:=as.numeric(NA)]
	dat07[, disttype:=as.factor(NA)]
	dat07[, nearPOS:=as.numeric(NA)]
	levels(dat07$disttype) <- c('outlier', 'notoutlier')

	if('distclass' %in% colnames(dat40)) dat40[,':=' (distclass=NULL, disttype=NULL, nearPOS=NULL)]
	dat40[, distclass:=as.numeric(NA)]
	dat40[, disttype:=as.factor(NA)]
	dat40[, nearPOS:=as.numeric(NA)]
	levels(dat40$disttype) <- c('outlier', 'notoutlier')
	if('distclass' %in% colnames(datMod)) datMod[,':=' (distclass=NULL, disttype=NULL, nearPOS=NULL)]
	datMod[, distclass:=as.numeric(NA)]
	datMod[, disttype:=as.factor(NA)]
	datMod[, nearPOS:=as.numeric(NA)]
	levels(datMod$disttype) <- c('outlier', 'notoutlier')

	# from each outlier
	# note this will skip any outliers whose 5000 bp overlap any earlier values
		# Norway
	dim(dat14) # 
	dim(dat11)
	dim(dat07)
	for(i in which(dat14$outlier)){
		dists <- dat14[CHROM==dat14$CHROM[i], abs(POS-dat14$POS[i])] # distances from focal SNP. Exact same for dat11 and dat07
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		inds <- dat14$CHROM==dat14$CHROM[i] # indices into dat14 for this chromosome (same for dat11 and dat07)
		inds[inds] <- dists<=thresh # label rows with distance < thresh on this chromosome
		if(any(!is.na(dat14$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat14$distclass[inds])), 'loci already had distance values from another outlier')) # warn if we're close to another locus
		} else {
			dat14[inds, distclass:=distclasses[dists<=thresh]] # add distance class into data.table
			dat14[inds, disttype:='outlier'] # label type of locus
			dat14[inds, nearPOS:=dat14$POS[i]] # label locus it is near
			dat11[inds, distclass:=distclasses[dists<=thresh]] # same for dat11
			dat11[inds, disttype:='outlier']
			dat11[inds, nearPOS:=dat14$POS[i]]
			dat07[inds, distclass:=distclasses[dists<=thresh]] # same for dat07
			dat07[inds, disttype:='outlier']
			dat07[inds, nearPOS:=dat14$POS[i]]
		
			# add monomorphic positions
			mono <- allsites[CHROM==dat14$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat14$POS[dat14$CHROM==dat14$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat14$POS[i])]
			mono <- mono[dists <= thresh,] # trim to nearby SNPs
			mono[, ':=' (PI=0, disttype='outlier', nearPOS=dat14$POS[i], kmer25=as.numeric(NA), dpLofFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat14
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[,dists:=NULL] # remove dists column
		
			dat14 <- rbind(dat14, mono)
			dat11 <- rbind(dat11, mono)
			dat07 <- rbind(dat07, mono)
		}		
	}
	dim(dat14) # 1636722 (union<0.2 no dpfilter no mapfilter) but may vary somewhat run to run? not quite sure why
	dim(dat11)
	dim(dat07)

		# Canada
	dim(dat40) # 
	dim(datMod)
	for(i in which(dat40$outlier)){
		dists <- dat40[CHROM==dat40$CHROM[i], abs(POS-dat40$POS[i])] # distances from focal SNP. Exact same for datMod
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		inds <- dat40$CHROM==dat40$CHROM[i] # indices into dat40 for this chromosome (same for datMod)
		inds[inds] <- dists<=thresh # label rows with distance < 5000 on this chromosome
		if(any(!is.na(dat40$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat40$distclass[inds])), 'loci already had distance values from another outlier')) # warn if we're close to another locus
		} else {
			dat40[inds, distclass:=distclasses[dists<=thresh]] # add distance class into data.table
			dat40[inds, disttype:='outlier'] # label type of locus
			dat40[inds, nearPOS:=dat40$POS[i]] # label locus it is near
			datMod[inds, distclass:=distclasses[dists<=thresh]] # same for datMod
			datMod[inds, disttype:='outlier']
			datMod[inds, nearPOS:=dat40$POS[i]]
		
			# add monomorphic positions
			mono <- allsites[CHROM==dat40$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat40$POS[dat40$CHROM==dat40$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat40$POS[i])]
			mono <- mono[dists <= thresh,]
			mono[, ':=' (PI=0, disttype='outlier', nearPOS=dat40$POS[i], kmer25=as.numeric(NA), dpCanFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat40
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[,dists:=NULL] # remove dists column
		
			dat40 <- rbind(dat40, mono)
			datMod <- rbind(datMod, mono)
		}		
	}
	dim(dat40) # 1117315 (union<0.2 no dpfilter no mapfilter)
	dim(datMod)
	

	# and from each nonoutlier
		# Norway
	dim(dat14) #
	dim(dat11)
	dim(dat07)
	for(i in which(dat14$notoutlier)){
		cat(which(i==which(dat14$notoutlier))) # count to nnotoutl. takes 60 min for 1000 loci?
		
		dists <- dat14[CHROM==dat14$CHROM[i], abs(POS-dat14$POS[i])]
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		inds <- dat14$CHROM==dat14$CHROM[i]
		inds[inds] <- dists<=thresh
		if(any(!is.na(dat14$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat14$distclass[inds])), 'loci already had distance values from another outlier'))
		} else {
			dat14[inds, distclass:=distclasses[dists<=thresh]] # add distance class into data.table
			dat14[inds, disttype:='notoutlier'] # label type of locus
			dat14[inds, nearPOS:=dat14$POS[i]] # label locus it is near
			dat11[inds, distclass:=distclasses[dists<=thresh]] # add distance class into data.table
			dat11[inds, disttype:='notoutlier'] # label type of locus
			dat11[inds, nearPOS:=dat14$POS[i]] # label locus it is near
			dat07[inds, distclass:=distclasses[dists<=thresh]] # add distance class into data.table
			dat07[inds, disttype:='notoutlier'] # label type of locus
			dat07[inds, nearPOS:=dat14$POS[i]] # label locus it is near

			# add monomorphic positions
			mono <- allsites[CHROM==dat14$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat14$POS[dat14$CHROM==dat14$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat14$POS[i])]
			mono <- mono[dists <= thresh,]
			mono[, ':=' (PI=0, disttype='notoutlier', nearPOS=dat14$POS[i], kmer25=as.numeric(NA), dpLofFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat14
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[,dists:=NULL] # remove dists column
		
			dat14 <- rbind(dat14, mono)
			dat11 <- rbind(dat11, mono)
			dat07 <- rbind(dat07, mono)
		}
	}
	dim(dat14) # 2477709 (union<0.2 no dpfilter no mapfilter 100 nnotoutl)
	dim(dat11)
	dim(dat07)

		# Canada
	dim(dat40) # 
	dim(datMod)
	for(i in which(dat40$notoutlier)){ # count to nnotoutl
		cat(which(i==which(dat40$notoutlier)))
		
		dists <- dat40[CHROM==dat40$CHROM[i], abs(POS-dat40$POS[i])]
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		inds <- dat40$CHROM==dat40$CHROM[i]
		inds[inds] <- dists<=thresh
		if(any(!is.na(dat40$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat40$distclass[inds])), 'loci already had distance values from another outlier'))
		} else {
			dat40[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			dat40[inds, disttype:='notoutlier'] # label type of locus
			dat40[inds, nearPOS:=dat40$POS[i]] # label locus it is near
			datMod[inds, distclass:=distclasses[dists<=5000]] # same for datMod
			datMod[inds, disttype:='notoutlier']
			datMod[inds, nearPOS:=dat40$POS[i]] # label locus it is near

			# add monomorphic positions
			mono <- allsites[CHROM==dat40$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat40$POS[dat40$CHROM==dat40$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat40$POS[i])]
			mono <- mono[dists <= thresh,]
			mono[, ':=' (PI=0, disttype='notoutlier', nearPOS=dat40$POS[i], kmer25=as.numeric(NA), dpCanFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat40
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[,dists:=NULL] # remove dists column
		
			dat40 <- rbind(dat40, mono)
			datMod <- rbind(datMod, mono)
		}
	}
	dim(dat40) # 1961704 (union<0.2 no dpfilter no mapfilter)
	dim(datMod)


# calculate averages for each locus
snps14 <- dat14[!is.na(PI) & !is.na(nearPOS),.(nsite=length(PI), nsnp=sum(PI>0)), by=list(CHROM, nearPOS, disttype, distclass)] # average within distance classes for outlier or not outlier (but not neither)

snps11 <- dat11[!is.na(PI) & !is.na(nearPOS),.(nsite=length(PI), nsnp=sum(PI>0)), by=list(CHROM, nearPOS, disttype, distclass)] # average within distance classes for outlier or not outlier (but not neither)

snps07 <- dat07[!is.na(PI) & !is.na(nearPOS),.(nsite=length(PI), nsnp=sum(PI>0)), by=list(CHROM, nearPOS, disttype, distclass)] # average within distance classes for outlier or not outlier (but not neither)

snps40 <- dat40[!is.na(PI) & !is.na(nearPOS),.(nsite=length(PI), nsnp=sum(PI>0)), by=list(CHROM, nearPOS, disttype, distclass)] # average within distance classes for outlier or not outlier (but not neither)

snpsMod <- datMod[!is.na(PI) & !is.na(nearPOS),.(nsite=length(PI), nsnp=sum(PI>0)), by=list(CHROM, nearPOS, disttype, distclass)] # average within distance classes for outlier or not outlier (but not neither)


setkey(snps14, disttype, CHROM, nearPOS, distclass)
setkey(snps11, disttype, CHROM, nearPOS, distclass)
setkey(snps07, disttype, CHROM, nearPOS, distclass)
setkey(snps40, disttype, CHROM, nearPOS, distclass)
setkey(snpsMod, disttype, CHROM, nearPOS, distclass)


# examine
snps14
snps07

snps14[distclass<400,.(snpdens.mean=mean(nsnp/nsite), snpdens.sd=sd(nsnp/nsite)), by=c('disttype', 'distclass')]
snps11[distclass<400,.(snpdens.mean=mean(nsnp/nsite), snpdens.sd=sd(nsnp/nsite)), by=c('disttype', 'distclass')]
snps07[distclass<400,.(snpdens.mean=mean(nsnp/nsite), snpdens.sd=sd(nsnp/nsite)), by=c('disttype', 'distclass')]

snpsMod[distclass<400,.(snpdens.mean=mean(nsnp/nsite), snpdens.sd=sd(nsnp/nsite)), by=c('disttype', 'distclass')]
snps40[distclass<400,.(snpdens.mean=mean(nsnp/nsite), snpdens.sd=sd(nsnp/nsite)), by=c('disttype', 'distclass')]


# label by population
snps14[,pop:='LOF_S_14']
snps11[,pop:='LOF_S_11']
snps07[,pop:='LOF_07']
snps40[,pop:='CAN40']
snpsMod[,pop:='CANMod']

# merge
snps <- rbind(snps07, snps11, snps14, snps40, snpsMod)

# write out
if(dpfilter==TRUE & mapfilter==TRUE) outfile <- paste('analysis/snpdens_outliers_', outliertype, '.csv', sep='')
if(dpfilter==FALSE & mapfilter==TRUE) outfile <- paste('analysis/snpdens_outliers_', outliertype, '_nodpfilter.csv', sep='')
if(dpfilter==FALSE & mapfilter==FALSE) outfile <- paste('analysis/snpdens_outliers_', outliertype, '_nodpfilter_nomapfilter.csv', sep='')
outfile
write.csv(snps, file=outfile, row.names=FALSE)


#####################
## examine
# to run on macbook
####################
snps <- fread('analysis/snpdens_outliers_union_nodpfilter_nomapfilter.csv', ); outliertype <- 'union_nodpfilter_nomapfilter'
setkey(snps, disttype, CHROM, nearPOS, distclass)


# which loci have elevated snps in modern?
snps[distclass==50 & nsnp/nsite>0.04, .(CHROM, nearPOS, disttype, nsite, nsnp, dens=round(nsnp/nsite,3), pop)]


###############
# plot
# to run on macbook
###############
require(RColorBrewer)
require(data.table)


snps <- fread('analysis/snpdens_outliers_union_nodpfilter_nomapfilter.csv', ); outliertype <- 'union_nodpfilter_nomapfilter'

setkey(snps, disttype, CHROM, nearPOS, distclass)


# plot lines on top of each other
cols3 <- c('#00000055', '#00000055', '#00000055', '#00000055', '#00000055') # Lof07, Lof11, Lof14, Can40, CanMod (grey)
xmx <- 500

quartz(width=10, height=6)
# pdf(width=8, height=4, file=paste('figures/snp_density_outliers_bylocus_', outliertype, '.pdf', sep=''))
par(mfcol=c(3,4), mai=c(0.6, 0.6, 0.2, 0.1))
plot(0,0, type='n', xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='LOF_07 outlier')
snps[disttype=='outlier' & pop=='LOF_07', lines(distclass, nsnp/nsite, col=cols3[1]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='LOF_S_11 outlier')
snps[disttype=='outlier' & pop=='LOF_S_11', lines(distclass, nsnp/nsite, col=cols3[2]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='LOF_S_14 outlier')
snps[disttype=='outlier' & pop=='LOF_S_14', lines(distclass, nsnp/nsite, col=cols3[3]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n', xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='LOF_07 not outlier')
snps[disttype=='notoutlier' & pop=='LOF_07', lines(distclass, nsnp/nsite, col=cols3[1]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='LOF_S_11 not outlier')
snps[disttype=='notoutlier' & pop=='LOF_S_11', lines(distclass, nsnp/nsite, col=cols3[2]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='LOF_S_14 not outlier')
snps[disttype=='notoutlier' & pop=='LOF_S_14', lines(distclass, nsnp/nsite, col=cols3[3]), by=c('CHROM', 'nearPOS')]


plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='CAN40 outlier')
snps[disttype=='outlier' & pop=='CAN40', lines(distclass, nsnp/nsite, col=cols3[1]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='CANMod outlier')
snps[disttype=='outlier' & pop=='CANMod', lines(distclass, nsnp/nsite, col=cols3[2]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n',bty='n', xaxt='n', yaxt='n', xlab='', ylab='', main='')

plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='CAN40 not outlier')
snps[disttype=='notoutlier' & pop=='CAN40', lines(distclass, nsnp/nsite, col=cols3[1]), by=c('CHROM', 'nearPOS')]

plot(0,0, type='n',xlim=c(0,xmx), ylim=c(0, snps[,max(nsnp/nsite)]), xlab='Distance (bp)', ylab='SNP density (# snps/#sites)', main='CANMod not outlier')
snps[disttype=='notoutlier' & pop=='CANMod', lines(distclass, nsnp/nsite, col=cols3[2]), by=c('CHROM', 'nearPOS')]

dev.off()

# plot mean snp density and #snps
cols3 <- c('#e7d4e8', '#af8dc3', '#762a83', '#d9f0d3', '#1b7837') # Lof07, Lof11, Lof14, Can40, CanMod (PRGn colorbrewer)
bins <- snps[, .(nsnp = mean(nsnp), snpdens = mean(nsnp/nsite), nsnpsd = sd(nsnp), snpdenssd = sd(nsnp/nsite), n=.N), by=c('pop', 'disttype', 'distclass')]
bins[,':=' (nsnpse=nsnpsd/n, snpdensse=snpdenssd/n)]

xmx <- 500
lwd <- 2

quartz(width=8, height=4)
# pdf(width=8, height=4, file=paste('figures/snp_density_outliers_', outliertype, '.pdf', sep=''))
par(mfrow=c(2,2), mai=c(0.6, 1, 0.2, 0.1), mgp=c(2,0.7,0), tcl=-0.3)
plot(0,0, type='n', xlim=c(0,xmx), ylim=c(0, bins[,max(snpdens)]), xlab='Distance (bp)', ylab='SNP density', main='')
bins[disttype=='outlier' & pop=='LOF_07', lines(distclass, snpdens, col=cols3[1], lwd=lwd)]
bins[disttype=='outlier' & pop=='LOF_S_11', lines(distclass, snpdens, col=cols3[2], lwd=lwd)]
bins[disttype=='outlier' & pop=='LOF_S_14', lines(distclass, snpdens, col=cols3[3], lwd=lwd)]

	bins[disttype=='outlier' & pop=='LOF_07', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[1]), by=distclass]
	bins[disttype=='outlier' & pop=='LOF_S_11', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[2]), by=distclass]
	bins[disttype=='outlier' & pop=='LOF_S_14', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[3]), by=distclass]

bins[disttype=='notoutlier' & pop=='LOF_07', lines(distclass, snpdens, col=cols3[1], lwd=lwd, lty=2)]
bins[disttype=='notoutlier' & pop=='LOF_S_11', lines(distclass, snpdens, col=cols3[2], lwd=lwd, lty=2)]
bins[disttype=='notoutlier' & pop=='LOF_S_14', lines(distclass, snpdens, col=cols3[3], lwd=lwd, lty=2)]

	bins[disttype=='notoutlier' & pop=='LOF_07', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[1]), by=distclass]
	bins[disttype=='notoutlier' & pop=='LOF_S_11', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[2]), by=distclass]
	bins[disttype=='notoutlier' & pop=='LOF_S_14', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[3]), by=distclass]

legend('topright', legend=c('Norway 1907', 'Norway 2011', 'Norway 2014', 'Norway outlier', 'Not outlier'), col=c(cols3[1:3], 'black', 'black'), lty=c(rep(1,3), 1, 3), bty='n', cex=0.8)


plot(0,0, type='n', xlim=c(0,xmx), ylim=c(0, bins[,max(snpdens)]), xlab='Distance (bp)', ylab='SNP density', main='')
bins[disttype=='outlier' & pop=='CAN40', lines(distclass, snpdens, col=cols3[4], lwd=lwd)]
bins[disttype=='outlier' & pop=='CANMod', lines(distclass, snpdens, col=cols3[5], lwd=lwd)]

	bins[disttype=='outlier' & pop=='CAN40', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[4]), by=distclass]
	bins[disttype=='outlier' & pop=='CANMod', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[5]), by=distclass]

bins[disttype=='notoutlier' & pop=='CAN40', lines(distclass, snpdens, col=cols3[4], lwd=lwd, lty=2)]
bins[disttype=='notoutlier' & pop=='CANMod', lines(distclass, snpdens, col=cols3[5], lwd=lwd, lty=2)]

	bins[disttype=='notoutlier' & pop=='CAN40', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[4]), by=distclass]
	bins[disttype=='notoutlier' & pop=='CANMod', lines(c(distclass, distclass), c(snpdens-snpdensse, snpdens+snpdensse), col=cols3[5]), by=distclass]

legend('topright', legend=c('Canada 1940', 'Canada 2013', 'Canada outlier', 'Not outlier'), col=c(cols3[4:5], 'black', 'black'), lty=c(rep(1,2), 1, 3), bty='n', cex=0.8)


plot(0,0, type='n', xlim=c(0,xmx), ylim=c(0, bins[,max(nsnp)]), xlab='Distance (bp)', ylab='# SNPs', main='')
bins[disttype=='outlier' & pop=='LOF_07', lines(distclass, nsnp, col=cols3[1], lwd=lwd)]
bins[disttype=='outlier' & pop=='LOF_S_11', lines(distclass, nsnp, col=cols3[2], lwd=lwd)]
bins[disttype=='outlier' & pop=='LOF_S_14', lines(distclass, nsnp, col=cols3[3], lwd=lwd)]

	bins[disttype=='outlier' & pop=='LOF_07', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[1]), by=distclass]
	bins[disttype=='outlier' & pop=='LOF_S_11', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[2]), by=distclass]
	bins[disttype=='outlier' & pop=='LOF_S_14', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[3]), by=distclass]

bins[disttype=='notoutlier' & pop=='LOF_07', lines(distclass, nsnp, col=cols3[1], lwd=lwd, lty=2)]
bins[disttype=='notoutlier' & pop=='LOF_S_11', lines(distclass, nsnp, col=cols3[2], lwd=lwd, lty=2)]
bins[disttype=='notoutlier' & pop=='LOF_S_14', lines(distclass, nsnp, col=cols3[3], lwd=lwd, lty=2)]

	bins[disttype=='notoutlier' & pop=='LOF_07', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[1]), by=distclass]
	bins[disttype=='notoutlier' & pop=='LOF_S_11', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[2]), by=distclass]
	bins[disttype=='notoutlier' & pop=='LOF_S_14', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[3]), by=distclass]

plot(0,0, type='n', xlim=c(0,xmx), ylim=c(0, bins[,max(nsnp)]), xlab='Distance (bp)', ylab='# SNPs', main='')
bins[disttype=='outlier' & pop=='CAN40', lines(distclass, nsnp, col=cols3[4], lwd=lwd)]
bins[disttype=='outlier' & pop=='CANMod', lines(distclass, nsnp, col=cols3[5], lwd=lwd)]

	bins[disttype=='outlier' & pop=='CAN40', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[4]), by=distclass]
	bins[disttype=='outlier' & pop=='CANMod', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[5]), by=distclass]

bins[disttype=='notoutlier' & pop=='CAN40', lines(distclass, nsnp, col=cols3[4], lwd=lwd, lty=2)]
bins[disttype=='notoutlier' & pop=='CANMod', lines(distclass, nsnp, col=cols3[5], lwd=lwd, lty=2)]

	bins[disttype=='notoutlier' & pop=='CAN40', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[4]), by=distclass]
	bins[disttype=='notoutlier' & pop=='CANMod', lines(c(distclass, distclass), c(nsnp-nsnpse, nsnp+nsnpse), col=cols3[5]), by=distclass]

dev.off()