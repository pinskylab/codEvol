# compare pi near outlier loci and near non-outlier loci
# account for which loci can be called as SNPs/not SNPs


######################
# calculate pi in windows
# to run on a cod node!
######################
require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
outl <- fread("zcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz") # 270MB file
#	outl

# read in pi data on cod node (note: includes inversions)
dat14 <- fread('analysis/LOF_S_14.sites.pi')
dat11 <- fread('analysis/LOF_S_11.sites.pi')
dat07 <- fread('analysis/LOF_07.sites.pi')
dat40 <- fread('analysis/CAN40.sites.pi')
datMod <- fread('analysis/CANMod.sites.pi')

# read in all callable sites information
allsites <- fread('zcat all_sites_data_29_06_18/AllSites.kept.sites.kmer25.gz') # gzipped file
setkey(allsites, CHROM, POS)

# merge in outlier information to dat
# use q<0.3 for more loci in pi calculations
nrow(dat14) # 1402283
nrow(dat11)
nrow(dat07)
nrow(dat40) # 1018417
nrow(datMod)

dat14 <- merge(dat14, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=p.Lof.adj3<0.3)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)

dat11 <- merge(dat11, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=p.Lof.adj3<0.3)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)
	
dat07 <- merge(dat07, outl[,.(CHROM, POS, kmer25, dpLofFlag, outlier=p.Lof.adj3<0.3)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)

dat40 <- merge(dat40, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=p.Can.adj3<0.3)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)

datMod <- merge(datMod, outl[,.(CHROM, POS, kmer25, dpCanFlag, outlier=p.Can.adj3<0.3)], by.x=c('CHROM', 'POS'), by.y=c('CHROM', 'POS'), all.x=TRUE)

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

dat14 <- dat14[kmer25==TRUE & dpLofFlag==TRUE,]
dat11 <- dat11[kmer25==TRUE & dpLofFlag==TRUE,]
dat07 <- dat07[kmer25==TRUE & dpLofFlag==TRUE,]
dat40 <- dat40[kmer25==TRUE & dpCanFlag==TRUE,]
datMod <- datMod[kmer25==TRUE & dpCanFlag==TRUE,]

nrow(dat14) # 408771
nrow(dat11)
nrow(dat07)
nrow(dat40) # 313787
nrow(datMod)

# how many outlier comparisons?
# WHY SO MANY NAs IN THE OUTLIER COLUMN? APPARENTLY WEREN'T ANY FOR THE LD CALCULATIONS
dat14[,summary(outlier)]
dat11[,summary(outlier)]
dat07[,summary(outlier)]
dat40[,summary(outlier)]
datMod[,summary(outlier)]

# set up 1000 not-outlier loci
	dat14[,notoutlier:=FALSE]
	dat11[,notoutlier:=FALSE]
	dat07[,notoutlier:=FALSE]
	dat40[,notoutlier:=FALSE]
	datMod[,notoutlier:=FALSE]

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
	notoutl <- sample(setdiff(paste(allloci$CHROM, allloci$POS), paste(exclude2$CHROM, exclude2$POS)), 1000) # sample 1000 loci to be the "not outliers"
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
	notoutl <- sample(setdiff(paste(allloci$CHROM, allloci$POS), paste(exclude2$CHROM, exclude2$POS)), 1000) # sample 1000 loci to be the "not outliers"
	dat40[paste(CHROM, POS) %in% notoutl, notoutlier:=TRUE] # label the not-outliers
	datMod[paste(CHROM, POS) %in% notoutl, notoutlier:=TRUE] # label the not-outliers


	dat14[,summary(notoutlier)]
	dat11[,summary(notoutlier)]
	dat07[,summary(notoutlier)]
	dat40[,summary(notoutlier)]
	datMod[,summary(notoutlier)]

# calculate distance, add monomorphic positions, and set up distance bins from each outlier
	stp1 <- 10 # step size for small distances
	thresh <- 50 # use stp2 above this distance
	stp2 <- 1000

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
	# note this will skip any outliers whose 5000 bp overlap any earlier values
		# Norway
	dim(dat14) # 408771
	dim(dat11)
	dim(dat07)
	for(i in which(dat14$outlier)){
		dists <- dat14[CHROM==dat14$CHROM[i], abs(POS-dat14$POS[i])] # distances from focal SNP. Exact same for dat11 and dat07
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		distclasses[dists>thresh] <- floor(dists[dists>thresh]/stp2)*stp2+stp2/2 # second step size
		inds <- dat14$CHROM==dat14$CHROM[i] # indices into dat14 for this chromosome (same for dat11 and dat07)
		inds[inds] <- dists<=5000 # label rows with distance < 5000 on this chromosome
		if(any(!is.na(dat14$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat14$distclass[inds])), 'loci already had distance values from another outlier')) # warn if we're close to another locus
		} else {
			dat14[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			dat14[inds, disttype:='outlier'] # label type of locus
			dat11[inds, distclass:=distclasses[dists<=5000]] # same for dat11
			dat11[inds, disttype:='outlier']
			dat07[inds, distclass:=distclasses[dists<=5000]] # same for dat07
			dat07[inds, disttype:='outlier']
		
			# add monomorphic positions
			mono <- allsites[CHROM==dat14$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat14$POS[dat14$CHROM==dat14$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat14$POS[i])]
			mono <- mono[dists <= 5000,]
			mono[, ':=' (PI=0, disttype='outlier', kmer25=as.numeric(NA), dpLofFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat14
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[dists>thresh, distclass := floor(dists[dists>thresh]/stp2)*stp2+stp2/2]
			mono[,dists:=NULL] # remove dists column
		
			dat14 <- rbind(dat14, mono)
			dat11 <- rbind(dat11, mono)
			dat07 <- rbind(dat07, mono)
		}		
	}
	dim(dat14) #470947
	dim(dat11)
	dim(dat07)

		# Canada
	dim(dat40) # 313787
	dim(datMod)
	for(i in which(dat40$outlier)){
		dists <- dat40[CHROM==dat40$CHROM[i], abs(POS-dat40$POS[i])] # distances from focal SNP. Exact same for datMod
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		distclasses[dists>thresh] <- floor(dists[dists>thresh]/stp2)*stp2+stp2/2 # second step size
		inds <- dat40$CHROM==dat40$CHROM[i] # indices into dat40 for this chromosome (same for datMod)
		inds[inds] <- dists<=5000 # label rows with distance < 5000 on this chromosome
		if(any(!is.na(dat40$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat40$distclass[inds])), 'loci already had distance values from another outlier')) # warn if we're close to another locus
		} else {
			dat40[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			dat40[inds, disttype:='outlier'] # label type of locus
			datMod[inds, distclass:=distclasses[dists<=5000]] # same for datMod
			datMod[inds, disttype:='outlier']
		
			# add monomorphic positions
			mono <- allsites[CHROM==dat40$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat40$POS[dat40$CHROM==dat40$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat40$POS[i])]
			mono <- mono[dists <= 5000,]
			mono[, ':=' (PI=0, disttype='outlier', kmer25=as.numeric(NA), dpCanFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat40
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[dists>thresh, distclass := floor(dists[dists>thresh]/stp2)*stp2+stp2/2]
			mono[,dists:=NULL] # remove dists column
		
			dat40 <- rbind(dat40, mono)
			datMod <- rbind(datMod, mono)
		}		
	}
	dim(dat40) # 317928
	dim(datMod)
	

	# and from each nonoutlier
		# Norway
	dim(dat14) #470947
	dim(dat11)
	dim(dat07)
	for(i in which(dat14$notoutlier)){
		cat(which(i==which(dat14$notoutlier)))
		
		dists <- dat14[CHROM==dat14$CHROM[i], abs(POS-dat14$POS[i])]
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		distclasses[dists>thresh] <- floor(dists[dists>thresh]/stp2)*stp2+stp2/2 # second step size
		inds <- dat14$CHROM==dat14$CHROM[i]
		inds[inds] <- dists<=5000
		if(any(!is.na(dat14$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat14$distclass[inds])), 'loci already had distance values from another outlier'))
		} else {
			dat14[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			dat14[inds, disttype:='notoutlier'] # label type of locus
			dat11[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			dat11[inds, disttype:='notoutlier'] # label type of locus
			dat07[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			dat07[inds, disttype:='notoutlier'] # label type of locus

			# add monomorphic positions
			mono <- allsites[CHROM==dat14$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat14$POS[dat14$CHROM==dat14$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat14$POS[i])]
			mono <- mono[dists <= 5000,]
			mono[, ':=' (PI=0, disttype='notoutlier', kmer25=as.numeric(NA), dpLofFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat14
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[dists>thresh, distclass := floor(dists[dists>thresh]/stp2)*stp2+stp2/2]
			mono[,dists:=NULL] # remove dists column
		
			dat14 <- rbind(dat14, mono)
			dat11 <- rbind(dat11, mono)
			dat07 <- rbind(dat07, mono)
		}
	}
	dim(dat14) # 2284776
	dim(dat11)
	dim(dat07)

		# Canada
	dim(dat40) # 317928
	dim(datMod)
	for(i in which(dat40$notoutlier)){
		cat(which(i==which(dat40$notoutlier)))
		
		dists <- dat40[CHROM==dat40$CHROM[i], abs(POS-dat40$POS[i])]
		distclasses <- floor(dists/stp1)*stp1+stp1/2 # calculate a distance class
		distclasses[dists>thresh] <- floor(dists[dists>thresh]/stp2)*stp2+stp2/2 # second step size
		inds <- dat40$CHROM==dat40$CHROM[i]
		inds[inds] <- dists<=5000
		if(any(!is.na(dat40$distclass[inds]))){
			print(paste('i=', i, ',', sum(!is.na(dat40$distclass[inds])), 'loci already had distance values from another outlier'))
		} else {
			dat40[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			dat40[inds, disttype:='notoutlier'] # label type of locus
			datMod[inds, distclass:=distclasses[dists<=5000]] # add distance class into data.table
			datMod[inds, disttype:='notoutlier'] # label type of locus

			# add monomorphic positions
			mono <- allsites[CHROM==dat40$CHROM[i], .(CHROM, POS)] # trim to all loci on this chromosome
			mono <- mono[!(POS %in% dat40$POS[dat40$CHROM==dat40$CHROM[i]]),] # trim out position on this chromosome already in the dataframe (SNPs or monomorphic positions added for nearby SNPs)
			mono[, dists := abs(POS-dat40$POS[i])]
			mono <- mono[dists <= 5000,]
			mono[, ':=' (PI=0, disttype='notoutlier', kmer25=as.numeric(NA), dpCanFlag=as.logical(NA), outlier=as.logical(NA), notoutlier=as.logical(NA))] # add PI=0 and other columns for binding to dat40
			mono[, distclass := floor(dists/stp1)*stp1+stp1/2]
			mono[dists>thresh, distclass := floor(dists[dists>thresh]/stp2)*stp2+stp2/2]
			mono[,dists:=NULL] # remove dists column
		
			dat40 <- rbind(dat40, mono)
			datMod <- rbind(datMod, mono)
		}
	}
	dim(dat40) # 2154018
	dim(datMod)



# calculate averages within bins (like Bario et al. 2016 eLife)
bins14 <- dat14[!is.na(PI),.(piave=mean(PI), pisd=sd(PI), pil95=quantile(PI,probs=0.025), piu95=quantile(PI,probs=0.975), pil75=quantile(PI,probs=0.125), piu75=quantile(PI,probs=0.875), pin=.N), by=list(disttype, distclass)] # average within distance classes for outlier or not outlier (but not neither)

bins11 <- dat11[!is.na(PI),.(piave=mean(PI), pisd=sd(PI), pil95=quantile(PI,probs=0.025), piu95=quantile(PI,probs=0.975), pil75=quantile(PI,probs=0.125), piu75=quantile(PI,probs=0.875), pin=.N), by=list(disttype, distclass)]

bins07 <- dat07[!is.na(PI),.(piave=mean(PI), pisd=sd(PI), pil95=quantile(PI,probs=0.025), piu95=quantile(PI,probs=0.975), pil75=quantile(PI,probs=0.125), piu75=quantile(PI,probs=0.875), pin=.N), by=list(disttype, distclass)]

bins40 <- dat40[!is.na(PI),.(piave=mean(PI), pisd=sd(PI), pil95=quantile(PI,probs=0.025), piu95=quantile(PI,probs=0.975), pil75=quantile(PI,probs=0.125), piu75=quantile(PI,probs=0.875), pin=.N), by=list(disttype, distclass)]

binsMod <- datMod[!is.na(PI),.(piave=mean(PI), pisd=sd(PI), pil95=quantile(PI,probs=0.025), piu95=quantile(PI,probs=0.975), pil75=quantile(PI,probs=0.125), piu75=quantile(PI,probs=0.875), pin=.N), by=list(disttype, distclass)]



bins14[,pise:=pisd/pin]
bins11[,pise:=pisd/pin]
bins07[,pise:=pisd/pin]
bins40[,pise:=pisd/pin]
binsMod[,pise:=pisd/pin]

setkey(bins14, disttype, distclass)
setkey(bins11, disttype, distclass)
setkey(bins07, disttype, distclass)
setkey(bins40, disttype, distclass)
setkey(binsMod, disttype, distclass)

# examine
bins14[,.(disttype,distclass,piave, pisd, pin, pise)]
bins07[,.(disttype,distclass,piave, pisd, pin, pise)]

# label by population
bins14[,pop:='LOF_S_14']
bins11[,pop:='LOF_S_11']
bins07[,pop:='LOF_07']
bins40[,pop:='CAN40']
binsMod[,pop:='CANMod']

# merge
bins <- rbind(bins07, bins11, bins14, bins40, binsMod)

# write out
write.csv(bins, file='analysis/pi_outliers.csv', row.names=FALSE)

###############
# plot
# to run on macbook
###############
require(RColorBrewer)
require(data.table)

bins <- fread('analysis/pi_outliers.csv', )
cols <- brewer.pal(5, 'Set1')
cex=0.5
lwd=2
ylims <- bins[!is.na(disttype),range(piave)]

# plot binned data
quartz(width=6, height=3)
# pdf(width=6, height=3, file='figures/pi_decay_outliers.pdf')
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Lof
bins[pop=='LOF_S_14' & disttype=='outlier',plot(distclass, piave, ylim=ylims, type='l', xlab='Distance (bp)', ylab='Pi per site', cex=cex, main='Lof', log='x', col=cols[1], lwd=lwd)]
bins[pop=='LOF_S_11' & disttype=='outlier',lines(distclass, piave, type='l', cex=cex, col=cols[2], lwd=lwd)]
bins[pop=='LOF_07' & disttype=='outlier',lines(distclass, piave, type='l', cex=cex, col=cols[3], lwd=lwd)]

bins[pop=='LOF_S_14' & disttype=='nooutlier',lines(distclass, piave, type='l', cex=cex, col=cols[1], lty=3, lwd=lwd)]
bins[pop=='LOF_S_11' & disttype=='nooutlier',lines(distclass, piave, type='l', cex=cex, col=cols[2], lty=3, lwd=lwd)]
bins[pop=='LOF_07' & disttype=='notoutlier',lines(distclass, piave, type='l', cex=cex, col=cols[3], lty=3, lwd=lwd)]

bins[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[1]), by=distclass] # 95% CI on the mean.
bins[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[3]), by=distclass]

bins[pop=='LOF_S_14' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[1]), by=distclass]
bins[pop=='LOF_S_11' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[2]), by=distclass]
bins[pop=='LOF_07' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[3]), by=distclass]

legend('topright', legend=c('LOF_S_14', 'LOF_S_11', 'LOF_07', 'Outlier', 'Not outlier'), col=c(cols[1:3], 'black', 'black'), lty=c(rep(1,4), 3), bty='n', cex=0.4)

	# Canada
bins[pop=='CANMod' & disttype=='outlier',plot(distclass, piave, ylim=c(0,1), type='l', xlab='Distance (bp)', ylab='Pi per site', cex=cex, main='Can', log='x', col=cols[4], lwd=lwd)]
bins[pop=='CAN40' & disttype=='outlier',lines(distclass, piave, type='l', cex=cex, col=cols[5], lwd=lwd)]

bins[pop=='CANMod' & disttype=='notoutlier',lines(distclass, piave, type='l', cex=cex, col=cols[4], lty=3, lwd=lwd)]
bins[pop=='CAN40' & disttype=='notoutlier',lines(distclass, piave, type='l', cex=cex, col=cols[5], lty=3, lwd=lwd)]

bins[pop=='CANMod' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[4]), by=distclass]
bins[pop=='CAN40' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[5]), by=distclass]
bins[pop=='CANMod' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[4]), by=distclass]
bins[pop=='CAN40' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols[5]), by=distclass]

legend('topright', legend=c('CANMod', 'CAN40', 'Outlier', 'Not outlier'), col=c(cols[4:5], 'black', 'black'), lty=c(1,1,1,3), bty='n', cex=0.4)

dev.off()