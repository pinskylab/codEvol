# load functions
source('scripts/wfs.r')
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel)
	require(data.table)
	require(abc)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
	require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# load data
load('analysis/wfs_sims.rdata') # load the simulations: 'out' array. slow to load...
obs <- fread('analysis/LOF_07_to_LOF_S_14.w_obs_stats.txt') # the observations of Fsi and Fsd
dat <- fread('analysis/LOF_07_to_LOF_S_14.wfabc', skip=2) # the data on samples sizes and observed allele frequencies

# prep data
sampsize <- dat[seq(1,nrow(dat),by=2),] # sample sizes in # of chromosomes
alcnts <- dat[seq(2,nrow(dat),by=2),] # allele counts in # chromosomes


# calculate mean and sd for each summary stat in out (simulations)
# to allow scaling of target stats and simulation stats
	uniqsamps <- apply(do.call('rbind', strsplit(dimnames(out)[[3]], split=',')), 2, as.numeric) # extract unique sample sizes. could also have extracted (probably more easily) from targ
	meansds <- data.table(alcnt1=uniqsamps[,1], alcnt2=uniqsamps[,2])
	meansds[,f1samp.mean:=numeric(nrow(meansds))]
	meansds[,f1samp.sd:=numeric(nrow(meansds))]
	meansds[,fsd.mean:=numeric(nrow(meansds))]
	meansds[,fsd.sd:=numeric(nrow(meansds))]
	meansds[,fsi.mean:=numeric(nrow(meansds))]
	meansds[,fsi.sd:=numeric(nrow(meansds))]
	for(i in 1:nrow(meansds)){ # this takes 22%+ of a cod node RAM... why? perhaps because mean and sd on vectors of length 10M
		print(i)
		j <- which(dimnames(out)[[3]] == paste(meansds[i,.(alcnt1, alcnt2)], collapse=',')) # have to find the slice of out that matches the sample sizes for this row
		meansds[i,f1samp.mean:=mean(out['f1samp',,j], na.rm=TRUE)]
		meansds[i,f1samp.sd:=sd(out['f1samp',,j], na.rm=TRUE)]
		meansds[i,fsd.mean:=mean(out['fsdprime',,j], na.rm=TRUE)]
		meansds[i,fsd.sd:=sd(out['fsdprime',,j], na.rm=TRUE)]
		meansds[i,fsi.mean:=mean(out['fsiprime',,j], na.rm=TRUE)]
		meansds[i,fsi.sd:=sd(out['fsiprime',,j], na.rm=TRUE)]
	}	
	meansds[f1samp.sd==0, f1samp.sd:=1] # set sd = 1 for cases were sd=0, to prevent Inf
	meansds[fsd.sd==0, fsd.sd:=1] # set sd = 1 for cases were sd=0, to prevent Inf
	meansds[fsi.sd==0, fsi.sd:=1] # set sd = 1 for cases were sd=0, to prevent Inf

	# write out
	save(meansds, file='analysis/wfs_sims_meansds.rdata')

# set up target statistics
	load('analysis/wfs_sims_meansds.rdata')
	# stats for each locus
	targ <- cbind(1:nrow(sampsize), sampsize[,.(V1, V2)], alcnts[,V1]/sampsize[,V1], obs[,.(Fsd, Fsi)]) # alcnt1, alcnt2, f1samp, fsd, fsi
	setnames(targ, 1:6, c('locusnum', 'alcnt1', 'alcnt2', 'f1samp', 'fsd', 'fsi'))

	# add mean and sd to targ to allow easy normalization
	setkey(targ, alcnt1, alcnt2)
	setkey(meansds, alcnt1, alcnt2)
	targ <- meansds[targ,]

	# normalize the target stats
	targ[,f1samp.n:=(f1samp-f1samp.mean)/f1samp.sd]
	targ[,fsd.n:=(fsd-fsd.mean)/fsd.sd]
	targ[,fsi.n:=(fsi-fsi.mean)/fsi.sd]

	# re-order rows and columns
	setkey(targ, locusnum)
	setcolorder(targ, c('locusnum', 'alcnt1', 'alcnt2', 'f1samp.n', 'fsd.n', 'fsi.n', 'f1samp', 'fsd', 'fsi', 'f1samp.mean', 'fsd.mean', 'fsi.mean', 'f1samp.sd', 'fsd.sd', 'fsi.sd'))

	# write out
	save(targ, file='analysis/wfs_targ.rdata')	


# export out file from abc sims to ff data type for memory-efficient use across cluster
# also scale simulation statistics to mean 0 sd 1
	load('analysis/wfs_sims_meansds.rdata')
	options('fftempdir') # check where this is to monitor temp file creation in the file system /tmp/Rtmpp__

	# data.table to hold means, sds, and number of na values in out
	checkscl <- meansds[,.(alcnt1, alcnt2)]
	checkscl[,f1samp.mean:=numeric(nrow(checkscl))]
	checkscl[,f1samp.sd:=numeric(nrow(checkscl))]
	checkscl[,f1samp.na:=numeric(nrow(checkscl))]
	checkscl[,fsd.mean:=numeric(nrow(checkscl))]
	checkscl[,fsd.sd:=numeric(nrow(checkscl))]
	checkscl[,fsd.na:=numeric(nrow(checkscl))]
	checkscl[,fsi.mean:=numeric(nrow(checkscl))]
	checkscl[,fsi.sd:=numeric(nrow(checkscl))]
	checkscl[,fsi.na:=numeric(nrow(checkscl))]

	mean2 <- function(x) sum(x, na.rm=TRUE)/length(x) # faster version of mean

	for(i in 1:nrow(meansds)){
		print(paste('iteration', i))
		j <- which(dimnames(out)[[3]] == paste(meansds[i,.(alcnt1, alcnt2)], collapse=',')) # have to find the slice of out that matches the sample sizes for this row of meansds
		thisout <- out[,,j]

		# scale everything
		thisout['f1samp',] <- (thisout['f1samp',]-meansds[i,f1samp.mean])/meansds[i,f1samp.sd]
		thisout['fsdprime',] <- (thisout['fsdprime',]-meansds[i,fsd.mean])/meansds[i,fsd.sd]
		thisout['fsiprime',] <- (thisout['fsiprime',]-meansds[i,fsi.mean])/meansds[i,fsi.sd]
		
		# check scaling
		checkscl[i,f1samp.mean := mean2(thisout['f1samp',])]
		checkscl[i,f1samp.sd := sd(thisout['f1samp',], na.rm=TRUE)]
		checkscl[i,fsd.mean := mean2(thisout['fsdprime',])]
		checkscl[i,fsd.sd := sd(thisout['fsdprime',], na.rm=TRUE)]
		checkscl[i,fsi.mean := mean2(thisout['fsiprime',])]
		checkscl[i,fsi.sd := sd(thisout['fsiprime',], na.rm=TRUE)]
		checkscl[i,f1samp.na := sum(is.na(thisout['f1samp',]))]
		checkscl[i,fsd.na := sum(is.na(thisout['fsdprime',]))]
		checkscl[i,fsi.na := sum(is.na(thisout['fsiprime',]))]

		print(checkscl[i,])

		# make ff object
		thisout.ff <- ff(thisout, dim=dim(thisout), dimnames=dimnames(thisout)) # create in tempdir
		
		# save to permanent file (semi-permanent.. but in my temp directory)
		ffsave(thisout.ff, file=paste('analysis/temp/wfs_sims_ff', paste(meansds[i,.(alcnt1, alcnt2)], collapse=','), sep=''))

		# remove
		delete(thisout.ff)
		rm(thisout.ff)
	}
	
	save(checkscl, file='analysis/temp/wfs_sims_checkscl.rdata')

	# look for suspect normalizations
	checkscl[f1samp.mean > 1e-12 | f1samp.mean < -1e-12 | f1samp.sd > 1.0001 | f1samp.sd < 0.9999,]
	checkscl[fsd.mean > 1e-12 | fsd.mean < -1e-12 | fsd.sd > 1.0001 | fsd.sd < 0.9999,]
	checkscl[fsi.mean > 1e-12 | fsi.mean < -1e-12 | fsi.sd > 1.0001 | fsi.sd < 0.9999,]
	