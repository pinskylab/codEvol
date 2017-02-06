# ABC on each locus from ABC lookup table to estimate s
# run after wfs_make_sims.r

# load functions
source('scripts/wfs.r')
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel)
	require(data.table)
	require(abc)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	require(iterators, lib.loc="/projects/cees/lib/R_packages/") # used with foreach
	require(foreach, lib.loc="/projects/cees/lib/R_packages/") # for parallel procssing using foreach()
	require(doMC, lib.loc="/projects/cees/lib/R_packages/")
	require(doParallel, lib.loc="/projects/cees/lib/R_packages/") # needed for foreach
	require(doSNOW, lib.loc="/projects/cees/lib/R_packages/") # set up SOCK cluster
	require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
	require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(locfit, lib.loc="/projects/cees/lib/R_packages/") # needed by abc
	require(abc.data, lib.loc="/projects/cees/lib/R_packages/") # needed by abc
	require(abc, lib.loc="/projects/cees/lib/R_packages/")
}


# load data
load('analysis/wfs_sims.rdata') # load the simulations: 'out' array. slow to load...
obs <- fread('analysis/LOF_07_to_LOF_S_14.w_obs_stats.txt') # the observations of Fsi and Fsd
dat <- fread('analysis/LOF_07_to_LOF_S_14.wfabc', skip=2) # the data on samples sizes and observed allele frequencies

# parameters
tol <- 0.001 # 10k from 10M sims

# prep data
sampsize <- dat[seq(1,nrow(dat),by=2),] # sample sizes in # of chromosomes
alcnts <- dat[seq(2,nrow(dat),by=2),] # allele counts in # chromosomes

# abc function
# out.ff is an ff matrix (one slice of out[,,i])
runabc <- function(locusnum, thistarg, out.ff, tol){

	# with f1samp as a summary stat
	# calc euclidean distance
	sum1 <- (out.ff['f1samp',] - thistarg[,f1samp.n])^2
	sum1 <- sum1 + (out.ff['fsdprime',] - thistarg[,fsd.n])^2
	sum1 <- sum1 + (out.ff['fsiprime',] - thistarg[,fsi.n])^2
	dst <- sqrt(sum1)
	
	# includes the effect of gwt in the tolerance
	gwt <- !is.na(dst) # for now, simply define gwt as value that aren't NA
	dst[!gwt] <- floor(max(dst[gwt])+10)
	
	# wt1 defines the region we're interested in 
	abstol <- quantile(dst,tol)
	wt1 <- dst < abstol
	if(sum(wt1)==0){ # can happen if min(dst) == abstol
		wt1 <- dst <= abstol
	}
 
	# set up return matrix
	ret <- matrix(data=c(rep(locusnum, sum(wt1)), out.ff['f1samp',wt1], out.ff['fsdprime',wt1], out.ff['fsiprime',wt1], out.ff['ne',wt1], out.ff['f1',wt1], out.ff['s',wt1]), ncol=7, byrow=FALSE)
	colnames(ret) <- c('locus', 'f1samp', 'fsdprime', 'fsiprime', 'ne', 'f1', 's')

	return(ret)
}

####################################################
### run calculations in parallel (parallel package)
####################################################
	# start cluster
#if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
#	cl <- makeCluster(detectCores()-1) # set up cluster on my mac (or another computer), using all but one core
#	clusterEvalQ(cl, library(abc.data))
#	clusterEvalQ(cl, library(locfit))
#	clusterEvalQ(cl, library(abc))
#}
#if(grepl('hpc.uio.no', Sys.info()["nodename"])){
#	cl <- makeCluster(20) # set up cluster on a cod node
#	clusterEvalQ(cl, library(abc.data, lib.loc='/projects/cees/lib/R_packages/'))
#	clusterEvalQ(cl, library(locfit, lib.loc='/projects/cees/lib/R_packages/'))
#	clusterEvalQ(cl, library(abc, lib.loc='/projects/cees/lib/R_packages/'))
#}
#
#clusterExport(cl, c('runabc'))
#
#	# do calculations
#lociperpart <- 400 # how many loci in each part (each part will be written to file)
#parts <- 1:ceiling(nrow(sampsize)/lociperpart) # 1 to number of parts
#print(paste(length(parts), "parts"))
##partinds <- rep(parts, rep(lociperpart, length(parts))) # index for each row to each part
#partinds <- vector('list', length(parts)); for(i in 1:length(partinds)) partinds[[i]] <- ((i-1)*lociperpart+1):min(c(i*lociperpart, nrow(sampsize))) # indices for the loci in each part
#ndigits <- nchar(as.character(nrow(sampsize))) # used for formatting file names

#for(partnum in parts){
##for(partnum in 1:2){
#	print(paste('partnum', partnum, 'of', length(parts)))
#
#	thisdat <- cbind(1:nrow(sampsize)[partinds[[partnum]]], sampsize[partinds[[partnum]],], alcnts[partinds[[partnum]],], obs[partinds[[partnum]],]) # the 
#	temp <- parApply(cl, thisdat, MARGIN=1, FUN=runabc, out=out) # a slow start: have to load out onto each node. output is annoying: a 60000 x nrow(obs) matrix. Rows are f1samp, fsdprime, fsiprime, ne, f1, s, f1samp, fsdprime, fsiprime, ne, f1, s, etc. (samples from the posterior). Columns are each locus.
#
#	# reformat CHECK THIS PART CAREFULLY
#	posts <- data.table(f1samp=t(temp[seq(1,nrow(temp),by=6),]), fsdprime=t(temp[seq(2,nrow(temp),by=6),]), fsiprime=t(temp[seq(3,nrow(temp),by=6),]), ne=t(temp[seq(4,nrow(temp),by=6),]), f1 = t(temp[seq(5,nrow(temp),by=6),]), s = t(temp[seq(6,nrow(temp),by=6),]))
#
#	minloc <- formatC(min(partinds[[partnum]]), width=ndigits, flag='0')
#	maxloc <- formatC(max(partinds[[partnum]]), width=ndigits, flag='0')
#	outfile <- paste('analysis/temp/wfs_abc_locus', minloc, '-', maxloc, '.csv.gz', sep='')
#	write.csv(posts, file=gzfile(outfile), row.names=FALSE) # write directly to gzipped file
#}
#

###################################################
### run calculations in parallel with foreach CURRENTLY USES TOO MUCH MEMORY
###################################################


# set up target statistics
	# stats for each locus
	targ <- cbind(sampsize[,.(V1, V2)], alcnts[,V1]/sampsize[,V1], obs[,.(Fsd, Fsi)]) # alcnt1, alcnt2, f1samp, fsd, fsi
	setnames(targ, 1:5, c('alcnt1', 'alcnt2', 'f1samp', 'fsd', 'fsi'))
	
	# scale by the appropriate summary statistic
	# need mean and sd for each summary stat in out
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
	
	# add mean and sd to targ to allow easy normalization
	setkey(targ, alcnt1, alcnt2)
	setkey(meansds, alcnt1, alcnt2)
	targ <- meansds[targ,]
	
	# normalize the target stats
	targ[,f1samp.n:=(f1samp-f1samp.mean)/f1samp.sd]
	targ[,fsd.n:=(fsd-fsd.mean)/fsd.sd]
	targ[,fsi.n:=(fsi-fsi.mean)/fsi.sd]

	# write out?
	save(targ, file='analysis/targ.rdata')	

# export out file from abc sims to ff data type for memory-efficient use across cluster
	outszs <- dimnames(out)[[3]] # sample sizes

	# will need to do this for all 193 combinations of sample sizes
	i = 1
		j <- which(dimnames(out)[[3]] == paste(meansds[i,.(alcnt1, alcnt2)], collapse=',')) # have to find the slice of out that matches the sample sizes for this row
		out.ff <- ff(out[,,j], dim=dim(out[,,j]), dimnames=dimnames(out[,,j]), filename=paste(getwd(), '/analysis/temp/wfs_sims_ff', paste(meansds[i,.(alcnt1, alcnt2)], collapse=','), sep=''), finalizer='close')

		# scale everything
		out.ff['f1samp',] <- (out.ff['f1samp',]-meansds[i,f1samp.mean])/meansds[i,f1samp.sd]
		out.ff['fsdprime',] <- (out.ff['fsdprime',]-meansds[i,fsd.mean])/meansds[i,fsd.sd]
		out.ff['fsiprime',] <- (out.ff['fsiprime',]-meansds[i,fsi.mean])/meansds[i,fsi.sd]
	}

## parts before this can be in a data prep script

# set up data chunks to operate on and write out separately
lociperpart <- 2 # how many loci in each part (each part will be written to file)
parts <- 1:ceiling(nrow(targ)/lociperpart) # 1 to number of parts
print(paste(length(parts), "parts"))
#partinds <- rep(parts, rep(lociperpart, length(parts))) # index for each row to each part
partinds <- vector('list', length(parts)); for(i in 1:length(partinds)) partinds[[i]] <- ((i-1)*lociperpart+1):min(c(i*lociperpart, nrow(targ))) # indices for the loci in each part
ndigits <- nchar(as.character(nrow(targ))) # used for formatting file names

# set up cluster
#registerDoMC(2) # register cores. doesn't clean up after itself (creates 15 new cores on every loop)
.libPaths('/projects/cees/lib/R_packages/')
cl <- makeCluster(2, type='SOCK')
registerDoSNOW(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths()) # add libpaths to each node

#for(partnum in parts){
for(partnum in 1:2){
	print(paste('partnum', partnum, 'of', length(parts)))
	thisdat <- iter(cbind((1:nrow(sampsize))[partinds[[partnum]]], targ[partinds[[partnum]],]), by='row') # the iterator object: this part of the loop's chunk of data
	
	results <- foreach(i=thisdat, .combine=rbind, .packages=c('ff', 'doParallel', 'data.table')) %dopar% {
		# have to pick the correct ff object here (there will be 193)
#		j <- which(dimnames(out)[[3]] == paste(sampsize[i,], collapse=',')) # have to find the slice of out that matches the sample sizes for this row
		# then run abc
		runabc(locusnum=i[,V1], thistarg=i[,.(f1samp.n, fsd.n, fsi.n)], out.ff, tol=tol) 
	}

	minloc <- formatC(min(partinds[[partnum]]), width=ndigits, flag='0')
	maxloc <- formatC(max(partinds[[partnum]]), width=ndigits, flag='0')
	outfile <- paste('analysis/temp/wfs_abc_locus', minloc, '-', maxloc, '.csv.gz', sep='')
	write.csv(results, file=gzfile(outfile), row.names=FALSE) # write directly to gzipped file
}

# save
#save(posts, file='analysis/wfs_abc_posts.rdata')

# stop cluster
stopCluster(cl)
