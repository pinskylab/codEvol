# ABC on each locus from ABC lookup table to estimate s
# run after wfs_make_sims.r and wfs_process_sims.r

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
}


# load data
load('analysis/wfs_targ.rdata') # targets normalized
load('analysis/wfs_sims_meansds.rdata') # useful as a list of sample sizes, which help when loading the ff files


# parameters
tol <- 0.001 # 10k from 10M sims
lociperpart <- 10 # how many loci in each part (each part will be written to file)
ncores <- 10 # number of cores to use

# abc function
# thisout.ff is an ff matrix (one slice of out[,,i])
runabc <- function(locusnum, thistarg, thisout.ff, tol){

	# with f1samp as a summary stat
	# calc euclidean distance
	sum1 <- (thisout.ff['f1samp',] - thistarg[,f1samp.n])^2
	sum1 <- sum1 + (thisout.ff['fsdprime',] - thistarg[,fsd.n])^2
	sum1 <- sum1 + (thisout.ff['fsiprime',] - thistarg[,fsi.n])^2
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
	ret <- matrix(data=c(rep(locusnum, sum(wt1)), thisout.ff['f1samp',wt1], thisout.ff['fsdprime',wt1], thisout.ff['fsiprime',wt1], thisout.ff['ne',wt1], thisout.ff['f1',wt1], thisout.ff['s',wt1]), ncol=7, byrow=FALSE)
	colnames(ret) <- c('locus', 'f1samp', 'fsdprime', 'fsiprime', 'ne', 'f1', 's')

	return(ret)
}

# needs work! clunky because thisout.ff hasn't been reprocessed
runabcf1f2 <- function(locusnum, thistarg, thisout.ff, tol, f1samp.sd, f1samp.mean){

	# with f1samp as a summary stat
	# calc euclidean distance
	sum1 <- (thisout.ff['f1samp',]*f1samp.sd+f1samp.mean - thistarg[,f1samp])^2
	sum1 <- sum1 + (thisout.ff['f2samp',] - thistarg[,f2samp])^2
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
	ret <- matrix(data=c(rep(locusnum, sum(wt1)), thisout.ff['f1samp',wt1]*f1samp.sd+f1samp.mean, thisout.ff['f2samp',wt1], thisout.ff['fsdprime',wt1], thisout.ff['fsiprime',wt1], thisout.ff['ne',wt1], thisout.ff['f1',wt1], thisout.ff['s',wt1]), ncol=8, byrow=FALSE)
	colnames(ret) <- c('locus', 'f1samp', 'f2samp', 'fsdprime', 'fsiprime', 'ne', 'f1', 's')

	return(ret)
}



###################################################
### run calculations in parallel with foreach
###################################################

locnum <- 7153
sampsize <- '38,48'
meanind <- which(meansds$alcnt1 == 38 & meansds$alcnt2==48)

ffnm <- paste('analysis/temp/wfs_sims_ff', sampsize, sep='')
ffload(ffnm, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff

i <- targ[locnum,] # the input object
	
results.fsi <- runabc(locusnum=i[,locusnum], thistarg=i[,.(f1samp.n, fsd.n, fsi.n)], thisout.ff=thisout.ff, tol=tol) # run abc analysis with f1samp and fsi/fsd (the old way)

results.f1f2 <- runabcf1f2(locusnum=i[,locusnum], thistarg=i[,.(f1samp, f2samp)], thisout.ff=thisout.ff, tol=tol, f1samp.mean=meansds$f1samp.mean[meanind], f1samp.sd=meansds$f1samp.sd[meanind]) # run abc with unscaled f1samp and f2samp (the possible new way)


	summary(results.fsi)
	summary(results.f1f2)

	# remove ff temporary files
	delete(thisout.ff)
	rm(thisout.ff)
