# ABC on each locus from ABC lookup table to estimate s
# run after wfs_make_sims.r and wfs_process_sims.r
# This version set up to run a single sample size, taken as command line arguments

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args)<2) {
  stop("Have to specify myalcnt1 and myalcnt2", call.=FALSE)
} else if (length(args)==2){
	maxcores <- 16 # default maximum cores (for abel)
} else if (length(args)>2) {
	maxcores <- args[3]
}
myalcnt1 <- args[1]
myalcnt2 <- args[2]

print(paste('myalcnt1', myalcnt1, 'myalcnt2', myalcnt2, 'maxcores', maxcores))
print(Sys.info()["nodename"])

# load functions: assume this is run on a cod or abel node
require(parallel, lib.loc="/projects/cees/lib/R_packages/")
require(iterators, lib.loc="/projects/cees/lib/R_packages/") # used with foreach
require(foreach, lib.loc="/projects/cees/lib/R_packages/") # for parallel procssing using foreach()
require(doMC, lib.loc="/projects/cees/lib/R_packages/")
require(doParallel, lib.loc="/projects/cees/lib/R_packages/") # needed for foreach
require(doSNOW, lib.loc="/projects/cees/lib/R_packages/") # set up SOCK cluster
require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
require(data.table, lib.loc="/projects/cees/lib/R_packages/")


# load data
load('analysis/wfs_targ.rdata') # targets normalized
load('analysis/wfs_sims_meansds.rdata') # useful as a list of sample sizes, which help when loading the ff files


# parameters
tol <- 0.001 # 10k from 10M sims
lociperpart <- 100 # how many loci in each part (each part will be written to file)

# abc function
# thisout.ff is an ff matrix (one slice of out[,,i])
runabc <- function(locusnum, thistarg, thisout.ff, tol, rej=TRUE){

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

	if(rej){ # rejection method
		# set up return matrix
		ret <- matrix(data=c(rep(locusnum, sum(wt1)), thisout.ff['f1samp',wt1], thisout.ff['fsdprime',wt1], thisout.ff['fsiprime',wt1], thisout.ff['ne',wt1], thisout.ff['f1',wt1], thisout.ff['s',wt1]), ncol=7, byrow=FALSE)
		colnames(ret) <- c('locus', 'f1samp', 'fsdprime', 'fsiprime', 'ne', 'f1', 's')

		return(ret)
	} else {
		print(paste('Regression-adjusting locus', locusnum))
		regwt <- 1-dst[wt1]^2/abstol^2
		if(sd(regwt)==0){
			regwt <- rep(1,length(regwt)) # to catch all(regwt==0)
			warning(paste('No variation in regwt in the tolerance region for locus', locusnum))
		}
		fit1 <- lsfit(t(thisout.ff[c('f1samp', 'fsdprime', 'fsiprime'),wt1]),t(thisout.ff[c('ne', 'f1', 's'),wt1]),wt=regwt) # may return error that x matrix is colinear. probably means the model is useless.
		predmean.ne <- fit1$coeff[,1] %*% c(1,as.numeric(thistarg))
		predmean.f1 <- fit1$coeff[,2] %*% c(1,as.numeric(thistarg))
		predmean.s <- fit1$coeff[,3] %*% c(1,as.numeric(thistarg))
		
		ret <- matrix(data=c(rep(locusnum, sum(wt1)), thisout.ff['f1samp',wt1], thisout.ff['fsdprime',wt1], thisout.ff['fsiprime',wt1], fit1$residuals[,1]+predmean.ne, fit1$residuals[,2]+predmean.f1, fit1$residuals[,3]+predmean.s), ncol=7, byrow=FALSE) # the parameters are regression-adjusted values
		colnames(ret) <- c('locus', 'f1samp', 'fsdprime', 'fsiprime', 'ne', 'f1', 's')
	
		return(ret)
	}
}



###################################################
### run calculations in parallel with foreach
###################################################

# find if any loci of this sample size have already been run (for checkpointing)
existingfiles <- list.files(path='analysis/temp', pattern=paste('wfs_abcregr_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''))
print(paste('found', length(existingfiles), 'existing files relevant to this sample size'))
existingrngs <- strsplit(gsub(paste('wfs_abcregr_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus|.csv.gz', sep=''), '', existingfiles), split='-') # extract just the locus ranges
existingrngs <- lapply(existingrngs, as.numeric)
existingloci <- numeric(0)
if(length(existingrngs)>0) for(i in 1:length(existingrngs)) existingloci <- c(existingloci, existingrngs[[i]][1]:existingrngs[[i]][2])

# find the chunk of loci of the specified sample size to operate on
samplepartinds <- targ[alcnt1==myalcnt1 & alcnt2==myalcnt2,locusnum] # locus numbers in each sample size chunk
orig.n <- length(samplepartinds)
samplepartinds <- setdiff(samplepartinds, existingloci) # remove loci already run (from existing files)
sampleparts.n <- length(samplepartinds) # number of loci in the samplesize chunk
print(paste('Will run abc for', sampleparts.n, 'loci of', orig.n, 'loci originally'))

ndigits <- nchar(as.character(nrow(targ))) # used for formatting file names

if(sampleparts.n > 0){
	# set up cluster
	ncores <- min(maxcores, sampleparts.n) # don't set up any more cores than we need for this batch of loci
	print(paste('Using', ncores, 'cores'))
	.libPaths('/projects/cees/lib/R_packages/')
	cl <- makeCluster(rep('localhost', ncores), type='SOCK') # make the cluster on the localhost
	registerDoSNOW(cl) # register the cluster so accessible to foreach
	clusterCall(cl, function(x) .libPaths(x), .libPaths()) # add libpaths to each node


	# set up data chunks to operate on and write out separately within this chunk of sample sizes
	parts <- 1:ceiling(sampleparts.n/lociperpart) # 1 to number of parts for this samplepart
	partinds <- vector('list', length(parts))
	for(i in 1:length(partinds)) partinds[[i]] <- samplepartinds[((i-1)*lociperpart+1):min(c(i*lociperpart, sampleparts.n))] # indices for the loci in each part

	# load the appropriate file of abc simulations
	ffnm <- paste('analysis/temp/wfs_sims_ff', paste(myalcnt1, myalcnt2, sep=','), sep='')
	ffload(ffnm, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff

	# loop through each chunk of data, writing each to file
	for(partnum in parts){
		print(paste('partnum', partnum, 'of', length(parts)))
		thisdat <- iter(targ[locusnum %in% partinds[[partnum]],], by='row') # the iterator object: this part of the loop's chunk of data

		results <- foreach(i=thisdat, .combine=rbind, .packages=c('ff', 'doParallel', 'data.table')) %dopar% {
			runabc(locusnum=i[,locusnum], thistarg=i[,.(f1samp.n, fsd.n, fsi.n)], thisout.ff=thisout.ff, tol=tol, rej=FALSE) # run abc analysis. returns vector of selected, regression-adjusted simulation parameters
		}

		minloc <- formatC(min(partinds[[partnum]]), width=ndigits, flag='0')
		maxloc <- formatC(max(partinds[[partnum]]), width=ndigits, flag='0')
		outfile <- paste('analysis/temp/wfs_abcregr_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus', minloc, '-', maxloc, '.csv.gz', sep='')
		write.csv(results, file=gzfile(outfile), row.names=FALSE) # write directly to gzipped file. has normalized f1samp, fsdprime, fsiprime, plus ne, f1, s
		print(paste('wrote', outfile))
	}

	# remove ff temporary files
	delete(thisout.ff)
	rm(thisout.ff)

	# stop cluster
	stopCluster(cl)
} else {
	print('Exiting because all loci had already been run')
}