# ABC model choice on each locus from ABC lookup table to estimate s
# run after wfs_make_sims.r/wfs_process_sims.r and wfs_make_sims_null.r
# This version set up to run a single sample size, taken as command line arguments

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args)<2) {
  stop("Have to specify myalcnt1 and myalcnt2", call.=FALSE)
} else if (length(args)==2){
	maxcores <- 16 # default maximum cores (for abel)
} else if (length(args)>2) {
	maxcores <- as.numeric(args[3])
}
myalcnt1 <- as.numeric(args[1])
myalcnt2 <- as.numeric(args[2])

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

# parameters
tol <- 0.001 # 10k from 10M sims
lociperpart <- 100 # how many loci in each part (each part will be written to file)


# ABC model choice
# using f1samp and f2samp
# assumes neither have been normalized in thisout.ff
# expects thisout.ff1 is the alternative model, thisout.ff2 is the null model
runabc_mc <- function(locusnum, thistarg, thisout.ff1, thisout.ff2, tol){

	# calc euclidean distance for alternative model
	sum1 <- (thisout.ff1['f1samp',] - thistarg[,f1samp])^2
	sum1 <- sum1 + (thisout.ff1['f2samp',] - thistarg[,f2samp])^2
	dst1 <- sqrt(sum1)

	# calc euclidean distance for null model
	sum2 <- (thisout.ff2['f1samp',] - thistarg[,f1samp])^2
	sum2 <- sum1 + (thisout.ff2['f2samp',] - thistarg[,f2samp])^2
	dst2 <- sqrt(sum2)
	
	# includes the effect of gwt in the tolerance
	gwt <- !is.na(dst1) # for now, simply define gwt as value that aren't NA
	dst1[!gwt] <- floor(max(dst1[gwt])+10)
	gwt <- !is.na(dst2) # for now, simply define gwt as value that aren't NA
	dst2[!gwt] <- floor(max(dst2[gwt])+10)
	
	# wt1 defines the region we're interested in 
	abstol <- quantile(c(dst1, dst2),tol)
	wt1 <- dst1 < abstol
	wt2 <- dst2 < abstol
	if(sum(c(wt1,wt2))==0){ # can happen if min(dst) == abstol
		wt1 <- dst1 <= abstol
		wt2 <- dst2 <= abstol
	}
	
	# figure out models chosen
	mod <- c(rep(1,sum(wt1)), rep(2, sum(wt2)))
 
	# set up return matrix
	ret <- matrix(data=c(rep(locusnum, sum(wt1)+sum(wt2)), mod, c(thisout.ff1['f1samp',wt1], thisout.ff2['f1samp',wt2]), c(thisout.ff1['f2samp',wt1], thisout.ff2['f2samp',wt2]), c(thisout.ff1['f1',wt1], thisout.ff2['f1',wt2]), c(thisout.ff1['f2',wt1], thisout.ff2['f2',wt2]), c(thisout.ff1['ne',wt1], thisout.ff2['ne',wt2]), c(thisout.ff1['s',wt1], thisout.ff2['s',wt2])), ncol=8, byrow=FALSE)
	colnames(ret) <- c('locus', 'model', 'f1samp', 'f2samp', 'f1', 'f2', 'ne', 's')

	return(ret)
}




###################################################
### run calculations in parallel with foreach
###################################################

# find if any loci of this sample size have already been run (for checkpointing)
existingfiles <- list.files(path='analysis/temp', pattern=paste('wfs_abcmc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''))
print(paste('found', length(existingfiles), 'existing files relevant to this sample size'))
existingrngs <- strsplit(gsub(paste('wfs_abcmc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus|.csv.gz', sep=''), '', existingfiles), split='-') # extract just the locus ranges
existingrngs <- lapply(existingrngs, as.numeric)
existingloci <- numeric(0)
if(length(existingrngs)>0) for(i in 1:length(existingrngs)) existingloci <- c(existingloci, existingrngs[[i]][1]:existingrngs[[i]][2])

# find the chunk of loci of the specified sample size to operate on
samplepartinds <- targ[alcnt1==myalcnt1 & alcnt2==myalcnt2,locusnum] # locus numbers in each sample size chunk
orig.n <- length(samplepartinds)
samplepartinds <- setdiff(samplepartinds, existingloci) # remove loci already run (from existing files)
sampleparts.n <- length(samplepartinds) # number of loci in the samplesize chunk
print(paste('Will run abc model choice for', sampleparts.n, 'loci of', orig.n, 'loci originally'))

ndigits <- nchar(as.character(nrow(targ))) # used for formatting file names

if(sampleparts.n > 0){
	# set up cluster
	ncores <- min(c(maxcores, sampleparts.n)) # don't set up any more cores than we need for this batch of loci
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
		# full model
	ffnm1 <- paste('analysis/temp/wfs_sims_ff', paste(myalcnt1, myalcnt2, sep=','), sep='')
	ffload(ffnm1, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff
	thisout.ff1 <- thisout.ff

		# remove ff object (temp file still exists)
		rm(thisout.ff)

		# null model
	ffnm2 <- paste('analysis/temp/wfs_simsnull_ff', paste(myalcnt1, myalcnt2, sep=','), sep='')
	ffload(ffnm2, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff
	thisout.ff2 <- thisout.ff

		# remove ff object (temp file still exists)
		rm(thisout.ff)

	# loop through each chunk of data, writing each to file
	for(partnum in parts){
		print(paste('partnum', partnum, 'of', length(parts)))
		thisdat <- iter(targ[locusnum %in% partinds[[partnum]],], by='row') # the iterator object: this part of the loop's chunk of data

		results <- foreach(i=thisdat, .combine=rbind, .packages=c('ff', 'doParallel', 'data.table')) %dopar% {
			runabc_mc(locusnum=i[,locusnum], thistarg=i[,.(f1samp, f2samp)], thisout.ff1=thisout.ff1, thisout.ff2=thisout.ff2, tol=tol) # run abc model choice analysis. returns vector of selected simulation parameters
		}

		minloc <- formatC(min(partinds[[partnum]]), width=ndigits, flag='0')
		maxloc <- formatC(max(partinds[[partnum]]), width=ndigits, flag='0')
		outfile <- paste('analysis/temp/wfs_abcmc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus', minloc, '-', maxloc, '.csv.gz', sep='')
		write.csv(results, file=gzfile(outfile), row.names=FALSE) # write directly to gzipped file. has normalized f1samp, fsdprime, fsiprime, plus ne, f1, s
		print(paste('wrote', outfile))
	}

	# remove ff temporary files
	delete(thisout.ff1)
	rm(thisout.ff1)
	delete(thisout.ff2)
	rm(thisout.ff2)

	# stop cluster
	stopCluster(cl)
} else {
	print('Exiting because all loci had already been run')
}