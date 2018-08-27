# OBSOLETE: NOW USE wfs_abc_function.r
# ABC on each locus from ABC lookup table to estimate s
# run after wfs_make_simsCAN.r
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


# load target and sample size data
targ <- fread('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', header=TRUE)
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locusnum := 1:nrow(targ)]


# parameters
tol <- 0.001 # 10k from 10M sims
lociperpart <- 100 # how many loci in each part (each part will be written to file)

# abc function
# thisout.ff is an ff matrix (one slice of out[,,i])
# assumes that thisout.ff['f1samp' and 'f2samp' are not normalized
runabc <- function(locusnum, thistarg, thisout.ff, tol){

	# with f2samp-f1samp as a summary stat
	# calc euclidean distance
	d1 <- thisout.ff['f2samp',] - thisout.ff['f1samp',]
	d2 <- thistarg[,f2samp] - thistarg[,f1samp]
	dst <- sqrt((d1-d2)^2)

	# includes the effect of gwt in the tolerance
	gwt <- !is.na(dst) # for now, simply define gwt as values that aren't NA
	dst[!gwt] <- floor(max(dst[gwt])+10)

	# wt1 defines the region we're interested in 
	abstol <- quantile(dst,tol)
	wt1 <- dst < abstol
	if(sum(wt1)==0){ # can happen if min(dst) == abstol
		wt1 <- dst <= abstol
	}

	# set up return matrix
	ret <- matrix(data=c(rep(locusnum, sum(wt1)), thisout.ff['f1samp',wt1], thisout.ff['f2samp',wt1], thisout.ff['ne',wt1], thisout.ff['f1',wt1], thisout.ff['f2',wt1], thisout.ff['s',wt1]), ncol=7, byrow=FALSE)
	colnames(ret) <- c('locus', 'f1samp', 'f2samp', 'ne', 'f1', 'f2', 's')

	return(ret)
}



###################################################
### run calculations in parallel with foreach
###################################################

# find if any loci of this sample size have already been run (for checkpointing)
existingfiles <- list.files(path='analysis/temp', pattern=paste('wfs_abc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''))
print(paste('found', length(existingfiles), 'existing files relevant to this sample size'))
existingrngs <- strsplit(gsub(paste('wfs_abc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus|.csv.gz', sep=''), '', existingfiles), split='-') # extract just the locus ranges
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
	ffnm <- paste('analysis/temp/wfs_simsCAN_ff', paste(myalcnt1, myalcnt2, sep=','), sep='')
	ffload(ffnm, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff

	# loop through each chunk of data, writing each to file
	for(partnum in parts){
		print(paste('partnum', partnum, 'of', length(parts)))
		thisdat <- iter(targ[locusnum %in% partinds[[partnum]],], by='row') # the iterator object: this part of the loop's chunk of data

		results <- foreach(i=thisdat, .combine=rbind, .packages=c('ff', 'doParallel', 'data.table')) %dopar% {
			runabc(locusnum=i[,locusnum], thistarg=i[,.(f1samp, f2samp)], thisout.ff=thisout.ff, tol=tol) # run abc analysis. returns vector of selected simulation parameters
		}

		minloc <- formatC(min(partinds[[partnum]]), width=ndigits, flag='0')
		maxloc <- formatC(max(partinds[[partnum]]), width=ndigits, flag='0')
		outfile <- paste('analysis/temp/wfs_abc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus', minloc, '-', maxloc, '.csv.gz', sep='')
		write.csv(results, file=gzfile(outfile), row.names=FALSE) # write directly to gzipped file. has f1samp, f2samp, plus ne, f1, s
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