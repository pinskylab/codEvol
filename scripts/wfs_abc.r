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



###################################################
### run calculations in parallel with foreach
###################################################

# set up cluster
#registerDoMC(2) # register cores. doesn't clean up after itself (creates new cores on every loop without removing old ones)
.libPaths('/projects/cees/lib/R_packages/')
cl <- makeCluster(ncores, type='SOCK')
registerDoSNOW(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths()) # add libpaths to each node

# set up chunks of the same sample size to operate on and write out separately
samplepartinds <- vector('list', nrow(meansds)) # hold the locus numbers for each sample size
names(samplepartinds) <- meansds[,paste(alcnt1, alcnt2, sep=',')]
for(i in 1:length(samplepartinds)) samplepartinds[[i]] <- targ[alcnt1==meansds[i,alcnt1] & alcnt2==meansds[i,alcnt2],locusnum] # locus numbers in each sample size chunk
sampleparts.n <- sapply(samplepartinds, length) # number of loci in each samplesize chunk

ndigits <- nchar(as.character(nrow(targ))) # used for formatting file names

for(k in 1:length(samplepartinds)){
	print(paste('samplepart', k, 'of', length(samplepartinds), ':', names(samplepartinds)[k], sampleparts.n[k], 'loci', Sys.time()))

	# set up data chunks to operate on and write out separately within this chunk of sample sizes
	parts <- 1:ceiling(length(samplepartinds[[k]])/lociperpart) # 1 to number of parts for this samplepart
	partinds <- vector('list', length(parts))
	for(i in 1:length(partinds)) partinds[[i]] <- samplepartinds[[k]][((i-1)*lociperpart+1):min(c(i*lociperpart, length(samplepartinds[[k]])))] # indices for the loci in each part

	ffnm <- paste('analysis/temp/wfs_sims_ff', names(samplepartinds)[k], sep='')
	ffload(ffnm, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff

	for(partnum in parts){
		print(paste('partnum', partnum, 'of', length(parts)))
		thisdat <- iter(targ[locusnum %in% partinds[[partnum]],], by='row') # the iterator object: this part of the loop's chunk of data
	
		results <- foreach(i=thisdat, .combine=rbind, .packages=c('ff', 'doParallel', 'data.table')) %dopar% {
			runabc(locusnum=i[,locusnum], thistarg=i[,.(f1samp.n, fsd.n, fsi.n)], thisout.ff=thisout.ff, tol=tol) # run abc analysis. returns vector of selected simulation parameters
		}

		minloc <- formatC(min(partinds[[partnum]]), width=ndigits, flag='0')
		maxloc <- formatC(max(partinds[[partnum]]), width=ndigits, flag='0')
		outfile <- paste('analysis/temp/wfs_abc_sampsize', names(samplepartinds)[k], '_locus', minloc, '-', maxloc, '.csv.gz', sep='')
		write.csv(results, file=gzfile(outfile), row.names=FALSE) # write directly to gzipped file. has normalized f1samp, fsdprime, fsiprime, plus ne, f1, s
		print(paste('wrote', outfile))
	}

	# remove ff temporary files
	delete(thisout.ff)
	rm(thisout.ff)
}

# stop cluster
stopCluster(cl)

