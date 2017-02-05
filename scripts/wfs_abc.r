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


# prep data
sampsize <- dat[seq(1,nrow(dat),by=2),] # sample sizes in # of chromosomes
alcnts <- dat[seq(2,nrow(dat),by=2),] # allele counts in # chromosomes

# abc function
# out.ff is an ff matrix (one slice of out[,,i])
runabc <- function(allobsrow, out.ff){
	locusnum <- as.numeric(allobsrow)[1]
	sampsizerow <- as.numeric(allobsrow)[2:3]
	alcntsrow <- as.numeric(allobsrow)[4:5]
	obsrow <- as.numeric(allobsrow)[6:11]

#	j <- which(dimnames(out)[[3]] == paste(sampsizerow, collapse=',')) # find the slice of out that has simulations with the equivalent sample size
#	if(length(j)==1){
	if(length(dim(out.ff))==2){
		targ <- c(alcntsrow[1]/sampsizerow[1], obsrow[2:1]) # f1samp, fsd, fsi
		pams <- as.data.table(t(out.ff[c(1,2,3),])) # ne, f1, s
		stats <- as.data.frame(t(out[c(6,8,9),])) # f1samp, fsd, fsi
#		pams <- as.data.frame(t(out[c(1,2,3),,j])) # ne, f1, s
#		stats <- as.data.frame(t(out[c(6,8,9),,j])) # f1samp, fsd, fsi


		# with f1samp as a summary stat
#		post <- abc(target=targ, param=pams, sumstat=stats, tol=0.001, method='rejection')

		# do by hand
    
		# calc euclidean distance

			sum1 <- 0
			for(j in 1:nss){
				sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
		   }
		   dst <- sqrt(sum1)
		# includes the effect of gwt in the tolerance
			dst[!gwt] <- floor(max(dst[gwt])+10)
	

		# wt1 defines the region we're interested in 
			abstol <- quantile(dst,tol)
			wt1 <- dst < abstol
 


		ret <- matrix(data=c(rep(locusnum, nrow(post$ss)), post$ss[,1], post$ss[,2], post$ss[,3], post$unadj.values[,1], post$unadj.values[,2], post$unadj.values[,3]), ncol=7, byrow=FALSE)
		colnames(ret) <- c('locus', 'f1samp', 'fsdprime', 'fsiprime', 'ne', 'f1', 's')

		return(ret)
	} else {
		warning(paste("Didn't have sample size sims for row", locusnum))
	}
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
# set up cluster
#registerDoMC(15) # register cores. seems to use too much memory (creates 15 new cores on every loop)
cl <- makeCluster(15, type='SOCK')
registerDoSNOW(cl)

# set up data chunks to operate on and write out separately
lociperpart <- 15 # how many loci in each part (each part will be written to file)
parts <- 1:ceiling(nrow(sampsize)/lociperpart) # 1 to number of parts
print(paste(length(parts), "parts"))
#partinds <- rep(parts, rep(lociperpart, length(parts))) # index for each row to each part
partinds <- vector('list', length(parts)); for(i in 1:length(partinds)) partinds[[i]] <- ((i-1)*lociperpart+1):min(c(i*lociperpart, nrow(sampsize))) # indices for the loci in each part
ndigits <- nchar(as.character(nrow(sampsize))) # used for formatting file names

# set up target statistics
targ <- cbind(alcnts[,V1]/sampsize[,V1], obs[,.(Fsd, Fsi)]) # f1samp, fsd, fsi
	# scale by the appropriate summary statistic
	# need mean and sd for each summary stat in out
	uniqsamps <- apply(do.call('rbind', strsplit(dimnames(out)[[3]], split=',')), 2, as.numeric) # extract unique sample sizes
	meansds <- data.table(alcnt1=uniqsamps[,1], alcnt2=uniqsamps[,2])
	meansds[,f1samp.mean:=NA]
	meansds[,f1samp.sd:=NA]
	meansds[,fsd.mean:=NA]
	meansds[,fsd.sd:=NA]
	meansds[,fsi.mean:=NA]
	meansds[,fsi.sd:=NA]
	for(i in 1:nrow(meansds)){ 
		print(i)
		j <- which(dimnames(out)[[3]] == paste(meansds[i,.(alcnt1, alcnt2)], collapse=',')) # have to find the slice of out that matches the sample sizes for this row
		meansds[i,f1samp.mean:=mean(out['f1samp',,j])]
		meansds[i,f1samp.sd:=sd(out['f1samp',,j])]
		meansds[i,fsd.mean:=mean(out['fsdprime',,j])]
		meansds[i,fds.sd:=sd(out['fsdprime',,j])]
		meansds[i,fsi.mean:=mean(out['fsiprime',,j])]
		meansds[i,fsi.sd:=sd(out['fsiprime',,j])]
	}
	
for(i in 1:nrow(out.ff)){
	targ[j] <- normalise(targ[j],sumstat[,j][gwt])
}

# export out file from abc sims to ff data type for memory-efficient use across cluster
outszs <- dimnames(out)[[3]] # sample sizes

out.ff <- ff(out[,,5], dim=dim(out[,,5]), filename=paste(getwd(), '/analysis/temp/wfs_sims_ff', paste(sampsize[1,], collapse=','), sep=''), finalizer='delete')
	# scale everything
	for(i in 1:nrow(out.ff)){
		out.ff[i,] <- normalise(sumstat[,j],sumstat[,j][gwt])
	}

#for(partnum in parts){
for(partnum in 1:2){
print(paste('partnum', partnum, 'of', length(parts)))
	thisdat <- iter(cbind((1:nrow(sampsize))[partinds[[partnum]]], targ[partinds[[partnum]],]), by='row') # the iterator object: this part of the loop's chunk of data
	
	results <- foreach(i=thisdat, .combine=rbind, .packages=c('ff', 'abc', 'abc.data')) %dopar% {
		j <- which(dimnames(out)[[3]] == paste(sampsize[i,], collapse=',')) # have to find the slice of out that matches the sample sizes for this row
		# have to pick the ff object here
		runabc(i, out.ff) 
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
