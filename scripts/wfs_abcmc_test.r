# ABC on each locus from ABC lookup table to estimate s
# run after wfs_make_sims.r and wfs_process_sims.r

# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
	require(ff)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
	require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# load data
load('analysis/wfs_targ.rdata') # target statistics


# parameters
tol <- 0.001 # 10k from 10M sims


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

locnum <- 73; sampsize <- '50,44'

# load data
	# full model
ffnm1 <- paste('analysis/temp/wfs_sims_ff', sampsize, sep='')
ffload(ffnm1, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff
thisout.ff1 <- thisout.ff
thisout.ff1 <- thisout.ff1[rownames(thisout.ff1) %in% c('ne', 'f1', 's', 'gen', 'f2', 'f1samp', 'f2samp'),] # slow!

	# remove ff object (temp file still exists)
	rm(thisout.ff)

	# null model
ffnm2 <- paste('analysis/temp/wfs_simsnull_ff', sampsize, sep='')
ffload(ffnm2, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff
thisout.ff2 <- thisout.ff

	# remove ff object (temp file still exists)
	rm(thisout.ff)

# check whether f1samp was normalized
summary(thisout.ff1['f1samp',]) 
summary(thisout.ff2['f1samp',])


thistarg <- targ[locnum,] # the input object
	
# run abc
results <- runabc_mc(locusnum=thistarg[,locusnum], thistarg=thistarg[,.(f1samp, f2samp)], thisout.ff1=thisout.ff1, thisout.ff2=thisout.ff2, tol=tol) # run abc model choice


	summary(results)
	sum(results[,2]==1)/nrow(results)
	sum(results[,2]==2)/nrow(results)

	inds1 <- results[,2]==1
	inds2 <- results[,2]==2
	summary(results[inds1,])
	summary(results[inds2,])

	# remove ff temporary files
	delete(thisout.ff1)
	rm(thisout.ff1)
	delete(thisout.ff2)
	rm(thisout.ff2)
