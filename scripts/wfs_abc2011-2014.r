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
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(locfit, lib.loc="/projects/cees/lib/R_packages/") # needed by abc
	require(abc.data, lib.loc="/projects/cees/lib/R_packages/") # needed by abc
	require(abc, lib.loc="/projects/cees/lib/R_packages/")
}


# load data
load('analysis/wfs_sims_notrim.rdata') # load the simulations: 'out' array. may be slow to load...
obs <- fread('analysis/LOF_S_11_LG03_to_LOF_S_14_LG03_notrim.w_obs_stats.txt') # the observations of Fsi and Fsd
dat <- fread('analysis/LOF_S_11_LG03_to_LOF_S_14_LG03_notrim.wfabc', skip=2) # the data on samples sizes and observed allele frequencies


# prep data
sampsize <- dat[seq(1,nrow(dat),by=2),] # sample sizes in # of chromosomes
alcnts <- dat[seq(2,nrow(dat),by=2),] # allele counts in # chromosomes

# abc function
runabc <- function(allobsrow, out){
	sampsizerow <- as.numeric(allobsrow)[1:2]
	alcntsrow <- as.numeric(allobsrow)[3:4]
	obsrow <- as.numeric(allobsrow)[5:10]

	j <- which(dimnames(out)[[3]] == paste(sampsizerow, collapse=',')) # find the slice of out that has simulations with the equivalent sample size
	if(length(j)==1){
		targ <- c(alcntsrow[1]/sampsizerow[1], obsrow[2:1]) # f1samp, fsd, fsi
		pams <- as.data.frame(t(out[c(1,2,3),,j])) # ne, f1, s
		stats <- as.data.frame(t(out[c(6,8,9),,j])) # f1samp, fsd, fsi
#		stats$fsdprime[is.na(stats$fsdprime)] <- 0
#		stats$fsiprime[is.na(stats$fsiprime)] <- 0

		# with f1samp as a summary stat
		post <- abc(target=targ, param=pams, sumstat=stats, tol=0.001, method='rejection')
		#post3 <- abc(target=targ, param=pams, sumstat=stats, tol=0.001, method='ridge', transf=c('log', 'logit', 'none'), logit.bounds=matrix(c(0,0,0,1,1,1),ncol=2)) # fails: infinite or missing values in 'x'
		#post3 <- abc(target=targ, param=pams, sumstat=stats, tol=0.01, method='ridge', transf=c('log', 'logit', 'none'), logit.bounds=matrix(c(0,0,0,1,1,1),ncol=2)) 
		#post4 <- abc(target=targ, param=pams, sumstat=stats, tol=0.001, method='neuralnet', transf=c('log', 'logit', 'none'), logit.bounds=matrix(c(0,0,0,1,1,1),ncol=2)) # can be slow
	
		# without f1samp
		#post_alt <- abc(target=targ[2:3], param=pams, sumstat=stats[,2:3], tol=0.001, method='rejection') # about the same answer for s, despite not filtering on f1samp
		#post2_alt <- abc(target=targ[2:3], param=pams, sumstat=stats[,2:3], tol=0.001, method='loclinear', transf=c('log', 'logit', 'none'), logit.bounds=matrix(c(0,0,0,1,1,1),ncol=2)) # fails with colinearity
		#post3_alt <- abc(target=targ[2:3], param=pams, sumstat=stats[,2:3], tol=0.01, method='ridge', transf=c('log', 'logit', 'none'), logit.bounds=matrix(c(0,0,0,1,1,1),ncol=2)) # fails with infinite values
		#post4_alt <- abc(target=targ[2:3], param=pams, sumstat=stats[,2:3], tol=0.001, method='neuralnet', transf=c('log', 'logit', 'none'), logit.bounds=matrix(c(0,0,0,1,1,1),ncol=2)) # slow
	
#		sums <- summary(post, print=FALSE)
#		ret <- c('ne_mean'=sums[4,1], 'ne_l95'=sums[2,1], 'ne_u95'=sums[6,1], 'f1_mean'=sums[4,2], 'f1_l95'=sums[2,2], 'f1_u95'=sums[6,2], 's_mean'=sums[4,3], 's_l95'=sums[2,3], 's_u95'=sums[6,3])
		ret <- matrix(data=c(post$ss[,1], post$ss[,2], post$ss[,3], post$unadj.values[,1], post$unadj.values[,2], post$unadj.values[,3]), nrow=6, byrow=TRUE)
		rownames(ret) <- c('f1samp', 'fsdprime', 'fsiprime', 'ne', 'f1', 's')

		return(ret)
	}
}

# start cluster
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	cl <- makeCluster(detectCores()-1) # set up cluster on my mac (or another computer), using all but one core
	clusterEvalQ(cl, library(data.table))
	clusterEvalQ(cl, library(abc.data))
	clusterEvalQ(cl, library(locfit))
	clusterEvalQ(cl, library(abc))
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	cl <- makeCluster(20) # set up 20-core cluster on a cod node
	clusterEvalQ(cl, library(data.table, lib.loc='/projects/cees/lib/R_packages/'))
	clusterEvalQ(cl, library(abc.data, lib.loc='/projects/cees/lib/R_packages/'))
	clusterEvalQ(cl, library(locfit, lib.loc='/projects/cees/lib/R_packages/'))
	clusterEvalQ(cl, library(abc, lib.loc='/projects/cees/lib/R_packages/'))
}

clusterExport(cl, c('runabc'))


# run calculations in parallel
temp <- parApply(cl, cbind(sampsize, alcnts, obs), MARGIN=1, FUN=runabc, out=out) # a slow start: have to load out onto each node. output is annoying: a 60000 x nrow(obs) matrix. Rows are f1samp, fsdprime, fsiprime, ne, f1, s, f1samp, fsdprime, fsiprime, ne, f1, s, etc. (samples from the posterior). Columns are each locus.

save(temp, file='analysis/temp/temp.rdata')

# reformat
posts <- list(f1samp=t(temp[seq(1,nrow(temp),by=6),]), fsdprime=t(temp[seq(2,nrow(temp),by=6),]), fsiprime=t(temp[seq(3,nrow(temp),by=6),]), ne=t(temp[seq(4,nrow(temp),by=6),]), f1 = t(temp[seq(5,nrow(temp),by=6),]), s = t(temp[seq(6,nrow(temp),by=6),]))

# save
save(posts, file='analysis/wfs_abc_posts2011-2014.rdata')

# stop cluster
stopCluster(cl)
