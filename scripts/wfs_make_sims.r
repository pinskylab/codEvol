# load functions
source('scripts/wfs.r')
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel)
	require(data.table)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# data
nes <- read.table('analysis/LOF_07_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	print(summary(nes)) # to check
nchrs <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE) # to figure out the sample sizes to simulate
setnames(nchrs, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))

# parameters
nsims <- 100000
c1s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_1] # the sample sizes to simulate (first sample)
c2s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_2]
	print(c1s)
	print(c2s)

# start cluster
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	cl <- makeCluster(detectCores()-1) # set up cluster on my mac (or another computer), using all but one core
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	cl <- makeCluster(30) # set up 30-core cluster on a cod node
}
clusterExport(cl, c('sumstats', 'wfs'))

# make simulations (parallel way)
	out <- array(dim=c(9,nsims,length(c1s)), dimnames=list(var=c('ne', 'f1', 's', 'gen', 'f2', 'f1samp', 'f2samp', 'fsdprime', 'fsiprime'), sim=1:nsims, sampsize=paste(c1s, c2s, sep=',')))

	for(i in 1:length(c1s)){
		out[,,i] <- parSapply(cl, 1:nsims, FUN=wfs, f1min=0, f1max=1, smin=-1, smax=1, c1=c1s[i], c2=c2s[i], gen=11, ne=nes, h=0.5, simplify=TRUE)
	}

# make simulations (data.table way: failed when reach 200G RAM)
#a <- rep(as.numeric(NA), nsims*length(c1s))
#out <- data.table(i=1:length(a), c1=rep(c1s,rep(nsims,length(c1s))), c2=rep(c2s,rep(nsims,length(c2s))), ne=a, f1=a, s=a, gen=a, f2=a, f1samp=a, f2samp=a, fsdprime=a, fsiprime=a) # takes 154G RAM with nsims=10M
#
#out[,c('ne', 'f1', 's', 'gen', 'f2', 'f1samp', 'f2samp', 'fsdprime', 'fsiprime') := wfs(i=1,f1min=0, f1max=1, smin=-1, smax=1, c1=c1, c2=c2, gen=11, ne=nes, h=0.5), by=i]


# write out
save(out, file='analysis/wfs_sims.rdata')

# stop the cluster
stopCluster(cl)

