# Make simulations for the ABC model for CAN samples
# best to run on cod node with nohup if many sample sizes to make

# parameters
nsims <- 500000 # number of sims to run for each starting frequency in each sample size
gen <- 8 # number of generations between samples
nefile <- 'analysis/Can_40_to_Can.w_Ne_bootstrap.txt'
freqfile <- 'data_2017.11.24/Frequency_table_CAN_40_TGA.txt'
ncores <- 30

# load functions
source('scripts/wfs_byf1samp.r')
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel)
	require(data.table)
	require(bit)
	require(ff)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
	require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# read in Ne data
nes <- read.table(nefile)[,1] # the values of Ne from wfabc_1
nchrs <- fread(freqfile, header=TRUE)
setnames(nchrs, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))

# trim to focal loci (outliers in CAN)
locs <- fread('analysis/outlier_annotation.csv') # the outliers
	locs <- locs[q3.Can !='' | q3.comb071114Can !='',.(CHROM, POS)]
	locs[, POS := as.numeric(POS)]
nrow(locs)
setkey(locs, CHROM, POS)
setkey(nchrs, CHROM, POS)
nchrs <- merge(nchrs, locs)
nrow(nchrs)


# find sample sizes to run
c1s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_1] # the sample sizes to simulate (first sample)
c2s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_2]
	print(c1s)
	print(c2s)

# check that nes >0 
print(summary(nes))
print(paste('Of', length(nes), 'Ne values', sum(nes > 0), 'are >0. Trimming all others out.'))
nes <- nes[nes>0]

# check temp directory for previous simulations (if Ne, gen and other parameters are the same)
torun <- paste(c1s, c2s, sep=',') # the list of sample sizes I need to run
existing <- list.files(path='analysis/temp', pattern='wfs_simsCAN_ff.*\\.ffData') # the existing simulations
existing <- gsub('wfs_simsCAN_ff|.ffData', '', existing) # trim file names to just sample sizes

keep <- !(torun %in% existing) # the sample sizes that still need to be run
c1s <- c1s[keep]
c2s <- c2s[keep]
length(c1s) # how many to run?

# run simulations for missing sample sizes
if(length(c1s)>0){

	# start cluster
	if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
		cl <- makeCluster(detectCores()-1) # set up cluster on my mac (or another computer), using all but one core
	}
	if(grepl('hpc.uio.no', Sys.info()["nodename"])){
		cl <- makeCluster(ncores) # set up 30-core cluster on a cod node
	}
	clusterExport(cl, c('wfs_byf1samp'))


	# check where this is to monitor temp file creation in the file system /tmp/Rtmpp__
	options('fftempdir')


	# make simulations (parallel way)
	for(i in 1:length(c1s)){ # loop through each set of sample sizes
		print(paste('Sample size', i, 'of', length(c1s), 'to do.'))
		
		# make list of all possible starting sample frequencies and make nsims copies of each
		f1samps=rep((0:c1s[i])/c1s[i], rep(nsims, c1s[i]+1))

		thisout <- parSapply(cl, f1samps, FUN=wfs_byf1samp, smin=-1, smax=1, c1=c1s[i], c2=c2s[i], gen=gen, ne=nes, h=0.5, simplify=TRUE)

		thisout.ff <- ff(thisout, dim=dim(thisout), dimnames=dimnames(thisout)) # create in tempdir
	
		# save to permanent file (semi-permanent.. but in my temp directory)
		ffsave(thisout.ff, file=paste('analysis/temp/wfs_simsCAN_ff', paste(c1s[i], c2s[i], sep=','), sep=''))
		
		# remove the temp files
		delete(thisout.ff)
		rm(thisout.ff)

	}


	# stop the cluster
	stopCluster(cl)
}