# load functions
source('scripts/wfs.r')
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel)
	require(data.table)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
	require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# data
nes <- read.table('analysis/Can_40_to_Can_150k.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	print(summary(nes)) # to check
nchrs <- fread('data_11.07.17/Frequency_table_Can_40_Can_25k.txt', header=TRUE) # to figure out the sample sizes to simulate
setnames(nchrs, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))

# parameters
nsims <- 10000000
c1s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_1] # the sample sizes to simulate (first sample)
c2s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_2]
	print(c1s)
	print(c2s)
gen <- 8

# check that nes >0 
print(summary(nes))
print(paste('Of', length(nes), 'Ne values', sum(nes > 0), 'are >0. Trimming all others out.'))
nes <- nes[nes>0]
	length(nes)


# check temp directory for previous simulations (if Ne, gen and other parameters are the same)
torun <- paste(c1s, c2s, sep=',') # the list of sample sizes I need to run
existing <- list.files(path='analysis/temp', pattern='wfs_simsnullCAN_ff.*\\.ffData') # the existing simulations
existing <- gsub('wfs_simsnullCAN_ff|.ffData', '', existing) # trim file names to just sample sizes
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
		cl <- makeCluster(30) # set up 30-core cluster on a cod node
	}
	clusterExport(cl, c('wfs'))


	# check where this is to monitor temp file creation in the file system /tmp/Rtmpp__
	options('fftempdir')


	# make simulations (parallel way)
	for(i in 1:length(c1s)){
		print(paste('Sample size', i, 'of', length(c1s), 'to do.'))

		thisout <- parSapply(cl, 1:nsims, FUN=wfs, f1min=0, f1max=1, smin=0, smax=0, c1=c1s[i], c2=c2s[i], gen=gen, ne=nes, h=0.5, simplify=TRUE)
			
		thisout.ff <- ff(thisout, dim=dim(thisout), dimnames=dimnames(thisout)) # create in tempdir
	
		# save to permanent file (semi-permanent.. but in my temp directory)
		ffsave(thisout.ff, file=paste('analysis/temp/wfs_simsnullCAN_ff', paste(c1s[i], c2s[i], sep=','), sep=''))

		# remove the temp files
		delete(thisout.ff)
		rm(thisout.ff)

	}

	# stop the cluster
	stopCluster(cl)

}


