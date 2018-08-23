# Make simulations to estimate s with ABC
# Use data from 3 samples
# Use only outlier loci
# best to run on cod node with nohup if many sample sizes to make

# set parameters
pop <- 'Lof'; yrs<-'07-11-14'
ncores <- 30 # cod nodexw

# load functions
source('scripts/wfs_byf1samp_3samps.r')
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


# read in Ne data
if(pop=='Lof'){
	nes <- read.table('analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	dat14 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat14, 3:7, c('N_CHR_1', 'Freq_07', 'N_CHR_3', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
	dat11 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(dat11, 3:7, c('N_CHR_1', 'Freq_07', 'N_CHR_2', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
	setkey(dat11, CHROM, POS)
	setkey(dat14, CHROM, POS)

	nchrs <- dat11[dat14, .(CHROM, POS, N_CHR_1, N_CHR_2, N_CHR_3)]
	gen1 <- 11
	gen2 <- 11
}

# trim to focal loci (outliers in Lof)
locs <- fread('analysis/outlier_annotation.csv') # the outliers
	locs <- locs[q3.Lof071114 !='' | q3.comb071114Can !='',.(CHROM, POS)]
	locs[, POS := as.numeric(POS)]
nrow(locs)
setkey(locs, CHROM, POS)
setkey(nchrs, CHROM, POS)
nchrs <- merge(nchrs, locs)
nrow(nchrs)

# parameters
nsims <- 500000 # number of sims to run for each starting frequency in each sample size
c1s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2, N_CHR_3)),N_CHR_1] # the sample sizes to simulate (first sample)
c2s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2, N_CHR_3)),N_CHR_2]
c3s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2, N_CHR_3)),N_CHR_3]
	print(c1s)
	print(c2s)
	print(c3s)

# check that nes >0 
print(summary(nes))
print(paste('Of', length(nes), 'Ne values', sum(nes > 0), 'are >0. Trimming all others out.'))
nes <- nes[nes>0]

# check temp directory for previous simulations (if Ne, gen and other parameters are the same)
torun <- paste(c1s, c2s, c3s, sep=',') # the list of sample sizes I need to run
if(pop=='Lof'){
	existing <- list.files(path='analysis/temp', pattern='wfs_sims_ff.*\\.ffData') # the existing simulations
	existing <- gsub('wfs_sims_ff|.ffData', '', existing) # trim file names to just sample sizes
}

keep <- !(torun %in% existing) # the sample sizes that still need to be run
c1s <- c1s[keep]
c2s <- c2s[keep]
c3s <- c3s[keep]
length(c1s) # how many to run?


# run simulations for missing sample sizes
if(length(c1s)>0){

	# start cluster
	if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
		cl <- makeCluster(detectCores()-1) # set up cluster on my mac (or another computer), using all but one core
	}
	if(grepl('hpc.uio.no', Sys.info()["nodename"])){
		cl <- makeCluster(ncores) # set up cluster on a cod node
	}
	clusterExport(cl, c('wfs_byf1samp_3samps'))


	# check where this is to monitor temp file creation in the file system /tmp/Rtmpp__
	options('fftempdir')

	# make simulations (parallel way)
	for(i in 1:length(c1s)){ # loop through each set of sample sizes
		print(paste('Sample size', i, 'of', length(c1s), 'to do.'))
		
		# make list of all possible starting sample frequencies and make nsims copies of each
		f1samps=rep((0:c1s[i])/c1s[i], rep(nsims, c1s[i]+1))
		
		thisout <- parSapply(cl, f1samps, FUN=wfs_byf1samp_3samps, smin=-1, smax=1, c1=c1s[i], c2=c2s[i], c3=c3s[i], gen1=gen1, gen2=gen2, ne=nes, h=0.5, simplify=TRUE)

		thisout.ff <- ff(thisout, dim=dim(thisout), dimnames=dimnames(thisout)) # create in tempdir
	
		# save to permanent file (semi-permanent.. but in my temp directory)
		if(pop=='Lof'){
			ffsave(thisout.ff, file=paste('analysis/temp/wfs_sims_ff', paste(c1s[i], c2s[i], c3s[i], sep=','), sep=''))
		}
		
		# remove the temp files
		delete(thisout.ff)
		rm(thisout.ff)
	}

	# stop the cluster
	stopCluster(cl)
}
