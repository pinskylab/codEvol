# Make simulations for the null model
# best to run on cod node with nohup if many sample sizes to make
# Arguments (ONLY SOME ARE IMPLEMENTED RIGHT NOW, SEE wfs_make_sims_null_3times.r FOR EXAMPLE)
#	pop: Can or Lof
#	yr1: 07 or 11 (not needed for Can)
#	yr2: 11 or 14 (not needed for Can)
#	half: 1, 2, or 3 (which half of loci to run, or 3 for run all)
# 	ncores: number of cores to use
#	trimlowsampsize: 1 to trim loci with less than half of individuals genotypes, 0 to not trim
#	rerunlow: 0 for default, 1 to only run simulations for loci with low p-values in wfs_nullmodel_pos&pvals_07-11-14.rds, 2 to only run even more simulations for loci with low p-values in wfs_nullmodel_outliers_lowp_Lof_07-11-14.tsv.gz

# read arguments
args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    pop <- 'Lof'
    yr1 <- '07'
    yr2 <- '14'
	half <- 3 # supply default values
	ncores <- 30
	trimlowsampsize <- 1
	rerunlow <- 0
}else{
	for(i in 1:length(args)){
		eval(parse(text=args[[i]]))
	}
	if(!(pop %in% c('Lof', 'Can'))) stop('pop must be Lof or Can!')
	if(pop=='Lof' & !(yr1 %in% c('07', '11'))) stop('yr1 must be 07 or 11 for Lof!')
	if(pop=='Lof' & !(yr2 %in% c('11', '14'))) stop('yr2 must be 11 or 14 for Lof!')
	if(!(half %in% 1:3)) stop('half must be 1, 2, or 3!')
	if(!(ncores > 0)) stop('ncores must be positive')
	if(!(trimlowsampsize %in% 0:1)) stop('trimlowsampsize must be 0 or 1!')
	if(!(rerunlow %in% 0:2)) stop('rerunlow must be 0 or 1 or 2!')
}
print(paste('Arguments: pop=', pop, ', yr1=', yr1, ', yr2=', yr2, ', half=', half, ', ncores=', ncores, ', trimlowsampsize=', trimlowsampsize, ', rerunlow=', rerunlow, sep=''))


# load functions
source('scripts/wfs_byf1samp.r')
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


# read in Ne and frequency data
if(pop=='Lof'){
	nes <- read.table('analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	freqfile <- paste('data_2018.09.05/Frequency_table_Lof', yr1, '_Lof', yr2, '.txt', sep='')
	nchrs <- fread(freqfile, header=TRUE) # read in frequency table data
	gen <- 11
}
if(pop=='Can'){
	nes <- read.table('analysis/Can_40_to_Can.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	freqfile <- paste('data_2018.09.05/Frequency_table_CAN_40_TGA.txt', sep='')
	nchrs <- fread(freqfile, header=TRUE)
	gen <- 8
}

setnames(nchrs, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))

# trim out missing loci
	nrow(nchrs)
nchrs <- nchrs[N_CHR_1>0 & N_CHR_2>0,]
	print(nrow(nchrs))

# trim to loci with at least half of individuals genotyped?
if(trimlowsampsize==1){
	print('Trimming to loci with at least half of the individuals genotyped')
	nrow(nchrs)
	nchrs <- nchrs[N_CHR_1>=max(N_CHR_1)/2 & N_CHR_2>=max(N_CHR_2)/2,]
	print(nrow(nchrs))
}	

# trim only to loci with low p-values in wfs_nullmodel_pos&pvals_*.rds (from analysis of a previous run)?
if(rerunlow==1){
	print('Trimming to loci with p<=8e-6 (4 out of 500,000)')
	if(pop=='Lof'){
		pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-11-14.rds'))
	}
	if(pop=='Can'){
		pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_Can.rds')) # NOT SURE YET THIS WORKS FOR CAN
	}
	nchrs <- merge(nchrs, pvals[,.(CHROM, POS, p)]) # merge in p-values
	print(nrow(nchrs))		
	nchrs <- nchrs[p<=8e-6,]
	print(nrow(nchrs))
}

# trim only to loci with low p-values in wfs_nullmodel_outliers_lowp_*.tsv.gz (from analysis of a previous rerunlow=1 run)?
if(rerunlow==2){
	print('Trimming to loci with p<=1e-6 (4 out of 5,000,000)')
	
	if(pop=='Lof'){
		pvals <- fread("gunzip -c analysis/wfs_nullmodel_outliers_lowp_Lof_07-11-14.tsv.gz")
		nchrs <- merge(nchrs, pvals[,.(CHROM, POS, pLof071114low)], all.x=TRUE) # merge in p-values for Lof
		print(nrow(nchrs))		
		nchrs <- nchrs[pLof071114low<=1e-6,]
		print(nrow(nchrs))
	}
	if(pop=='Can'){
		pvals <- fread("gunzip -c analysis/wfs_nullmodel_outliers_lowp_Can.tsv.gz") # NOT SURE YET THIS WORKS FOR CAN
		nchrs <- merge(nchrs, pvals[,.(CHROM, POS, pCanlow)], all.x=TRUE) # merge in p-values
		print(nrow(nchrs))		
		nchrs <- nchrs[pCanlow<=1e-6,]
		print(nrow(nchrs))
	}
}


# parameters for the simulation
nsims <- 500000 # number of sims to run for each starting frequency in each sample size
nreps <- 1 # number of repetitions of the number of simulations
if(rerunlow==1) nreps <- 10 # number of reps to run for loci with low p-values (50M sims total)
if(rerunlow==2) nreps <- 100 # number of reps to run for loci with really low p-values (500M sims total)
print(paste('Running', nreps, 'reps of', nsims, 'simulations for', nreps*nsims, 'total per sample size.'))

c1s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_1] # the sample sizes to simulate (first sample)
c2s <- nchrs[!duplicated(paste(N_CHR_1, N_CHR_2)),N_CHR_2]
# 	print(c1s)
# 	print(c2s)
	length(c1s)
	length(c2s)

# expand to number of reps (if nreps>1)
print('Expanding sample sizes to number of repetitions')
reps <- rep(1:nreps, times=length(c1s)) # holds repetition counter
c1s <- rep(c1s, each=nreps)
c2s <- rep(c2s, each=nreps)
	print(length(c1s))
	print(length(c2s))
	print(length(reps))
	

# trim to first or second half (or keep all)
if(half==1){
	print('Using first half of sample sizes')
	c1s <- c1s[1:floor(length(c1s)/2)]
	c2s <- c2s[1:floor(length(c2s)/2)]
	reps <- reps[1:floor(length(reps)/2)] # also for reps
}
if(half==2){
	print('Using second half of sample sizes')
	c1s <- c1s[(floor(length(c1s)/2)+1):length(c1s)]
	c2s <- c2s[(floor(length(c2s)/2)+1):length(c2s)]
	reps <- reps[(floor(length(reps)/2)+1):length(reps)]
}
if(half==3){
	print('Using all sample sizes')
}
	print(length(c1s))
	print(length(c2s))
	print(length(reps))


# check that nes >0 
print(summary(nes))
print(paste('Of', length(nes), 'Ne values', sum(nes > 0), 'are >0. Trimming all others out.'))
nes <- nes[nes>0]

# check temp directory for previous simulations (if Ne, gen and other parameters are the same)
if(pop=='Lof' & rerunlow==0){
	torun <- paste(c1s, c2s, sep=',') # the list of sample sizes we need to run. these have all starting frequencies
	existing <- list.files(path='analysis/temp', pattern='wfs_simsnull_ff.*\\.ffData') # the existing simulations. note that this search string assumes the simulations were made by a previous version of the script that didn't append rep# to the end of the file name (see rerunlow==1 below). Changed March 5, 2019.
	existing <- gsub('wfs_simsnull_ff|.ffData', '', existing) # trim file names to just sample sizes

	print('Trimming out files already run for 1 nrep')
	keep <- !(torun %in% existing) # the sample sizes that still need to be run
	c1s <- c1s[keep]
	c2s <- c2s[keep]
	reps <- reps[keep]
}
if(pop=='Lof' & rerunlow %in% 1:2){ # if we need more than 1 rep of each sample size. these files only have the starting frequencies we need
	torun <- paste(paste(c1s, c2s, sep=','), reps, sep='_') # the list of sample sizes we need to run
	existing <- list.files(path='analysis/temp', pattern='wfs_simsnull_ff.+_[[:digit:]]+\\.ffData') # the existing simulations
	existing <- gsub('wfs_simsnull_ff|.ffData', '', existing) # trim file names to just sample sizes and rep#

	print('Trimming out files already run for nrep>1')
	keep <- !(torun %in% existing) # the sample sizes that still need to be run
	c1s <- c1s[keep]
	c2s <- c2s[keep]
	reps <- reps[keep]
}
if(pop=='Can'){
	torun <- paste(paste(c1s, c2s, sep=','), reps, sep='_') # the list of sample sizes we need to run
	existing <- list.files(path='analysis/temp', pattern='wfs_simsnullCAN_ff.+_[[:digit:]]+\\.ffData') # the existing simulations. works on new and old simulations because I renamed all older sim files to have _1 as of March 5 2019.
	existing <- gsub('wfs_simsnullCAN_ff|.ffData', '', existing) # trim file names to just sample sizes and rep#

	print('Trimming out files already run for nrep>1')
	keep <- !(torun %in% existing) # the sample sizes that still need to be run
	c1s <- c1s[keep]
	c2s <- c2s[keep]
	reps <- reps[keep]
}

print(length(c1s)) # how many to run?
print(length(c2s))
print(length(reps))

# run simulations for missing sample sizes
if(length(c1s)>0){

	# start cluster
	if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
		cl <- makeCluster(detectCores()-1) # set up cluster on my mac (or another computer), using all but one core
	}
	if(grepl('hpc.uio.no', Sys.info()["nodename"])){
		cl <- makeCluster(ncores) # set up cluster on a cod node
	}
	clusterExport(cl, c('wfs_byf1samp'))


	# check where this is to monitor temp file creation in the file system /tmp/Rtmpp__
	options('fftempdir')


	# make simulations (parallel way)
	for(i in 1:length(c1s)){ # loop through each set of sample sizes
		print(paste('Sample size', i, 'of', length(c1s), 'to do at', Sys.time()))
		
		# only run observed starting sample frequencies
		nchr1s <- nchrs[c1s[i]==N_CHR_1 & c2s[i]==N_CHR_2, unique(round(Freq_1*N_CHR_1))] # number of starting copies of allele for the retained loci
		f1samps=rep((nchr1s)/c1s[i], rep(nsims, length(nchr1s)))	
		
		thisout <- parSapply(cl, f1samps, FUN=wfs_byf1samp, smin=0, smax=0, c1=c1s[i], c2=c2s[i], gen=gen, ne=nes, h=0.5, simplify=TRUE)

		thisout.ff <- ff(thisout, dim=dim(thisout), dimnames=dimnames(thisout)) # create in tempdir
	
		# save to permanent file (semi-permanent.. in my temp directory)
		if(pop=='Lof'){
			ffsave(thisout.ff, file=paste('analysis/temp/wfs_simsnull_ff', paste(c1s[i], c2s[i], sep=','), '_', reps[i], sep=''))

		}
		if(pop=='Can'){
			ffsave(thisout.ff, file=paste('analysis/temp/wfs_simsnullCAN_ff', paste(c1s[i], c2s[i], sep=','), '_', reps[i], sep=''))
		}
		
		# remove the temp files
		delete(thisout.ff)
		rm(thisout.ff)

	}


	# stop the cluster
	stopCluster(cl)
}