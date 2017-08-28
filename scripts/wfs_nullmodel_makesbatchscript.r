# write out a shell script that will submit an sbatch job for every combination of sample sizes that we need

if(!grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table)
}
if(grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}

# settings
#pop <- 'Lof'; myyr1 <- '07'; myyr2 <- '11'; kmer <- 25
pop <- 'Lof'; myyr1 <- '07'; myyr2 <- '11'; kmer <- 150
#pop <- 'Lof'; myyr1 <- '11'; myyr2 <- '14'; kmer <- 25
#pop <- 'Can'; myyr1 <- '00'; myyr2 <- '00'; kmer <- 25 # myyr are placeholders since only one set of years for Canada
#pop <- 'Can'; myyr1 <- '00'; myyr2 <- '00'; kmer <- 150 # myyr are placeholders since only one set of years for Canada


# load observed data
if(pop == 'Lof'){
	targfile <- paste('data_29.06.17/Frequency_table_Lof', myyr1, '_Lof', myyr2, '_', kmer, 'k.txt', sep='')
}
if(pop == 'Can'){
	targfile <- paste('data_11.07.17/Frequency_table_Can_40_Can_', kmer, 'k.txt', sep='')
}
targ <- fread(targfile, header=TRUE)
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locusnum:=1:nrow(targ)] # add a locus number indicator


# how many loci at each sample size
setkey(targ, alcnt1, alcnt2)
nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2)]
	nrow(nloci) # 1907-2014: 176 (25kmer), 193 (150k) (how many sample sizes)
				# 1907-2011: 249 (25k), 267 (150k)
				# 2011-2014: 162 (25k)
				# Can: 		55 (25k), 55 (150k)

# sample sizes with >100 loci
nloci[nloci>5000,]

# check which sample sizes are already run
nloci[,todo:=1] # set up a flag for which we should run
nloci[,nlocitorun:=nloci] # number of loci left to run at each sample size

for(i in 1:nrow(nloci)){
	# check if the simulation files for this sample size exists
	myalcnt1 <- nloci[i,alcnt1]
	myalcnt2 <- nloci[i,alcnt2]
	
	exist <- file.exists(paste('analysis/temp/wfs_simsnull_ff', myalcnt1, ',', myalcnt2, '.RData', sep='')) # null model sims

	if(exist){
		# find the chunk of loci of the specified sample size to operate on
		samplepartinds <- targ[alcnt1==myalcnt1 & alcnt2==myalcnt2,locusnum] # locus numbers in each sample size chunk
		orig.n <- length(samplepartinds)

		# check for existing files (don't re-run these)
		existingfiles <- list.files(path='analysis/temp', pattern=paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''))
	#	print(paste('found', length(existingfiles), 'existing file(s) relevant to this sample size'))
		existingrngs <- strsplit(gsub(paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus|.csv.gz', sep=''), '', existingfiles), split='-') # extract just the locus ranges
		existingrngs <- lapply(existingrngs, as.numeric)
		existingloci <- numeric(0)
		if(length(existingrngs)>0) for(j in 1:length(existingrngs)) existingloci <- c(existingloci, existingrngs[[j]][1]:existingrngs[[j]][2])

		# find the set of loci remaining to run
		samplepartinds <- setdiff(samplepartinds, existingloci) # remove loci already run (from existing files)
		sampleparts.n <- length(samplepartinds) # number of loci in the samplesize chunk

		# flag which ones still need to be run
		if(sampleparts.n==0){
			nloci[i,todo:=0]
			nloci[i,nlocitorun:=0]
		} else {
			print(paste(myalcnt1, myalcnt2, 'Run', sampleparts.n, 'loci of', orig.n, 'loci originally'))
			nloci[i,nlocitorun:=sampleparts.n]
		}
	} else {
		nloci[i,todo:=0] # can't run this now because sim files don't exist
		print(paste(myalcnt1, myalcnt2, 'Cannot run because not all sim files exist'))
		warning(paste(myalcnt1, myalcnt2, 'Cannot run because not all sim files exist'))
	}
}
	
# how many sample sizes to run?
sum(nloci$todo)
nrow(nloci) # compare to number we could run
	
# any sample sizes with < all loci to run?
nloci[nlocitorun<nloci,]
	
# write out the shell script
outfile <- 'scripts/wfs_nullmodel_submit_all_sbatch1.sh'

cat('#!/bin/bash\n', file=outfile, append=FALSE) # header

for(i in which(nloci[,todo==1])){
	cat(paste('sbatch --job-name=nlmd', nloci[i,alcnt1], ',', nloci[i,alcnt2], ' --cpus-per-task=', min(16, nloci[i,nloci]), ' scripts/wfs_nullmodel_sbatch.sh ', nloci[i,alcnt1], ' ', nloci[i, alcnt2], ' ', pop, ' ', myyr1, ' ', myyr2, ' ', kmer, ' ', min(16, nloci[i,nloci]), '\n', sep=''), file=outfile, append=TRUE)
}
