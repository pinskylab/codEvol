# write out a shell script that will submit an sbatch job for every combination of sample sizes that we need

if(!grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table)
}
if(grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# load observed data
targ <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE)
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locusnum:=1:nrow(targ)] # add a locus number indicator


# how many loci at each sample size
setkey(targ, alcnt1, alcnt2)
nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2)]
	nrow(nloci) # 176

# sample sizes with >100 loci
nloci[nloci>5000,]

# check which sample sizes are already run
nloci[,todo:=1] # set up a flag for which we should run

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
		} else {
			print(paste(myalcnt1, myalcnt2, 'Run', sampleparts.n, 'loci of', orig.n, 'loci originally'))
			nloci[i,nloci:=sampleparts.n]
		}
	} else {
		nloci[i,todo:=0] # can't run this now because sim files don't exist
		print(paste(myalcnt1, myalcnt2, 'Cannot run because not all sim files exist'))
		warning(paste(myalcnt1, myalcnt2, 'Cannot run because not all sim files exist'))
	}
}
	
# how many files to run?
sum(nloci$todo)
	
# write out the shell script
outfile <- 'scripts/wfs_nullmodel_submit_all_sbatch1.sh'

cat('#!/bin/bash\n', file=outfile, append=FALSE) # header

for(i in which(nloci[,todo==1])){
	cat(paste('sbatch --job-name=nlmd', nloci[i,alcnt1], ',', nloci[i,alcnt2], ' --cpus-per-task=', min(16, nloci[i,nloci]), ' scripts/wfs_nullmodel_sbatch.sh ', nloci[i,alcnt1], ' ', nloci[i, alcnt2], ' ', min(16, nloci[i,nloci]), '\n', sep=''), file=outfile, append=TRUE)
}