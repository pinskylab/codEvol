# write out a shell script that will submit an sbatch job for every combination of sample sizes that we need

if(!grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table)
}
if(grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}

# settings
#pop <- 'Lof'; myyr1 <- '07'; myyr2 <- '14'
#pop <- 'Lof'; myyr1 <- '07'; myyr2 <- '11'
#pop <- 'Lof'; myyr1 <- '11'; myyr2 <- '14'
#pop <- 'Lof'; myyr1 <- '07'; myyr2 <- '1114' # all 3 time points
pop <- 'Can'; myyr1 <- '00'; myyr2 <- '00' # myyr are placeholders since only one set of years for Canada


# load observed data
if(pop == 'Lof'){
	if(myyr2 != '1114'){
		targfile <- paste('data_2018.09.05/Frequency_table_Lof', myyr1, '_Lof', myyr2, '.txt', sep='')
	}
	if(myyr2 == '1114'){
		targfile <- paste('data_2018.09.05/Frequency_table_Lof07_Lof11.txt', sep='')
		targfile2 <- paste('data_2018.09.05/Frequency_table_Lof07_Lof14.txt', sep='')
	}
}
if(pop == 'Can'){
	targfile <- paste('data_2018.09.05/Frequency_table_CAN_40_TGA.txt', sep='')
}
targ <- fread(targfile, header=TRUE)
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locusnum:=1:nrow(targ)] # add a locus number indicator

if(myyr2=='1114'){
	targ2 <- fread(targfile2, header=TRUE)
	setnames(targ2, 3:7, c('alcnt1', 'f1samp', 'alcnt3', 'f3samp', 'ABS_DIFF13'))
	setkey(targ, CHROM, POS)
	setkey(targ2, CHROM, POS)
	targ <- targ[targ2,.(locusnum, CHROM, POS, alcnt1, alcnt2, alcnt3, f1samp, f2samp, f3samp)]
}

# trim out missing loci
	nrow(targ)
if(myyr2!='1114') targ <- targ[alcnt1>0 & alcnt2>0,]
if(myyr2=='1114') targ <- targ[alcnt1>0 & alcnt2>0 & alcnt3>0,]
	print(nrow(targ))

# trim to loci with at least half of individuals genotyped
# 	nrow(targ)
# if(myyr2!='1114') targ <- targ[alcnt1>=max(alcnt1)/2 & alcnt2>=max(alcnt2)/2,]
# if(myyr2=='1114') targ <- targ[alcnt1>=max(alcnt1)/2 & alcnt2>=max(alcnt2)/2 & alcnt3>=max(alcnt3)/2,]
# 	print(nrow(targ))


# how many loci at each sample size
if(myyr2 != '1114'){
	setkey(targ, alcnt1, alcnt2)
	nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2)]
}
if(myyr2 == '1114'){
	setkey(targ, alcnt1, alcnt2, alcnt3)
	nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2, alcnt3)]
}
	nrow(nloci) # Lof 1907-2011-2014: 1717 (>=50% indivs genotyped) 6152 (all)
				# Lof 1907-2014: 
				# Can: 143 (>=50% genotyped) 497 (all)

# sample sizes with >5000 loci
nloci[nloci>5000, .(.N, max(nloci))]  # Lof 1907-2011-2014: 0 sample sizes
					# Lof 1907-2014:  samplesizes					
					# Can: 29 (>=50%) 49 (all)

# check which sample sizes are already run
nloci[,todo:=1] # set up a flag for which we should run
nloci[,nlocitorun:=nloci] # number of loci left to run at each sample size

for(i in 1:nrow(nloci)){
	cat(paste(i, " ", sep=''))
	# check if the simulation files for this sample size exists
	myalcnt1 <- nloci[i,alcnt1]
	myalcnt2 <- nloci[i,alcnt2]
	if(myyr2=='1114') myalcnt3 <- nloci[i,alcnt3]
	
	if(pop=='Lof'){
		exist <- file.exists(paste('analysis/temp/wfs_simsnull_ff', myalcnt1, ',', myalcnt2, '_1.RData', sep='')) # null model sims
		if(myyr2=='1114') exist <- file.exists(paste('analysis/temp/wfs_simsnull_ff', myalcnt1, ',', myalcnt2, ',', myalcnt3, '_1.RData', sep='')) # null model sims
	}
	if(pop=='Can'){
		exist <- file.exists(paste('analysis/temp/wfs_simsnullCAN_ff', myalcnt1, ',', myalcnt2, '_1.RData', sep='')) # null model sims	
	}

	if(exist){
		# find the chunk of loci of the specified sample size to operate on
		samplepartinds <- targ[alcnt1==myalcnt1 & alcnt2==myalcnt2,locusnum] # locus numbers in each sample size chunk
		if(myyr2=='1114') samplepartinds <- targ[alcnt1==myalcnt1 & alcnt2==myalcnt2 & alcnt3==myalcnt3,locusnum]
		orig.n <- length(samplepartinds)

		# check for existing files (don't re-run these)
		if(myyr2 != '1114'){
			existingfiles <- list.files(path='analysis/temp', pattern=paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''))
		#	print(paste('found', length(existingfiles), 'existing file(s) relevant to this sample size'))
			existingrngs <- strsplit(gsub(paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus|.csv.gz', sep=''), '', existingfiles), split='-') # extract just the locus ranges
		}
		
		if(myyr2=='1114'){ # 3 time points
			existingfiles <- list.files(path='analysis/temp', pattern=paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, myalcnt3, sep=','), '_locus*', sep=''))
		#	print(paste('found', length(existingfiles), 'existing file(s) relevant to this sample size'))
			existingrngs <- strsplit(gsub(paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, myalcnt3, sep=','), '_locus|.csv.gz', sep=''), '', existingfiles), split='-') # extract just the locus ranges
		}

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
			if(myyr2 != '1114') print(paste(i, 'of', nrow(nloci), ':', myalcnt1, myalcnt2, 'Run', sampleparts.n, 'loci of', orig.n, 'loci originally'))
			if(myyr2 == '1114') print(paste(i, 'of', nrow(nloci), ':', myalcnt1, myalcnt2, myalcnt3, 'Run', sampleparts.n, 'loci of', orig.n, 'loci originally'))
			nloci[i,nlocitorun:=sampleparts.n]
		}
	} else {
		nloci[i,todo:=0] # can't run this now because sim files don't exist
		if(myyr2!='1114'){
			print(paste(i, 'of', nrow(nloci), ':', myalcnt1, myalcnt2, 'Cannot run because not all sim files exist'))
			warning(paste(i, 'of', nrow(nloci), ':', myalcnt1, myalcnt2, 'Cannot run because not all sim files exist'))
		}
		if(myyr2=='1114'){
			print(paste(i, 'of', nrow(nloci), ':', myalcnt1, myalcnt2, myalcnt3, 'Cannot run because not all sim files exist'))
			warning(paste(i, 'of', nrow(nloci), ':', myalcnt1, myalcnt2, myalcnt3, 'Cannot run because not all sim files exist'))
		}
	}
}

# how many sample sizes to run?
sum(nloci$todo)
nrow(nloci) # compare to number we could run
	
# any sample sizes with < all loci to run?
nloci[nlocitorun<nloci & nlocitorun>0,]

# trim nloci to rows to do
nloci2 <- nloci[todo==1,]
nrow(nloci2)

# loci with large sample sizes?
nloci2[nlocitorun>1000,.(.N, max(nlocitorun))]

# break into batches of <= 400
nbatch <- ceiling(nrow(nloci2)/400)
nloci2[,batch:=0]
for(i in 1:nbatch){
	rows <- ((i-1)*400+1):min(i*400, nrow(nloci2))
	nloci2[rows,batch:=i]
}
print(paste(nbatch, 'batches to run'))
	
# write out the shell script for each batch
for(j in 1:nbatch){
	outfile <- paste('scripts/wfs_nullmodel_submit_all_sbatch', j, '.sh', sep='')
	print(paste('Wrote', outfile))

	cat('#!/bin/bash\n', file=outfile, append=FALSE) # header

	if(myyr2 != '1114'){
		for(i in which(nloci2[,batch==j])){
			cat(paste('sbatch --job-name=nlmd', nloci2[i,alcnt1], ',', nloci2[i,alcnt2], ' --cpus-per-task=', min(16, nloci2[i,nloci]), ' scripts/wfs_nullmodel_sbatch.sh ', nloci2[i,alcnt1], ' ', nloci2[i, alcnt2], ' ', pop, ' ', myyr1, ' ', myyr2, ' ', min(16, nloci2[i,nloci]), '\n', sep=''), file=outfile, append=TRUE)
		}
	}
	if(myyr2 == '1114'){ # 3 time points
		for(i in which(nloci2[,batch==j])){
			cat(paste('sbatch --job-name=nlmd', nloci2[i,alcnt1], ',', nloci2[i,alcnt2], ',', nloci2[i,alcnt3], ' --cpus-per-task=', min(16, nloci2[i,nloci]), ' scripts/wfs_nullmodel_sbatch_3times.sh ', nloci2[i,alcnt1], ' ', nloci2[i, alcnt2], ' ', pop, ' ', myyr1, ' ', myyr2, ' ', min(16, nloci2[i,nloci]), ' ', nloci2[i,alcnt3], '\n', sep=''), file=outfile, append=TRUE)
		}
	}
}