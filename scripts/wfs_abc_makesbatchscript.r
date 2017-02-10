# write out a shell script that will submit an sbatch job for every combination of sample sizes that we need

if(!grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table)
}
if(grepl('hpc.uio.no|login', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


load('analysis/wfs_targ.rdata') # targets normalized
load('analysis/wfs_sims_meansds.rdata') # useful as a list of sample sizes

# how many loci at each sample size
setkey(targ, alcnt1, alcnt2)
nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2)]
	nrow(nloci) # 185... why not 193?
	nrow(meansds)
	setdiff(meansds[,paste(alcnt1, alcnt2)], nloci[,paste(alcnt1, alcnt2)]) # 8 sample sizes not in targ (maybe because LG01, 02, etc were trimmed out?)

# sample sizes with >100 loci
nloci[nloci>5000,]

# check which sample sizes are already run
nloci[,todo:=1] # set up a flag for which we should run

for(i in 1:nrow(nloci)){
	myalcnt1 <- nloci[i,alcnt1]
	myalcnt2 <- nloci[i,alcnt2]
	# find the chunk of loci of the specified sample size to operate on
	samplepartinds <- targ[alcnt1==myalcnt1 & alcnt2==myalcnt2,locusnum] # locus numbers in each sample size chunk
	orig.n <- length(samplepartinds)

	existingfiles <- list.files(path='analysis/temp', pattern=paste('wfs_abc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''))
#	print(paste('found', length(existingfiles), 'existing file(s) relevant to this sample size'))
	existingrngs <- strsplit(gsub(paste('wfs_abc_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus|.csv.gz', sep=''), '', existingfiles), split='-') # extract just the locus ranges
	existingrngs <- lapply(existingrngs, as.numeric)
	existingloci <- numeric(0)
	if(length(existingrngs)>0) for(j in 1:length(existingrngs)) existingloci <- c(existingloci, existingrngs[[j]][1]:existingrngs[[j]][2])

	samplepartinds <- setdiff(samplepartinds, existingloci) # remove loci already run (from existing files)
	sampleparts.n <- length(samplepartinds) # number of loci in the samplesize chunk

	if(sampleparts.n==0){
		nloci[i,todo:=0]
	} else {
		print(paste(myalcnt1, myalcnt2, 'Run', sampleparts.n, 'loci of', orig.n, 'loci originally'))
	}
}
	
# write out the shell script
outfile <- 'scripts/wfs_abc_submit_all_sbatch2.sh'

cat('#!/bin/bash\n', file=outfile, append=FALSE) # header

for(i in which(nloci[,todo==1])){
	cat(paste('sbatch --job-name=abc', nloci[i,alcnt1], ',', nloci[i,alcnt2], ' --cpus-per-task=', min(16, nloci[i,nloci]), ' scripts/wfs_abc_sbatch.sh ', nloci[i,alcnt1], ' ', nloci[i, alcnt2], ' ', min(16, nloci[i,nloci]), '\n', sep=''), file=outfile, append=TRUE)
}