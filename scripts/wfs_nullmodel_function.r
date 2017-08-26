# Calculate probability of null model producing results as extreme as our observations
# run after wfs_make_sims.r/wfs_process_sims.r and wfs_make_sims_null.r
# This version set up to run for Lof or Can for a single sample size, taken as command line arguments
# Set up to run for a specified pair of years (07, 11, or 14) if Lof is chosen
# runs on cod node as a script with arguments (myalcnt1 myalcnt2 pop myyr1 myyr2 kmer maxcores)

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args)<6) {
  stop("Have to specify myalcnt1, myalcnt2, pop, myyr1, myyr2, kmer", call.=FALSE)
} else if (length(args)==6){
	maxcores <- 16 # default maximum cores (for abel)
} else if (length(args)>6) {
	maxcores <- as.numeric(args[7])
}

myalcnt1 <- as.numeric(args[1])
myalcnt2 <- as.numeric(args[2])
pop <- args[3]
myyr1 <- args[4]
myyr2 <- args[5]
kmer <- as.numeric(args[6])

if(!(pop %in% c('Lof', 'Can'))){
	stop('pop must be one of Lof or Can', call.=FALSE)
}

if(!(myyr1 %in% c('07', '11')) & pop == 'Lof'){
	stop('myyr1 must be one of 07 or 11 if pop is Lof', call.=FALSE)
}

if(!(myyr2 %in% c('11', '14')) & pop == 'Lof'){
	stop('myyr2 must be one of 11 or 14 if pop is Lof', call.=FALSE)
}

if(!(kmer %in% c(25, 150))){
	stop('kmer must be one of 25 or 150', call.=FALSE)
}

print(paste('myalcnt1', myalcnt1, 'myalcnt2', myalcnt2, 'pop', pop, 'myyr1', myyr1, 'myyr2', myyr2, 'kmer', kmer, 'maxcores', maxcores))
print(Sys.info()["nodename"])


# set parameters
lociperpart <- 20000 # how many loci in each part (each part will be written to file)


# load functions: assume this is run on a cod or abel node
require(parallel, lib.loc="/projects/cees/lib/R_packages/")
require(iterators, lib.loc="/projects/cees/lib/R_packages/") # used with foreach
require(foreach, lib.loc="/projects/cees/lib/R_packages/") # for parallel procssing using foreach()
require(doMC, lib.loc="/projects/cees/lib/R_packages/")
require(doParallel, lib.loc="/projects/cees/lib/R_packages/") # needed for foreach
require(doSNOW, lib.loc="/projects/cees/lib/R_packages/") # set up SOCK cluster
require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
require(data.table, lib.loc="/projects/cees/lib/R_packages/")


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

# set up file names
existingfilespattern <- paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep='') # regexp to test for existing p-value files
existingfilesgsubpattern <- paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus|.csv.gz', sep='') # regexp to strip out all but locus numbers from existing file names
if(pop == 'Lof'){
	ffnm <- paste('analysis/temp/wfs_simsnull_ff', paste(myalcnt1, myalcnt2, sep=','), sep='') # ff file name for simulations
}
if(pop == 'Can'){
	ffnm <- paste('analysis/temp/wfs_simsnullCAN_ff', paste(myalcnt1, myalcnt2, sep=','), sep='')
}

# Null model test: how likely are results this extreme?
# locusnum: the locus number
# thistarg: f1samp and f2samp (first and second sample allele frequencies), as a data.table
# thisout.ff: ff object from the null model simulations, 7x10million, with rows for f1samp and f2samp. 
# tol: how close simulations and observation have to be to have the "same" starting frequency or to have the same change in frequency
nullmodtest <- function(locusnum, thistarg, thisout.ff, tol=1/100){

	# find extreme simulations
	stinds <- abs(thisout.ff['f1samp',] - thistarg[,f1samp])<tol # null model simulations that had the same starting sample frequency
	delta <- thistarg[,f2samp] - thistarg[,f1samp]
	if(delta>0){
		extinds <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) >= delta - tol # simulations that also had as or greater allele frequency change
	} else {
		extinds <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) <= delta + tol # simulations that also had as or more extreme allele frequency change
	}
	
 	# calculate p-value
 	p <- sum(extinds)/sum(stinds)
 
	return(c(locusnum=locusnum, p=p, n=sum(stinds)))
}




###################################################
### run calculations in parallel with foreach
###################################################

# find if any loci of this sample size have already been run (for checkpointing)
existingfiles <- list.files(path='analysis/temp', pattern=existingfilespattern)
print(paste('found', length(existingfiles), 'existing files relevant to this sample size'))
existingrngs <- strsplit(gsub(existingfilesgsubpattern, '', existingfiles), split='-') # extract just the locus ranges
existingrngs <- lapply(existingrngs, as.numeric)
existingloci <- numeric(0)
if(length(existingrngs)>0) for(i in 1:length(existingrngs)) existingloci <- c(existingloci, existingrngs[[i]][1]:existingrngs[[i]][2]) # we have run all loci between the first and last locunum listed in the file name

# find the chunk of loci of the specified sample size to operate on
samplepartinds <- targ[alcnt1==myalcnt1 & alcnt2==myalcnt2,locusnum] # locus numbers in each sample size chunk
orig.n <- length(samplepartinds)
samplepartinds <- setdiff(samplepartinds, existingloci) # remove loci already run (from existing files)
sampleparts.n <- length(samplepartinds) # number of loci in the samplesize chunk
print(paste('Will run null model test for', sampleparts.n, 'loci of', orig.n, 'loci originally'))

ndigits <- nchar(as.character(nrow(targ))) # used for formatting file names

if(sampleparts.n > 0){
	# set up cluster
	ncores <- min(c(maxcores, sampleparts.n)) # don't set up any more cores than we need for this batch of loci
	print(paste('Using', ncores, 'cores'))

	.libPaths('/projects/cees/lib/R_packages/')
	cl <- makeCluster(rep('localhost', ncores), type='SOCK') # make the cluster on the localhost
	registerDoSNOW(cl) # register the cluster so accessible to foreach
	clusterCall(cl, function(x) .libPaths(x), .libPaths()) # add libpaths to each node


	# set up data chunks to operate on and write out separately within this chunk of sample sizes
	parts <- 1:ceiling(sampleparts.n/lociperpart) # 1 to number of parts for this samplepart
	partinds <- vector('list', length(parts))
	for(i in 1:length(partinds)) partinds[[i]] <- samplepartinds[((i-1)*lociperpart+1):min(c(i*lociperpart, sampleparts.n))] # indices for the loci in each part

	# load the appropriate file of abc simulations
		# null model
	ffload(ffnm, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff

	# loop through each chunk of data, writing each to file
	for(partnum in parts){
		print(paste('partnum', partnum, 'of', length(parts)))
		thisdat <- iter(targ[locusnum %in% partinds[[partnum]],], by='row') # the iterator object: this part of the loop's chunk of data

		results <- foreach(i=thisdat, .combine=rbind, .packages=c('ff', 'doParallel', 'data.table')) %dopar% {
			nullmodtest(locusnum=i[,locusnum], thistarg=i[,.(f1samp, f2samp)], thisout.ff=thisout.ff) # run null model test. returns vector of locusnum, p-value, and count of simulations with same starting sample freq
		}

		minloc <- formatC(min(partinds[[partnum]]), width=ndigits, flag='0')
		maxloc <- formatC(max(partinds[[partnum]]), width=ndigits, flag='0')
		outfile <- paste('analysis/temp/wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus', minloc, '-', maxloc, '.csv.gz', sep='')
		write.csv(results, file=gzfile(outfile), row.names=FALSE) # write directly to gzipped file. has locusnum, p, and n (number of simulations)
		print(paste('wrote', outfile))
	}

	# remove ff temporary files
	delete(thisout.ff)
	rm(thisout.ff)

	# stop cluster
	stopCluster(cl)
} else {
	print('Exiting because all loci had already been run')
}