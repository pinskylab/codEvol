# Null model test on one locus
# run after wfs_make_sims_null.r

# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
	require(ff)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
	require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# load data
load('analysis/wfs_targ.rdata') # target statistics

# Null model test: how likely are results this extreme?
nullmodtest <- function(locusnum, thistarg, thisout.ff){

	# find extreme simulations
	stinds <- thisout.ff['f1samp',] == thistarg[,f1samp] # null model simulations that had the same starting sample frequency
	delta <- thistarg[,f2samp] - thistarg[,f1samp]
	if(delta>0){
		extinds <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) >= delta # simulations that also had as or greater allele frequency change
	} else {
		extinds <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) <= delta # simulations that also had as or more extreme allele frequency change
	}
	
 	# calculate p-value
 	p <- sum(extinds)/sum(stinds)
 
	return(c(locusnum=locusnum, p=p, n=sum(stinds)))
}

###################################################
### run calculations in parallel with foreach
###################################################

locnum <- 73; sampsize <- '50,44'

# load data for null model
ffnm2 <- paste('analysis/temp/wfs_simsnull_ff', sampsize, sep='')
ffload(ffnm2, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff


thistarg <- targ[locnum,] # the input object
	
# run abc
results <- nullmodtest(locusnum=thistarg[,locusnum], thistarg=thistarg[,.(f1samp, f2samp)], thisout.ff=thisout.ff) # run null model test


	results

	# remove ff temporary files
	delete(thisout.ff)
	rm(thisout.ff)
