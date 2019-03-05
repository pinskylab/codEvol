# Calculate probability of null model producing results as extreme as our observations
# like wfs_nullmodel_function.r, but specifically for loci run with many simulations (those with low p-values)
# run after wfs_make_sims_null.r or wfs_make_sims_null_3times.r
# This version set up to run for Lof

pop <- 'Lof' # only one that works right now
yrs <- '071114' # only one that works right now

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
	if(yrs != '071114'){
		targfile <- paste('data_2018.09.05/Frequency_table_Lof', myyr1, '_Lof', myyr2, '.txt', sep='')
	}
	if(yrs == '071114'){
		targfile <- paste('data_2018.09.05/Frequency_table_Lof07_Lof11.txt', sep='')
		targfile2 <- paste('data_2018.09.05/Frequency_table_Lof07_Lof14.txt', sep='')
	}
}
# if(pop == 'Can'){
# 	targfile <- paste('data_2018.09.05/Frequency_table_CAN_40_TGA.txt', sep='')
# }
targ <- fread(targfile, header=TRUE)
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locusnum:=1:nrow(targ)] # add a locus number indicator

if(yrs=='071114'){
	targ2 <- fread(targfile2, header=TRUE)
	setnames(targ2, 3:7, c('alcnt1', 'f1samp', 'alcnt3', 'f3samp', 'ABS_DIFF13'))
	setkey(targ, CHROM, POS)
	setkey(targ2, CHROM, POS)
	targ <- targ[targ2,.(locusnum, CHROM, POS, alcnt1, alcnt2, alcnt3, f1samp, f2samp, f3samp)]
}

# trim only to loci with low p-values in a previous run
print('Trimming to loci with p<=8e-6') # 4 out of 500k sims
pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-11-14.rds'))
if(pop=='Lof'){
	targ <- merge(targ, pvals[,.(CHROM, POS, p)]) # merge in p-values for Lof
	print(nrow(targ)) # initial # loci
	targ <- targ[p<=8e-6,]
	print(nrow(targ)) # after trimming
} else {
	stop('Pop is not set to Lof!! We do not handle this case yet')
}

# trim only to loci with simulation files
simfiles <- list.files(path='analysis/temp', pattern='wfs_simsnull_ff.+_[[:digit:]]+\\.ffData') # the existing simulations
sims <- gsub('wfs_simsnull_ff|.ffData', '', simfiles) # trim file names to just sample sizes and rep #
sim_locs <- unique(gsub('_[[:digit:]]+', '', sims)) # trim to just the sample sizes
if(yrs!='071114') targ_locs <- targ[,paste(alcnt1, alcnt2, sep=',')]
if(yrs=='071114') targ_locs <- targ[,paste(alcnt1, alcnt2, alcnt3, sep=',')]
keep <- targ_locs %in% sim_locs
print(dim(targ)) # how many loci initially
sum(keep) # how many to keep
targ <- targ[keep,] # trim
print(dim(targ)) # how many kept (a double check)

# Null model test: how likely are results this extreme?
# locusnum: the locus number
# thistarg: f1samp and f2samp (first and second sample allele frequencies), as a data.table
# thisout.ff: ff object from the null model simulations, 7 x many, with rows for f1samp and f2samp. 
# tol: how close simulations and observation have to be to have the "same" starting frequency or to have the same change in frequency
nullmodtest <- function(thistarg, thisout.ff, tol=1/100){

	# find extreme simulations
	stinds <- abs(thisout.ff['f1samp',] - thistarg[,f1samp]) < tol # null model simulations that had the same starting sample frequency
	delta <- thistarg[,f2samp] - thistarg[,f1samp] # observed allele frequency change
	if(delta>0){
		extinds <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) >= (delta - tol) # simulations that also had as or greater allele frequency change
	} else {
		extinds <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) <= delta + tol # simulations that also had as or more extreme allele frequency change (for delta <= 0)
	}
	
	return(data.frame(e=sum(extinds), n=sum(stinds)))
}

# for 3 timepoints
nullmodtest3 <- function(thistarg, thisout.ff, tol=1/100){

	# find extreme simulations
	stinds <- abs(thisout.ff['f1samp',] - thistarg[,f1samp]) < tol # null model simulations that had the same starting sample frequency
	delta <- thistarg[,f2samp] - thistarg[,f1samp] # observed allele frequency change
	delta2 <- thistarg[,f3samp] - thistarg[,f1samp] # observed allele frequency change
	if(delta>0){
		extinds1 <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) >= (delta - tol) # simulations that also had as or greater allele frequency change
	}
	if(delta<=0){
		extinds1 <- (thisout.ff['f2samp',stinds] - thisout.ff['f1samp',stinds]) <= delta + tol # simulations that also had as or more extreme allele frequency change (for delta <= 0)
	}
	if(delta2>0){
		extinds2 <- (thisout.ff['f3samp',stinds] - thisout.ff['f1samp',stinds]) >= (delta2 - tol) # simulations that also had as or greater allele frequency change
	}
	if(delta2<=0){
		extinds2 <- (thisout.ff['f3samp',stinds] - thisout.ff['f1samp',stinds]) <= delta2 + tol # simulations that also had as or more extreme allele frequency change (for delta <= 0)
	}
	
	extinds <- extinds1 & extinds2 # is this the right way to combine the p-values?
	 
	return(data.frame(e=sum(extinds), n=sum(stinds)))
}




##############################################################
### run calculations for each row of targ with a low p-value
##############################################################
targ[, ':=' (pLof071114low=as.numeric(NA), pLof071114lowe=as.numeric(NA), pLof071114lown=as.numeric(NA))]

for(i in 1:nrow(targ)){
	print(paste('row', i, 'of', nrow(targ)))
	# find the appropriate files with null simulations
	if(yrs!='071114') targ_loc <- targ[i,paste(alcnt1, alcnt2, sep=',')]
	if(yrs=='071114') targ_loc <- targ[i,paste(alcnt1, alcnt2, alcnt3, sep=',')]
	ffnms <- list.files(path='analysis/temp', pattern=paste('wfs_simsnull_ff', targ_loc, '_[[:digit:]]+\\.ffData', sep=''), full.names=TRUE) # the existing simulations
		# could add in the file: pattern=paste('wfs_simsnull_ff', targ_loc, '\\.ffData', sep='')
	ffnms <- gsub('.ffData', '', ffnms)
	
	res <- data.frame(e=numeric(0),n=numeric(0)) # vector of number of simulations more extreme than observations and count of simulations with same starting sample freq
	
	# loop through the files of simulations
	for(j in 1:length(ffnms)){
		print(paste('     file', j, 'of', length(ffnms)))
		# load one file of simulations
		ffload(ffnms[j], overwrite=TRUE) # takes 30 sec or so. loads thisout.ff

		if(yrs != '071114') thisres <- nullmodtest(thistarg=targ[i,.(f1samp, f2samp)], thisout.ff=thisout.ff) # run null model test. returns vector of number of simulations more extreme than observations and count of simulations with same starting sample freq
		if(yrs == '071114') thisres <- nullmodtest3(thistarg=targ[i,.(f1samp, f2samp, f3samp)], thisout.ff=thisout.ff) # 3 time points

		# save results from this loop
		res <- rbind(res, thisres)

		# remove ff temporary files
		delete(thisout.ff)
		rm(thisout.ff)
	}
	
 	# calculate p-value (see North et al. 2002 Am J Hum Gen for why the +1s) 	
 	# add to targ
 	targ[i, ':=' (pLof071114low=(sum(res$e)+1)/(sum(res$n)+1), pLof071114lowe=sum(res$e), pLof071114lown=sum(res$n))]

}


# need to write out targ again
if(pop=='Lof'){
	outfile <- 'analysis/wfs_nullmodel_outliers_lowp_Lof_07-11-14.tsv.gz'
	write.table(targ, file=gzfile(outfile), sep='\t', row.names=FALSE, quote=FALSE)
}