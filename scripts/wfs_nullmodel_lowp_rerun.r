# Calculate probability of null model producing results as extreme as our observations
# like wfs_nullmodel_function.r, but specifically for loci run with many simulations (those with low p-values)
# run after wfs_make_sims_null.r or wfs_make_sims_null_3times.r
# This version set up to run for Lof

# pop <- 'Lof'; yrs <- '071114' # for Lof 3 time points
# pop <- 'Lof'; yrs <- '0711'
#pop <- 'Lof'; yrs <- '0714'
 pop <- 'Can'; yrs <- '00'

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
	if(yrs == '0711'){
		targfile <- 'data_2019_03_18/Frequency_table_Lof07_Lof11.txt'
	}
	if(yrs == '0714'){
		targfile <- 'data_2019_03_18/Frequency_table_Lof07_Lof14.txt'
	}
	if(yrs == '071114'){
		targfile <- 'data_2019_03_18/Frequency_table_Lof07_Lof11.txt'
		targfile2 <- 'data_2019_03_18/Frequency_table_Lof07_Lof14.txt'
	}
}
if(pop == 'Can'){
	targfile <- paste('data_2019_03_18/Frequency_table_CAN_40_TGA.txt', sep='')
}
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
print('Trimming to loci with p<=4/(n+1)') # 4 out of n+1 sims
if(pop=='Lof' & yrs=='0711'){
	pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-11.rds'))
	targ <- merge(targ, pvals[,.(CHROM, POS, n, p)], by=c('CHROM', 'POS')) # merge in p-values for Lof
}
if(pop=='Lof' & yrs=='0714'){
	pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-14.rds'))
	targ <- merge(targ, pvals[,.(CHROM, POS, n, p)], by=c('CHROM', 'POS')) # merge in p-values for Lof
}
if(pop=='Lof' & yrs=='1114'){
	pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_11-14.rds'))
	targ <- merge(targ, pvals[,.(CHROM, POS, n, p)], by=c('CHROM', 'POS')) # merge in p-values for Lof
}
if(pop=='Lof' & yrs=='071114'){
	pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-11-14.rds'))
	targ <- merge(targ, pvals[,.(CHROM, POS, n, p)], by=c('CHROM', 'POS')) # merge in p-values for Lof
}
if(pop=='Can'){
	pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_Can.rds'))
	targ <- merge(targ, pvals[,.(CHROM, POS, n, p)], by=c('CHROM', 'POS')) # merge in p-values for Lof
}
print(nrow(targ)) # initial # loci
targ <- targ[p<=4/(n+1),]
print(nrow(targ)) # after trimming


# trim only to loci with simulation files
if(pop=='Lof'){
	simfiles <- list.files(path='analysis/temp', pattern='wfs_simsnull_ff.+_[[:digit:]]+\\.ffData') # the existing simulations
	sims <- gsub('wfs_simsnull_ff|.ffData', '', simfiles) # trim file names to just sample sizes and rep #
}
if(pop=='Can'){
	simfiles <- list.files(path='analysis/temp', pattern='wfs_simsnullCAN_ff.+_[[:digit:]]+\\.ffData') # the existing simulations
	sims <- gsub('wfs_simsnullCAN_ff|.ffData', '', simfiles) # trim file names to just sample sizes and rep #
}
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

for(i in 1:nrow(targ)){
	print(paste('row', i, 'of', nrow(targ)))
	# find the appropriate files with null simulations
	if(yrs!='071114') targ_loc <- targ[i,paste(alcnt1, alcnt2, sep=',')]
	if(yrs=='071114') targ_loc <- targ[i,paste(alcnt1, alcnt2, alcnt3, sep=',')]
	if(pop=='Lof') ffnms <- list.files(path='analysis/temp', pattern=paste('wfs_simsnull_ff', targ_loc, '_[[:digit:]]+\\.ffData', sep=''), full.names=TRUE) # the existing simulations
	if(pop=='Can') ffnms <- list.files(path='analysis/temp', pattern=paste('wfs_simsnullCAN_ff', targ_loc, '_[[:digit:]]+\\.ffData', sep=''), full.names=TRUE) # the existing simulations
		# could add in the file: pattern=paste('wfs_simsnull_ff', targ_loc, '\\.ffData', sep='')
	ffnms <- gsub('.ffData', '', ffnms)
	
	res <- data.frame(e=numeric(0),n=numeric(0)) # vector of number of simulations more extreme than observations and count of simulations with same starting sample freq
	
	# loop through the files of simulations
	for(j in 1:length(ffnms)){
		print(paste('     file', j, 'of', length(ffnms)))
		# load one file of simulations
		ffload(ffnms[j], overwrite=TRUE, rootpath=getOption('fftempdir')) # takes 30 sec or so. loads thisout.ff

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
 	targ[i, ':=' (p=(sum(res$e)+1)/(sum(res$n)+1), n=sum(res$n))]


}

# insert new p-values and n sims into pvals dataframe (probably a more elegant way to do this)
for(i in 1:nrow(targ)){
	j <- pvals[,CHROM]==targ[i,CHROM] & pvals[,POS]==targ[i,POS]
	pvals[j, p:=targ[i,p]]
	pvals[j, n:=targ[i,n]]
}

# recalculate p.adj and p.adj3 (from wfs_nullmodel_combine.r and wfs_nullmodel_analysis.r)
# pvals[,p.adj := p.adjust(p, method='fdr')]
# pvals[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),p.adj3 := p.adjust(p, method='fdr')]

# examine
summary(pvals)
pvals[p==min(p),.(CHROM, POS, p, r=p*(n+1)-1, n=n)] # low p-values
pvals[,min(p*(n+1)-1)] # if <=3, have more to re-run
pvals[p*(n+1)-1 <= 3,.(CHROM, POS, nloci=.N, r=p*(n+1)-1, n=n)] # number left to run

# need to write out pos&pvals again
if(pop=='Lof' & yrs=='071114') outfile <- 'analysis/wfs_nullmodel_pos&pvals_07-11-14.rds'
if(pop=='Lof' & yrs=='0711') outfile <- 'analysis/wfs_nullmodel_pos&pvals_07-11.rds'
if(pop=='Lof' & yrs=='0714') outfile <- 'analysis/wfs_nullmodel_pos&pvals_07-14.rds'
if(pop=='Lof' & yrs=='1114') outfile <- 'analysis/wfs_nullmodel_pos&pvals_11-14.rds'
if(pop=='Can') outfile <- 'analysis/wfs_nullmodel_pos&pvals_Can.rds'
print(outfile)
# saveRDS(pvals[,.(CHROM, POS, n, p, p.adj, p.adj3)], file=outfile)
saveRDS(pvals, file=outfile)
