# Make simulations for the null model
# best to run on cod node with nohup if many sample sizes to make
# Arguments (ONLY SOME ARE IMPLEMENTED RIGHT NOW, SEE wfs_make_sims_null_3times.r FOR EXAMPLE)
#	pop: Can or Lof
#	yr1: 07 or 11 (not needed for Can)
#	yr2: 11_14 to run both 07-11 and 07-14 (not needed for Can)
#	trimlowsampsize: 1 to trim loci with less than half of individuals genotypes, 0 to not trim
#	rerunlow: 0 for default, 1 to only run simulations for loci with low p-values (<=4/(n+1)) in wfs_nullmodel_pos&pvals_*.rds

# read arguments
# args=(commandArgs(TRUE))
# if(length(args)==0){
#     print("No arguments supplied.")
#     pop <- 'Lof'
#     yr1 <- '07'
#     yr2 <- '14'
# 	trimlowsampsize <- 1
# 	rerunlow <- 0
# }else{
# 	for(i in 1:length(args)){
# 		eval(parse(text=args[[i]]))
# 	}
# 	if(!(pop %in% c('Lof', 'Can'))) stop('pop must be Lof or Can!')
# 	if(pop=='Lof' & !(yr1 %in% c('07', '11'))) stop('yr1 must be 07 or 11 for Lof!')
# 	if(pop=='Lof' & !(yr2 %in% c('11', '14'))) stop('yr2 must be 11 or 14 for Lof!')
# 	if(!(trimlowsampsize %in% 0:1)) stop('trimlowsampsize must be 0 or 1!')
# 	if(!(rerunlow %in% 0:1)) stop('rerunlow must be 0 or 1!')
# }

# pop<-'Lof'; yr1<-'07';yr2<-'11_14';trimlowsampsize<-0;rerunlow<-0 # initial run for all Lof loci 07-11 and 07-14
pop<-'Lof'; yr1<-'07';yr2<-'11_14';trimlowsampsize<-0;rerunlow<-1 # re-run low p-values for all Lof loci 07-11 and 07-14
# pop<-'Can'; yr1<-'00';yr2<-'00';trimlowsampsize<-0;rerunlow<-0 # initial run for all CAN loci
# pop<-'Can'; yr1<-'00';yr2<-'00';trimlowsampsize<-0;rerunlow<-1 # re-run low p-values for all CAN loci

print(paste('Arguments: pop=', pop, ', yr1=', yr1, ', yr2=', yr2, ', trimlowsampsize=', trimlowsampsize, ', rerunlow=', rerunlow, sep=''))


# load functions
source('scripts/wfs_byf1samp.r')
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


# read in Ne and frequency data
if(pop=='Lof' & yr1=='07' & yr2=='11'){
	nes <- read.table('analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	nchrs <- fread('data_2019_03_18/Frequency_table_Lof07_Lof11.txt', header=TRUE) # read in frequency table data
	gen <- 11
}
if(pop=='Lof' & yr1=='07' & yr2=='14'){
	nes <- read.table('analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	nchrs <- fread('data_2019_03_18/Frequency_table_Lof07_Lof14.txt', header=TRUE) # read in frequency table data
	gen <- 11
}
if(pop=='Lof' & yr1=='11' & yr2=='14'){
	nes <- read.table('analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	nchrs <- fread('data_2019_03_18/Frequency_table_Lof11_Lof14.txt', header=TRUE) # read in frequency table data
	gen <- 11
}
if(pop=='Lof' & yr1=='07' & yr2=='11_14'){ # both 07-11 and 07-14
	nes <- read.table('analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	nchrs11 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof11.txt', header=TRUE) # read in frequency table data
	nchrs14 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof14.txt', header=TRUE) # read in frequency table data
	nchrs11[,yr:=11]
	nchrs14[,yr:=14]
	nchrs <- rbind(nchrs11, nchrs14) # all possible sample size combinations across both pairs of years
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
	print('Trimming to loci with p<=4/(n+1)')
	if(pop=='Lof' & yr2=='11_14'){
		pvals11 <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-11.rds'))
		pvals14 <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-14.rds'))
		pvals11[,yr:=11]
		pvals14[,yr:=14]
		pvals <- rbind(pvals11, pvals14)
		nchrs <- merge(nchrs, pvals[,.(CHROM, POS, n, p, yr)], by=c('CHROM', 'POS', 'yr')) # merge in p-values
	}
	if(pop=='Can'){
		pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_Can.rds')) # NOT SURE YET THIS WORKS FOR CAN
		nchrs <- merge(nchrs, pvals[,.(CHROM, POS, n, p)], by=c('CHROM', 'POS')) # merge in p-values
	}
	print(nrow(nchrs))		
	nchrs <- nchrs[p<=4/(n+1),]
	print(nrow(nchrs))
	print(paste('retained p=', paste(signif(nchrs[,.(min(p), max(p))], digits=2), collapse=' to '), ', n=', paste(nchrs[,.(min(n), max(n))], collapse=' to '), sep=''))
}



# parameters for the simulation
nsims <- 500000 # number of sims to run for each starting frequency in each sample size
nreps <- 1 # number of repetitions of the number of simulations
if(rerunlow==1) nreps <- 10*nchrs[,max(n)]/nsims # run 10x more sims than already run
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
	


# check that nes >0 
print(summary(nes))
print(paste('Of', length(nes), 'Ne values', sum(nes > 0), 'are >0. Trimming all others out.'))
nes <- nes[nes>0]

# check temp directory for previous simulations (if Ne, gen and other parameters are the same)
if(pop=='Lof'){ 
	torun <- paste(paste(c1s, c2s, sep=','), reps, sep='_') # the list of sample sizes we need to run
	existing <- list.files(path='analysis/temp', pattern='wfs_simsnull_ff.+_[[:digit:]]+\\.ffData') # the existing simulations
	existing <- gsub('wfs_simsnull_ff|.ffData', '', existing) # trim file names to just sample sizes and rep#

	print('Trimming out files already run')
	keep <- !(torun %in% existing) # the sample sizes that still need to be run
	c1s <- c1s[keep]
	c2s <- c2s[keep]
	reps <- reps[keep]
}
if(pop=='Can'){
	torun <- paste(paste(c1s, c2s, sep=','), reps, sep='_') # the list of sample sizes we need to run
	existing <- list.files(path='analysis/temp', pattern='wfs_simsnullCAN_ff.+_[[:digit:]]+\\.ffData') # the existing simulations. works on new and old simulations because I renamed all older sim files to have _1 as of March 5 2019.
	existing <- gsub('wfs_simsnullCAN_ff|.ffData', '', existing) # trim file names to just sample sizes and rep#

	print('Trimming out files already run')
	keep <- !(torun %in% existing) # the sample sizes that still need to be run
	c1s <- c1s[keep]
	c2s <- c2s[keep]
	reps <- reps[keep]
}

print(length(c1s)) # how many to run?
print(length(c2s))
print(length(reps))

# combine into one data.table
torun <- data.table(c1s=c1s, c2s=c2s, reps=reps)

# break into batches of <= 400
nbatch <- ceiling(nrow(torun)/400)
torun[,batch:=0]
for(i in 1:nbatch){
	rows <- ((i-1)*400+1):min(i*400, nrow(torun))
	torun[rows,batch:=i]
}
print(paste(nbatch, 'batches to run'))
	
# write out the shell script for each batch
for(j in 1:nbatch){
	outfile <- paste('scripts/wfs_make_sims_null_submit_sbatch', j, '.sh', sep='')
	print(paste('Wrote', outfile))

	cat('#!/bin/bash\n', file=outfile, append=FALSE) # header

	for(i in which(torun[,batch==j])){
		cat(paste('sbatch --job-name=nlsm', torun[i,c1s], ',', torun[i,c2s], ' --cpus-per-task=1 scripts/wfs_make_sims_null_sbatch.sh ', pop, ' ', torun[i,c1s], ' ', torun[i, c2s], ' ', yr1, ' ', yr2, ' ', trimlowsampsize, ' ', rerunlow, ' ', torun[i, reps], ' ', nsims, '\n', sep=''), file=outfile, append=TRUE)
	}
}
