# Make null model simulations
# Set up to be called by wfs_make_sims_null_sbatch.sh
# arguments: pop c1 c2 yr1 yr2 trimlowsampsize rerunlow repnum nsims
#	pop: Can or Lof
#	c1: the first sample size (# chromosomes)
#	c2: the second sample size (# chromosomes)
#	yr1: 07 or 11 (not needed for Can)
#	yr2: 11 or 14 (not needed for Can)
#	trimlowsampsize: 1 to trim loci with less than half of individuals genotypes, 0 to not trim
#	rerunlow: 0 for default, 1 to only run simulations for loci with low p-values (<=4/(n+1)) in wfs_nullmodel_pos&pvals_07-11-14.rds
#	repnum: which repetition number to run
# 	nsims: usually 500k

# read arguments
args=(commandArgs(trailingOnly=TRUE))

pop <- args[1]
c1 <- as.numeric(args[2])
c2 <- as.numeric(args[3])
yr1 <- as.numeric(args[4])
yr2 <- as.numeric(args[5])
trimlowsampsize <- as.numeric(args[6])
rerunlow <- as.numeric(args[7])
repnum <- as.numeric(args[8])
nsims <- as.numeric(args[9])

if(!(pop %in% c('Lof', 'Can'))) stop('pop must be Lof or Can!')
if(pop=='Lof' & !(yr1 %in% c('07', '11'))) stop('yr1 must be 07 or 11 for Lof!')
if(pop=='Lof' & !(yr2 %in% c('11', '14'))) stop('yr2 must be 11 or 14 for Lof!')
if(!(trimlowsampsize %in% 0:1)) stop('trimlowsampsize must be 0 or 1!')
if(!(rerunlow %in% 0:1)) stop('rerunlow must be 0 or 1!')

print(paste('Arguments: pop=', pop, ', c1=', c1, ', c2=', c2, ', yr1=', yr1, ', yr2=', yr2, ', trimlowsampsize=', trimlowsampsize, ', repnum=', repnum, ', nsims=', nsims, sep=''))

# load functions: assume this is run on a cod or abel node
source('scripts/wfs_byf1samp.r')
require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
require(data.table, lib.loc="/projects/cees/lib/R_packages/")


# read in Ne and frequency data (latter for starting frequencies)
if(pop=='Lof'){
	nes <- read.table('analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	freqfile <- paste('data_2018.09.05/Frequency_table_Lof', yr1, '_Lof', yr2, '.txt', sep='')
	nchrs <- fread(freqfile, header=TRUE) # read in frequency table data
	gen <- 11
}
if(pop=='Can') {
	nes <- read.table('analysis/Can_40_to_Can.w_Ne_bootstrap.txt')[,1] # the values of Ne from wfabc_1
	freqfile <- paste('data_2018.09.05/Frequency_table_CAN_40_TGA.txt', sep='')
	nchrs <- fread(freqfile, header=TRUE)
	gen <- 8
}

setnames(nchrs, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))

# trim out missing loci
	print(paste(nrow(nchrs), 'loci initially'))
nchrs <- nchrs[N_CHR_1>0 & N_CHR_2>0,]
	print(paste(nrow(nchrs), 'loci after trimming out missing loci'))

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
	if(pop=='Lof'){
		pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_07-11-14.rds'))
	}
	if(pop=='Can') {
		pvals <- as.data.table(readRDS('analysis/wfs_nullmodel_pos&pvals_Can.rds'))
	}
	nchrs <- merge(nchrs, pvals[,.(CHROM, POS, n, p)], by=c('CHROM', 'POS')) # merge in p-values
	print(nrow(nchrs))		
	nchrs <- nchrs[p<=4/(n+1),]
	print(nrow(nchrs))
}


# check that nes >0 
print(summary(nes))
print(paste('Of', length(nes), 'Ne values', sum(nes > 0), 'are >0. Trimming all others out.'))
nes <- nes[nes>0]


# check where this is to monitor temp file creation in the file system /tmp/Rtmpp__
options('fftempdir')


# only run observed starting sample frequencies
nchr1s <- nchrs[c1==N_CHR_1 & c2==N_CHR_2, unique(round(Freq_1*N_CHR_1))] # number of starting copies of allele for the retained loci
f1samps=rep((nchr1s)/c1, rep(nsims, length(nchr1s)))	

print(paste(length(nchr1s), 'starting frequencies to run'))
print(paste(length(f1samps), 'simulations to run'))

# make simulations	
thisout <- sapply(f1samps, FUN=wfs_byf1samp, smin=0, smax=0, c1=c1, c2=c2, gen=gen, ne=nes, h=0.5, simplify=TRUE)

thisout.ff <- ff(thisout, dim=dim(thisout), dimnames=dimnames(thisout)) # create in tempdir

# save to permanent file (semi-permanent.. in my temp directory)
if(pop=='Lof') {
	outfile <- paste('analysis/temp/wfs_simsnull_ff', paste(c1, c2, sep=','), '_', repnum, sep='')
}
if(pop=='Can') {
	outfile=paste('analysis/temp/wfs_simsnullCAN_ff', paste(c1, c2, sep=','), '_', repnum, sep='')
}
ffsave(thisout.ff, file=outfile)
print(paste('saved', outfile))

# remove the temp files
delete(thisout.ff) # returns TRUE
rm(thisout.ff)


