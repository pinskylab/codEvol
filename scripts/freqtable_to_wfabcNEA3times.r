# Convert data from Frequency table from Bastiaan Star to input to WFABC_1
# Use 3 time points for NEA: 1907, 2011, and 2014

# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
	require(plyr)
	require(parallel)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
}

# set up parameters
ntimes <- 3
outfile <- 'analysis/LOF_07_to_LOF_S_11_to_LOF_S_14_25kmer_dp.wfabc' # w/ 25kmer trimming
#outfile <- 'analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.wfabc' # w/out 25kmer trimming
gens =c(0,11,11) # for 1907 vs. 2011 vs. 2014


# read in data
if(outfile == 'analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.wfabc'){
	print('Not using 25kmer filter')
	datNEA <- fread('data_2018.09.13/Frequency_table_Lof07_Lof14.txt', header=TRUE) # 25kmer and individual read depth filter not applied
		setnames(datNEA, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF0714'))
	datNEA11 <- fread('data_2018.09.13/Frequency_table_Lof07_Lof11.txt', header=TRUE)
		setnames(datNEA11, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF0711'))
}
if(outfile == 'analysis/LOF_07_to_LOF_S_11_to_LOF_S_14_25kmer_dp.wfabc'){
	print('Using 25kmer filter')
	datNEA <- fread('data_2018.09.05/Frequency_table_Lof07_Lof14.txt', header=TRUE) # 25kmer and individual read depth filter applied
		setnames(datNEA, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF0714'))
	datNEA11 <- fread('data_2018.09.05/Frequency_table_Lof07_Lof11.txt', header=TRUE)
		setnames(datNEA11, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF0711'))
}


# trim out inversions and Unplaced
datNEA <- datNEA[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
datNEA11 <- datNEA11[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

dim(datNEA)
dim(datNEA11)

# trim out only Unplaced
#datNEA <- datNEA[!(CHROM %in% c('Unplaced')),]
#datNEA11 <- datNEA11[!(CHROM %in% c('Unplaced')),]
#
#	dim(datNEA)
#	dim(datNEA11)

# merge
setkey(datNEA, CHROM, POS, N_CHR_07, Freq_07)
setkey(datNEA11, CHROM, POS, N_CHR_07, Freq_07)
dat <- datNEA[datNEA11, .(CHROM, POS, N_CHR_07, N_CHR_11, N_CHR_14, Freq_07, Freq_11, Freq_14)]
	dim(dat)
	dim(datNEA)
	dim(datNEA11)

# trim out missing loci
dat <- dat[N_CHR_07 > 0 & N_CHR_11 > 0 & N_CHR_14 > 0,]
	dim(dat)
	dim(datNEA)
	dim(datNEA11)
	
# figure out number of sites
nsitesNEA <- nrow(dat)

# figure out # chromosomes in each sample of the focal allele
dat[,alcnt07:=round(N_CHR_07*Freq_07)]
dat[,alcnt11:=round(N_CHR_11*Freq_11)]
dat[,alcnt14:=round(N_CHR_14*Freq_14)]

# make a data table to write out
sampsizesNEA <- dat[,.(N_CHR_07, N_CHR_11, N_CHR_14)] # sample sizes
sampsizesNEA[,index:=seq(1,length.out=nrow(dat),by=2)] # index to allow sorting with allele counts (interleaving)
alcntsNEA <- dat[,.(alcnt07, alcnt11, alcnt14)] # allele counts
alcntsNEA[,index:=seq(2,length.out=nrow(dat),by=2)] # index, offset by 1 compared to sample sizes
setnames(sampsizesNEA, 1:3, c('X1', 'X2', 'X3')) # set column names to allow rbinding
setnames(alcntsNEA, 1:3, c('X1', 'X2', 'X3'))
outNEA <- rbind(sampsizesNEA, alcntsNEA)
setkey(outNEA, index) # sort by index: each pair of rows is now sample sizes followed by allele counts

# write out
cat(paste(nsitesNEA, ' ', ntimes, '\n', sep=''), file=outfile) # header line 1
cat(paste(paste(gens, collapse=','), '\n', sep=''), file=outfile, append=TRUE) # write out generation of each sample
write.table(outNEA[,.(X1, X2, X3)], file=outfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',') # append the rows of data

