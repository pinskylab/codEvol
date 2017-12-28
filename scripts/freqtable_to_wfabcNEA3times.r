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
outfile <- 'analysis/LOF_07_to_LOF_S_11_to_LOF_S_14.wfabc'
gens =c(0,11,11) # for 1907 vs. 2011 vs. 2014

# read in data
	# NEA 1907-2014
datNEA <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE)
	setnames(datNEA, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF0714'))

	# NEA 1907-2011
datNEA11 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE)
	setnames(datNEA11, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF0711'))


# read in 25kmer filter
#loc25NEA <- fread('data_2017.11.24/Norway_25K_mer_positions.txt')

# trim to loci that meet 25kmer filter
#setkey(datNEA, CHROM, POS)
#setkey(datNEA11, CHROM, POS)
#setkey(loc25NEA, CHROM, POS)
#	nrow(datNEA)
#	nrow(datNEA11)
#datNEA <- datNEA[loc25NEA, nomatch=0] # nomatch=0 so that non-matching rows are dropped
#datNEA11 <- datNEA11[loc25NEA, nomatch=0] # nomatch=0 so that non-matching rows are dropped
#	nrow(datNEA)
#	nrow(datNEA11)
	
# trim out inversions and Unplaced
datNEA <- datNEA[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
datNEA11 <- datNEA11[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

	dim(datNEA)
	dim(datNEA11)

# merge
setkey(datNEA, CHROM, POS, N_CHR_07, Freq_07)
setkey(datNEA11, CHROM, POS, N_CHR_07, Freq_07)
dat <- datNEA[datNEA11, .(CHROM, POS, N_CHR_07, N_CHR_11, N_CHR_14, Freq_07, Freq_11, Freq_14)]
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

