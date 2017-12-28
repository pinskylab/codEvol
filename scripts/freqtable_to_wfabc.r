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
ntimes <- 2

# read in data
	# NEA 1907-2014
datNEA <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE)
	setnames(datNEA, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))
	outfileNEA <- 'analysis/LOF_07_to_LOF_S_14.wfabc'
	genNEA=11 # for 1907 vs. 2014. sample sizes

	# NEA 1907-2011
datNEA11 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE)
	setnames(datNEA11, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))
	outfileNEA11 <- 'analysis/LOF_07_to_LOF_S_11.wfabc'
	genNEA11=11 # for 1907 vs. 2011. sample sizes

	# CAN
datCAN <- fread('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', header=TRUE)
	setnames(datCAN, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))
	outfileCAN <- 'analysis/Can_40_to_Can.wfabc'
	genCAN=8 # for 1940 to contemporary. Guess 8 generations

# read in 25kmer filter
#loc25NEA <- fread('data_2017.11.24/Norway_25K_mer_positions.txt')
#loc25CAN <- fread('data_2017.11.24/Canada_25K_mer_positions.txt')
#
## trim to loci that meet 25kmer filter
#setkey(datNEA, CHROM, POS)
#setkey(datNEA11, CHROM, POS)
#setkey(datCAN, CHROM, POS)
#setkey(loc25NEA, CHROM, POS)
#setkey(loc25CAN, CHROM, POS)
#	nrow(datNEA)
#	nrow(datNEA11)
#	nrow(datCAN)
#datNEA <- datNEA[loc25NEA, nomatch=0] # nomatch=0 so that non-matching rows are dropped
#datNEA11 <- datNEA11[loc25NEA, nomatch=0] # nomatch=0 so that non-matching rows are dropped
#datCAN <- datCAN[loc25CAN, nomatch=0] # nomatch=0 so that non-matching rows are dropped
#	nrow(datNEA)
#	nrow(datNEA11)
#	nrow(datCAN)
	
# trim out inversions and Unplaced
datNEA <- datNEA[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
datNEA11 <- datNEA11[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
datCAN <- datCAN[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

	dim(datNEA)
	dim(datNEA11)
	dim(datCAN)
	
# figure out number of sites
nsitesNEA <- nrow(datNEA)
nsitesNEA11 <- nrow(datNEA11)
nsitesCAN <- nrow(datCAN)

# figure out # chromosomes in each sample of the focal allele
datNEA[,alcnt1:=round(N_CHR_1*Freq_1)]
datNEA[,alcnt2:=round(N_CHR_2*Freq_2)]
datNEA11[,alcnt1:=round(N_CHR_1*Freq_1)]
datNEA11[,alcnt2:=round(N_CHR_2*Freq_2)]
datCAN[,alcnt1:=round(N_CHR_1*Freq_1)]
datCAN[,alcnt2:=round(N_CHR_2*Freq_2)]

# make a data table to write out
sampsizesNEA <- datNEA[,.(N_CHR_1, N_CHR_2)] # sample sizes
sampsizesNEA[,index:=seq(1,length.out=nrow(datNEA),by=2)] # index to allow sorting with allele counts (interleaving)
alcntsNEA <- datNEA[,.(alcnt1, alcnt2)] # allele counts
alcntsNEA[,index:=seq(2,length.out=nrow(datNEA),by=2)] # index, offset by 1 compared to sample sizes
setnames(sampsizesNEA, 1:2, c('X1', 'X2')) # set column names to allow rbinding
setnames(alcntsNEA, 1:2, c('X1', 'X2'))
outNEA <- rbind(sampsizesNEA, alcntsNEA)
setkey(outNEA, index) # sort by index: each pair of rows is now sample sizes followed by allele counts

sampsizesNEA11 <- datNEA11[,.(N_CHR_1, N_CHR_2)] # sample sizes
sampsizesNEA11[,index:=seq(1,length.out=nrow(datNEA11),by=2)] # index to allow sorting with allele counts (interleaving)
alcntsNEA11 <- datNEA11[,.(alcnt1, alcnt2)] # allele counts
alcntsNEA11[,index:=seq(2,length.out=nrow(datNEA11),by=2)] # index, offset by 1 compared to sample sizes
setnames(sampsizesNEA11, 1:2, c('X1', 'X2')) # set column names to allow rbinding
setnames(alcntsNEA11, 1:2, c('X1', 'X2'))
outNEA11 <- rbind(sampsizesNEA11, alcntsNEA11)
setkey(outNEA11, index) # sort by index: each pair of rows is now sample sizes followed by allele counts

sampsizesCAN <- datCAN[,.(N_CHR_1, N_CHR_2)] # sample sizes
sampsizesCAN[,index:=seq(1,length.out=nrow(datCAN),by=2)] # index to allow sorting with allele counts (interleaving)
alcntsCAN <- datCAN[,.(alcnt1, alcnt2)] # allele counts
alcntsCAN[,index:=seq(2,length.out=nrow(datCAN),by=2)] # index, offset by 1 compared to sample sizes
setnames(sampsizesCAN, 1:2, c('X1', 'X2')) # set column names to allow rbinding
setnames(alcntsCAN, 1:2, c('X1', 'X2'))
outCAN <- rbind(sampsizesCAN, alcntsCAN)
setkey(outCAN, index) # sort by index: each pair of rows is now sample sizes followed by allele counts

# write out
cat(paste(nsitesNEA, ' ', ntimes, '\n', sep=''), file=outfileNEA) # header line 1
cat(paste('0,',genNEA, '\n', sep=''), file=outfileNEA, append=TRUE) # write out generation of each sample
write.table(outNEA[,.(X1, X2)], file=outfileNEA, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',') # append the rows of data

cat(paste(nsitesNEA11, ' ', ntimes, '\n', sep=''), file=outfileNEA11) # header line 1
cat(paste('0,',genNEA11, '\n', sep=''), file=outfileNEA11, append=TRUE) # write out generation of each sample
write.table(outNEA11[,.(X1, X2)], file=outfileNEA11, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',') # append the rows of data

cat(paste(nsitesCAN, ' ', ntimes, '\n', sep=''), file=outfileCAN) # header line 1
cat(paste('0,',genCAN, '\n', sep=''), file=outfileCAN, append=TRUE) # write out generation of each sample
write.table(outCAN[,.(X1, X2)], file=outfileCAN, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',') # append the rows of data
