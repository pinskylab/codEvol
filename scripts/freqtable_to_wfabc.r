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

# read in data (choose one)
datNEA <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE)
	setnames(datNEA, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))
	outfileNEA <- 'analysis/LOF_07_to_LOF_S_14.wfabc'
	genNEA=11 # for 1907 vs. 2014. sample sizes
datCAN <- fread('data_21_02_17/Frequency_table_Can_40_Can_TGA.txt', header=TRUE)
	setnames(datCAN, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF'))
	outfileCAN <- 'analysis/Can_40_to_Can_TGA.wfabc'
	genCAN=8 # for 1940 to contemporary. Guess 8 generations

	
# trim out inversions and Unplaced
#dat <- dat[!(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')),]
datNEA <- datNEA[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]
datCAN <- datCAN[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

	dim(datNEA)
	dim(datCAN)
	
# figure out number of sites
nsitesNEA <- nrow(datNEA)
nsitesCAN <- nrow(datCAN)

# figure out # chromosomes in each sample of the focal allele
datNEA[,alcnt1:=round(N_CHR_1*Freq_1)]
datNEA[,alcnt2:=round(N_CHR_2*Freq_2)]
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

cat(paste(nsitesCAN, ' ', ntimes, '\n', sep=''), file=outfileCAN) # header line 1
cat(paste('0,',genCAN, '\n', sep=''), file=outfileCAN, append=TRUE) # write out generation of each sample
write.table(outCAN[,.(X1, X2)], file=outfileCAN, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',') # append the rows of data
