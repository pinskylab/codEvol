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
#trim <- 'LG03'; outfile <- 'analysis/LOF_07_LG03_to_LOF_S_14_LG03.wfabc'
trim <- ''; outfile <- 'analysis/LOF_07_to_LOF_S_14.wfabc'

# read in data (choose one)
dat <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); nm='1907-2014'; gen=11 # for 1907 vs. 2014. sample sizes


# trim to one linkage group?
if(trim != ''){
	dat <- dat[dat$CHROM==trim,]
	print(paste('trimmed to', trim))
}
	dim(dat)
	
# trim out inversions and Unplaced
dat <- dat[!(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')),]
#dat <- dat[!(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

	dim(dat)
	
# figure out number of sites
nsites <- nrow(dat)

# figure out # chromosomes in each sample of the focal allele
dat[,alcnt1:=round(N_CHR_1*Freq_1)]
dat[,alcnt2:=round(N_CHR_2*Freq_2)]

# make a data table to write out
sampsizes <- dat[,.(N_CHR_1, N_CHR_2)] # sample sizes
sampsizes[,index:=seq(1,length.out=nrow(dat),by=2)] # index to allow sorting with allele counts (interleaving)
alcnts <- dat[,.(alcnt1, alcnt2)] # allele counts
alcnts[,index:=seq(2,length.out=nrow(dat),by=2)] # index, offset by 1 compared to sample sizes
setnames(sampsizes, 1:2, c('X1', 'X2')) # set column names to allow rbinding
setnames(alcnts, 1:2, c('X1', 'X2'))
out <- rbind(sampsizes, alcnts)
setkey(out, index) # sort by index: each pair of rows is now sample sizes followed by allele counts

# write out
cat(paste(nsites, ' ', ntimes, '\n', sep=''), file=outfile) # header line 1
cat(paste('0,',gen, '\n', sep=''), file=outfile, append=TRUE) # write out generation of each sample
write.table(out[,.(X1, X2)], file=outfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',') # append the rows of data
