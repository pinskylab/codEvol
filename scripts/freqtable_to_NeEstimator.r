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
outfile <- 'analysis/LOF_07_to_LOF_S_14.gen'

# read in data (choose one)
dat <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=56; c2=48; nm='1907-2014'; gen=11 # for 1907 vs. 2014. sample sizes


	
# trim out inversions and Unplaced
dat <- dat[!(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),]

	dim(dat)
	
# figure out number of sites
nsites <- nrow(dat)

# figure out # chromosomes in each sample of the focal allele
dat[,alcnt1:=round(N_CHR_1*Freq_1)]
dat[,alcnt2:=round(N_CHR_2*Freq_2)]

# calculate # homs, hets, and missing (HW proportions)
maxn1 <- dat[,max(N_CHR_1)/2] # num individuals at time 1
maxn2 <- dat[,max(N_CHR_2)/2] # num at time 2

dat[,nhom1a := round((N_CHR_1*Freq_1^2)/2)] # number of hom indivs
dat[,nhom1b := round(N_CHR_1*(1-Freq_1)^2/2)] # number of hom individuals (other allele)
dat[,nhet1 := round(N_CHR_1*2*Freq_1*(1-Freq_1)/2)] # number of het indivs
dat[,nmiss1 := maxn1-N_CHR_1/2] # number of missing individuals
dat[nhet1+nhom1a+nhom1b != (N_CHR_1/2),nhet1:=nhet1 + (N_CHR_1/2) - (nhet1+nhom1a+nhom1b)] # add missing individuals to hets, where missing

dat[,nhom2a := round((N_CHR_2*Freq_2^2)/2)] # number of hom indivs
dat[,nhom2b := round(N_CHR_2*(1-Freq_2)^2/2)] # number of hom individuals (other allele)
dat[,nhet2 := round(N_CHR_2*2*Freq_2*(1-Freq_2)/2)] # number of het indivs
dat[nhet2+nhom2a+nhom2b != (N_CHR_2/2),nhet2:=nhet2 + (N_CHR_2/2) - (nhet2+nhom2a+nhom2b)] # add missing individuals to hets, where missing
dat[,nmiss2 := maxn2-N_CHR_2/2] # number of missing individuals

	# check that numbers add up
	dat[nhet1+nhom1a+nhom1b < (N_CHR_1/2),] # any still missing?
	dat[nhet1+nhom1a+nhom1b > (N_CHR_1/2),] # any still excess?

	dat[nhet2+nhom2a+nhom2b < (N_CHR_2/2),] # any still missing?
	dat[nhet2+nhom2a+nhom2b > (N_CHR_2/2),] # any still excess?

	dat[nhet1+nhom1a+nhom1b+nmiss1 != maxn1,]
	dat[nhet2+nhom2a+nhom2b+nmiss2 != maxn2,]

# make a data table to write out
out1 <- data.table(pop=rep('1,', maxn1), l1 = c(
	rep('00', dat[1,nmiss1]), # missing data
	rep('11', dat[1,nhom1a]), # homozygotes
	rep('12', dat[1,nhet1]), # heterozygotes
	rep('22', dat[1,nhom1b]) # homozygotes
	))
for(i in 2:nrow(dat)){
	if(i %% 1000 == 0) print(i)
	colnm <- paste('l', i, sep='')
	out1[,(colnm):=c(
	rep('00', dat[i,nmiss1]), # missing data
	rep('11', dat[i,nhom1a]), # homozygotes
	rep('12', dat[i,nhet1]), # heterozygotes
	rep('22', dat[i,nhom1b]))] # homozygotes
}

out2 <- data.table(pop=rep('2,', maxn2), l1 = c(
	rep('00', dat[1,nmiss2]), # missing data
	rep('11', dat[1,nhom2a]), # homozygotes
	rep('12', dat[1,nhet2]), # heterozygotes
	rep('22', dat[1,nhom2b]) # homozygotes
	))
for(i in 2:nrow(dat)){
	if(i %% 1000 == 0) print(i)
	colnm <- paste('l', i, sep='')
	out2[,(colnm):=c(
	rep('00', dat[i,nmiss2]), # missing data
	rep('11', dat[i,nhom2a]), # homozygotes
	rep('12', dat[i,nhet2]), # heterozygotes
	rep('22', dat[i,nhom2b]))] # homozygotes
}


# write out
cat(paste('NEA cod sampled', nm, 'without LG01, LG02, LG07, LG12 or Unplaced scaffolds\n'), file=outfile) # header line 1
cat(paste(1:nrow(dat), collapse=','), file=outfile, append=TRUE) # write out locus names on one line
cat('\npop\n', file=outfile, append=TRUE)
write.table(out1, file=outfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t') # append the rows of data
cat('pop\n', file=outfile, append=TRUE)
write.table(out2, file=outfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t') # append the rows of data
