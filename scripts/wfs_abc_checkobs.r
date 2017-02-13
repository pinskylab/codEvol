# check observations


################################
# load functions and prep data
################################
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(RColorBrewer)
	require(data.table)
	require(hexbin)
	require(lattice)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(RColorBrewer, lib.loc="/projects/cees/lib/R_packages/")
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

# load data
dat <- fread('analysis/LOF_07_to_LOF_S_14.wfabc', skip=2) # the data on samples sizes and observed allele frequencies, from WFABC input file
locnms <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star


# prep data
	# extract sample frequencies
sampsize <- dat[seq(1,nrow(dat),by=2),] # sample sizes in # of chromosomes
obsfrqs <- dat[seq(2,nrow(dat),by=2),] # allele counts in # chromosomes
rm(dat)
setnames(obsfrqs, 1:2, c('count1', 'count2'))
setnames(sampsize, 1:2, c('size1', 'size2'))
obsfrqs <- cbind(obsfrqs, sampsize)
obsfrqs[,f1:=count1/size1]
obsfrqs[,f2:=count2/size2]
obsfrqs[,diff:=f2-f1]
obsfrqs[,locusnum:=1:nrow(obsfrqs)] # add locusnumber

	# trim locus names to match rest of data
#locnms <- locnms[CHROM=='LG03',] # if only looking at LG03
locnms <- locnms[!(locnms$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')),] # trim out inversions
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting

##############
## plots
##############
# fsi and fsd vs. f1samp and absdiff
i <- sample(1:nrow(obsfrqs), 2000)

pdf('analysis/figures/wfs_abc_obs_fs_vs_f1samp_diff.pdf')
par(mfrow=c(2,2))
targ[i,plot(f1samp, fsi)]
plot(obsfrqs$diff[i], targ$fsi[i], xlab='Diff', ylab='fsi')
targ[i,plot(f1samp, fsd)]
plot(obsfrqs$diff[i], targ$fsd[i], xlab='Diff', ylab='fsd')
dev.off()

	# also with lattice plot to do slices
myshingle <- equal.count(obsfrqs$f1[i], number=5, overlap=0.1)
xyplot(targ$fsi[i] ~ obsfrqs$diff[i]| myshingle)