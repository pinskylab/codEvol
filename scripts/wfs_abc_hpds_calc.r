# run after wfs_abc.r or wfs_abc_function.r (latter called by wfs_abc_sbatch.sh or wfs_abc_sbatch_3times.sh, which in turn are called by script output by wfs_abc_makesbatchscript.r)

# parameters
#pop <- 'Lof'; myyr1 <- '07'; myyr2 <- '1114' # all 3 time points
pop <- 'Can'; myyr1 <- '40'; myyr2 <- '13'

################################
# load functions and prep data
################################
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){ # on something other than cod node
	require(RColorBrewer)
	require(data.table)
	require(boa)
	require(hexbin)
	require(plyr) # for summarize
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){ # on cod node
	require(RColorBrewer, lib.loc="/projects/cees/lib/R_packages/")
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(boa, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
}

# mean, CIs for abc posteriors
# assumes data in column 2
mci <- function(x){
	b <- as.numeric(boa.hpd(x[,2], alpha=0.05))
	return(c(mean=mean(x[,2]), l95=b[1], u95=b[2]))
}

# mean, CIs, and p(x<>0) for abc posteriors
# assumes data in column 2
mcip <- function(x){ 
	b <- as.numeric(boa.hpd(x[,2], alpha=0.05))
	p <- sum(x[,2]>0)/nrow(x) # fraction of tail above 0
	if(p>0.5) p <- 1-p # convert to smaller of the two tails
	p <- 2*p # two-tailed test
	return(c(mean=mean(x[,2]), l95=b[1], u95=b[2], p=p))
}

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

####################
# load data
####################
if(pop == 'Lof'){
	if(myyr2 != '1114'){
		targfile <- paste('data_2017.11.24/Frequency_table_Lof', myyr1, '_Lof', myyr2, '.txt', sep='')
	}
	if(myyr2 == '1114'){
		targfile <- paste('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', sep='')
		targfile2 <- paste('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', sep='')
	}
}
if(pop == 'Can'){
	targfile <- paste('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', sep='')
}
hpds <- fread(targfile, header=TRUE)
setnames(hpds, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
hpds[,locusnum:=1:nrow(hpds)] # add a locus number indicator

if(myyr2!='1114'){
	hpds <- hpds[,.(locusnum, CHROM, POS, f1samp, f2samp)]
}
if(myyr2=='1114'){
	hpds2 <- fread(targfile2, header=TRUE)
	setnames(hpds2, 3:7, c('alcnt1', 'f1samp', 'alcnt3', 'f3samp', 'ABS_DIFF13'))
	setkey(hpds, CHROM, POS)
	setkey(hpds2, CHROM, POS)
	hpds <- hpds[hpds2,.(locusnum, CHROM, POS, f1samp, f2samp, f3samp)]
}


# trim to focal loci (outliers)
locs <- fread('analysis/outlier_annotation.csv') # the outliers
print(paste('Started with', nrow(hpds), 'loci'))
if(pop == 'Lof'){
	locs <- locs[q3.Lof071114 !='' | q3.comb071114Can !='',.(CHROM, POS)]
	locs[, POS := as.numeric(POS)]
	setkey(locs, CHROM, POS)
	setkey(hpds, CHROM, POS)
	hpds <- merge(hpds, locs)
}
if(pop == 'Can'){
	locs <- locs[q3.Can !='' | q3.comb071114Can !='',.(CHROM, POS)]
	locs[, POS := as.numeric(POS)]
	setkey(locs, CHROM, POS)
	setkey(hpds, CHROM, POS)
	hpds <- merge(hpds, locs)
}
print(paste('Trimmed to', nrow(hpds), 'loci by using only outliers in analysis/outlier_annotation.csv'))

# prep data
a <- rep(as.numeric(NA), nrow(hpds))
hpds[, ':=' (f1samp.mean=a, f1samp.l95=a, f1samp.u95=a, f2samp.mean=a, f2samp.l95=a, f2samp.u95=a, ne.mean=a, ne.l95=a, ne.u95=a, f1.mean=a, f1.l95=a, f1.u95=a, f2.mean=a, f2.l95=a, f2.u95=a, s.mean=a, s.l95=a, s.u95=a, s.p=a)]




#########################################################
# calculate HPDs
# loop over all the posterior files output by wfs_abc.r
#########################################################

# find files to read in
files <- list.files(path='analysis/temp', pattern='wfs_abc_sampsize*', full.names=TRUE)
print(paste(length(files), 'to process'))

# read in each file and calculate HPD: ~48 hours on cod node
for(i in 1:length(files)){
	print(paste(i, 'of', length(files), ':', files[i]))
	posts <- read.csv(gzfile(files[i]))
	theselocs <- sort(unique(posts$locus))
	inds <- which(hpds$locusnum %in% theselocs)
	
	if(length(inds)>0){
		temp <- ddply(.data=posts[,c('locus', 'f1samp')], .variables= ~locus, .fun=mci)
		hpds[inds, ':=' (f1samp.mean=temp$mean, f1samp.l95=temp$l95, f1samp.u95=temp$u95)]

		temp <- ddply(.data=posts[,c('locus', 'f2samp')], .variables= ~locus, .fun=mci)
		hpds[inds, ':=' (f2samp.mean=temp$mean, f2samp.l95=temp$l95, f2samp.u95=temp$u95)]

		temp <- ddply(.data=posts[,c('locus', 'ne')], .variables= ~locus, .fun=mci)
		hpds[inds, ':=' (ne.mean=temp$mean, ne.l95=temp$l95, ne.u95=temp$u95)]

		temp <- ddply(.data=posts[,c('locus', 'f1')], .variables= ~locus, .fun=mci)
		hpds[inds, ':=' (f1.mean=temp$mean, f1.l95=temp$l95, f1.u95=temp$u95)]

		temp <- ddply(.data=posts[,c('locus', 'f2')], .variables= ~locus, .fun=mci)
		hpds[inds, ':=' (f2.mean=temp$mean, f2.l95=temp$l95, f2.u95=temp$u95)]

		temp <- ddply(.data=posts[,c('locus', 's')], .variables= ~locus, .fun=mcip)
		hpds[inds, ':=' (s.mean=temp$mean, s.l95=temp$l95, s.u95=temp$u95, s.p=temp$p)]
	} else {
		warning(paste('WARNING: Could not find locusnum(s)', paste(theselocs, collapse=', ')))
	}
}


# write out
outfile <- paste('analysis/wfs_abc_hpds_', pop, '.csv', sep='')
outfile
write.csv(hpds, file=outfile)



