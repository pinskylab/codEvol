# run after wfs_abc.r

################################
# load functions and prep data
################################
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(RColorBrewer)
	require(data.table)
	require(boa)
	require(hexbin)
	require(plyr) # for summarize
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
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

# load data
dat <- fread('analysis/LOF_07_to_LOF_S_14.wfabc', skip=2) # the data on samples sizes and observed allele frequencies, from WFABC input file. for now only used to set up hpds dataframes


# prep data
obsfrqs <- dat[seq(2,nrow(dat),by=2),] # allele counts in # chromosomes. for now, only used to set up hpds data.frames.


#########################################################
# calculate HPDs
# loop over all the posterior files output by wfs_abc.r
#########################################################
startind <- 1 # to start looping over files part-way through (if desired)

# read in hpds if the file exists
if(file.exists('analysis/wfs_abc_hpds.rdata')){
	load('analysis/wfs_abc_hpds.rdata')
} else {
	# set up list to hold hpds (if starting from scratch)
	hpds <- list()
	a <- rep(NA, nrow(obsfrqs))
	locs <- 1:nrow(obsfrqs)
	hpds$f1samp <- data.frame(locus=locs, mean=a, l95=a, u95=a)
	hpds$f2samp <- data.frame(locus=locs, mean=a, l95=a, u95=a)
	hpds$ne <- data.frame(locus=locs, mean=a, l95=a, u95=a)
	hpds$f1 <- data.frame(locus=locs, mean=a, l95=a, u95=a)
	hpds$f2 <- data.frame(locus=locs, mean=a, l95=a, u95=a)
	hpds$s <- data.frame(locus=locs, mean=a, l95=a, u95=a, p=a)
}

# find files to read in
files <- list.files(path='analysis/temp', pattern='wfs_abc_sampsize*', full.names=TRUE)
print(paste(length(files), 'to process'))

# read in each file and calculate HPD: ~48 hours on cod node
for(i in startind:length(files)){
	print(paste(i, 'of', length(files), ':', files[i]))
	posts <- read.csv(gzfile(files[i]))
	theselocs <- sort(unique(posts$locus))
	inds <- hpds$f1samp$locus %in% theselocs
	hpds$f1samp[inds,] <- ddply(.data=posts[,c('locus', 'f1samp')], .variables= ~locus, .fun=mci)
	hpds$f2samp[inds,] <- ddply(.data=posts[,c('locus', 'f2samp')], .variables= ~locus, .fun=mci)
	hpds$ne[inds,] <- ddply(.data=posts[,c('locus', 'ne')], .variables= ~locus, .fun=mci)
	hpds$f1[inds,] <- ddply(.data=posts[,c('locus', 'f1')], .variables= ~locus, .fun=mci)
	hpds$f2[inds,] <- ddply(.data=posts[,c('locus', 'f2')], .variables= ~locus, .fun=mci)
	hpds$s[inds,] <- ddply(.data=posts[,c('locus', 's')], .variables= ~locus, .fun=mcip)
}


# write out
save(hpds, file='analysis/wfs_abc_hpds.rdata')



