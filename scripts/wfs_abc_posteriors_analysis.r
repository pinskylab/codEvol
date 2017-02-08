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



#########################################################
# calculate HPDs
# loop over all the posterior files output by wfs_abc.r
#########################################################

# set up list to hold hpds
hpds <- list()
a <- rep(NA, nrow(obsfrqs))
locs <- 1:nrow(obsfrqs)
hpds$f1samp <- data.frame(locus=locs, mean=a, l95=a, u95=a)
hpds$fsd <- data.frame(locus=locs, mean=a, l95=a, u95=a)
hpds$fsi <- data.frame(locus=locs, mean=a, l95=a, u95=a)
hpds$ne <- data.frame(locus=locs, mean=a, l95=a, u95=a)
hpds$f1 <- data.frame(locus=locs, mean=a, l95=a, u95=a)
hpds$s <- data.frame(locus=locs, mean=a, l95=a, u95=a, p=a)

# find files to read in
files <- list.files(path='analysis/temp', pattern='wfs_abc_sampsize*', full.names=TRUE)
print(paste(length(files), 'to process'))
for(i in 1:length(files)){
	print(files[i])
	posts <- read.csv(gzfile(files[i]))
	theselocs <- sort(unique(posts$locus))
	inds <- hpds$f1samp$locus %in% theselocs
	hpds$f1samp[inds,] <- ddply(.data=posts[,c('locus', 'f1samp')], .variables= ~locus, .fun=mci)
	hpds$fsd[inds,] <- ddply(.data=posts[,c('locus', 'fsdprime')], .variables= ~locus, .fun=mci)
	hpds$fsi[inds,] <- ddply(.data=posts[,c('locus', 'fsiprime')], .variables= ~locus, .fun=mci)
	hpds$ne[inds,] <- ddply(.data=posts[,c('locus', 'ne')], .variables= ~locus, .fun=mci)
	hpds$f1[inds,] <- ddply(.data=posts[,c('locus', 'f1')], .variables= ~locus, .fun=mci)
	hpds$s[inds,] <- ddply(.data=posts[,c('locus', 's')], .variables= ~locus, .fun=mcip)
}

	# all FDR correction for p-values
hpds$s$p.adj <- p.adjust(hpds$s$p, method='fdr')


# find loci potentially under selection
posinds <- hpds$s$p.adj < 0.05 & hpds$s$mean > 0 & !is.na(hpds$s$mean)
neginds <- hpds$s$p.adj < 0.05 & hpds$s$mean < 0 & !is.na(hpds$s$mean)

sum(posinds) # number of loci
sum(neginds)

# Examine candidates
print(locnms[posinds,], nrow=sum(posinds))
print(locnms[neginds,], nrow=sum(neginds))

hpds$s[posinds,]

# write out
out <- hpds$s
save(out, file='analysis/wfs_abc_hpds_s.rdata')

################################
# plots
################################
locnms <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star

load('analysis/wfs_abc_hpds_s.rdata')
hpds <- vector('list', 6)
hpds$s <- out

colrmp <- colorRamp(colors=brewer.pal(11, 'RdYlBu'))
getcol <- function(x, alpha=1){ rgbs <-colrmp(x); return(rgb(rgbs[,1], rgbs[,2], rgbs[,3], alpha*256, maxColorValue=256))}

# Initial frequency vs. posterior s
plot(hpds$f1$mean, hpds$s$mean, xlab='Initial frequency', ylab='Posterior mean s', col=getcol(hpds$s$mean/2+0.5))
points(hpds$f1$mean[neginds|posinds], hpds$s$mean[neginds|posinds], xlab='Initial frequency', ylab='Posterior mean s', pch=16)

# Obs init freq vs. posterior init freq
plot(obsfrqs$f1, hpds$f1$mean, col=getcol(hpds$s$mean/2+0.5, 0.2), xlab='Observed initial frequency', ylab='Posterior mean initial frequency')
points(obsfrqs$f1[neginds|posinds], hpds$f1$mean[neginds|posinds], pch=16)

# Initial frequency vs. change in frequency vs. posterior s
# overplot loci "under selection"
inds <- neginds | posinds
	
quartz(width=7, height=6)
layout(matrix(c(1,2),ncol=2),widths=c(4,1))
par(mai=c(1,1,0.2, 0.2), las=1, mgp=c(2.4, 0.6, 0), tcl=-0.2)
plot(obsfrqs$f1, obsfrqs$diff, xlab='Observed initial frequency', ylab='Observed change in frequency', col=getcol(hpds$s$mean/2+0.5, 0.2), ylim=c(-1,1)) # background
points(obsfrqs$f1[inds], obsfrqs$diff[inds], col='black', pch=16, cex=0.8)
plot(rep(1,100),seq(-1,1,length.out=100), col=getcol(seq(-1,1,length.out=100)/2+0.5), bty='n', xaxt='n', ylab='Posterior mean s', xlab='', pch=15, cex=0.7) # legend


# Examine MCMC sample distributions for individual loci
i <- which(obsfrqs$f1 > 0.8 & obsfrqs$f1 < 0.90 & obsfrqs$diff < -0.25) # initially high frequency loci under negative selection
j <- 2 # pick which one
plot(density(posts$s[i[j],]))
plot(density(posts$f1[i[j],])); abline

hexbinplot(posts$s[i[j],] ~ posts$f1[i[j],], colramp=rf, main=paste('locus', i[j]), xlab='posterior mean initial frequency', ylab='posterior mean s')

hexbinplot(posts$s[i[j],] ~ posts$ne[i[j],], colramp=rf, main=paste('locus', i[j]), xlab='posterior mean Ne', ylab='posterior mean s')


# Plot of s vs. genome position
i <- locnms$CHROM=='LG03'
plot(locnms$POS[i], hpds$s$mean,type='p', cex=0.5, ylim=c(-1,1))

plot(locnms$POS[i], hpds$s$mean,type='p', cex=0.5, xlim=c(4.3e6,4.5e6), ylim=c(-1,1))