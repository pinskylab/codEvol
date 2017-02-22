# run after wfs_abc.r
# Look at the posteriors for one locus

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

source('scripts/wfs_simp.r')

# mean, CIs for abc posteriors
# assumes data in column 2
mci <- function(x){
	b <- as.numeric(boa.hpd(x[,2], alpha=0.05))
	return(c(mean=mean(x[,2]), l95=b[1], u95=b[2]))
}

# mean, CIs, and min/max for abc posteriors
# assumes data in column 2
mcimm <- function(x){ 
	b <- as.numeric(boa.hpd(x[,2], alpha=0.05))
	mn <- min(x[,2])
	mx <- max(x[.2])
	return(c(mean=mean(x[,2]), l95=b[1], u95=b[2], min=mn, max=mx))
}

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

# load data
dat <- locnms <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_abc_hpds.rdata') # the HPD calculations

# trim data
dat <- dat[!(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')),] # trim out inversions
dat[,dat:=1:nrow(dat)] # add a locusnumber for plotting

# check
nrow(dat) == nrow(hpds$s)
nrow(dat) == nrow(hpds$f1samp)


#########################################################
# Examine posteriors
#########################################################
# find a locus
inds <- which(dat[,Freq_2-Freq_1 > 0.4 & N_CHR_1>=44 & N_CHR_2 >= 44] & hpds$s$mean > 0.5); length(inds)
cbind(dat[inds,], hpds$s[inds,], hpds$f1samp$mean[inds], hpds$f2samp$mean[inds])

inds <- which(dat[,Freq_2-Freq_1 < 0.05] & hpds$s$mean > 0.5); length(inds)
cbind(dat[inds,], hpds$s[inds,], hpds$f1samp$mean[inds], hpds$f2samp$mean[inds])


# find files to read in
myloc <- 57235; infile='analysis/temp/wfs_abc_sampsize44,46_locus055998-058015.csv.gz' # locus 57235, LG06 4709919: start at 0.55, go to 0.96
myloc <- 73; infile='analysis/temp/wfs_abc_sampsize50,44_locus000043-004882.csv.gz'; c1=50;c2=44 # locus 73, LG03 262521: start at 0.98, go to 1

# read in file and extract posterior samples
posts <- read.csv(gzfile(infile))
posts <- posts[posts$locus == myloc,]
	dim(posts)

# simulate with s=0
	# simulations
	nes <- sample(posts$ne, 1000)
	gen=11
	fs <- sample(posts$f1, 1000)

	sims <- data.frame(f1=fs)
	sims$f2 <- sims$f1samp <- sims$f2samp <- NA
	for(j in 1:length(fs)){
		sims[j,c('f2', 'f1samp', 'f2samp')] <- wfs_simp(j, f1=fs[j], s=0, c1=c1, c2=c2, gen=gen, ne=nes[j])
		
	}

	
# histograms
xlimdiff <- range(c(posts$f2 - posts$f1, posts$f2samp - posts$f1samp, sims$f2-sims$f1, sims$f2samp - sims$f1samp))
xlimdiff[1] <- floor(xlimdiff[1]*100)/100
xlimdiff[2] <- ceiling(xlimdiff[2]*100)/100

quartz(width=8, height=6)
# png(width=8, height=6, filename=paste('analysis/figures/wfs_abc_posteriors_locus', dat[myloc,CHROM], '_', dat[myloc,POS], '.png', sep=''), res=300, units='in')
par(mfrow=c(3,3), omi=c(0,0,0.4,0), mai=c(0.6, 0.6, 0.3, 0.1))

hist(posts$f1samp, col='grey', breaks=seq(0,1,by=0.01), main='Initial sample freq\nwith obs', xlab='Frequency', freq=FALSE)
abline(v=dat[myloc,Freq_1], lty=1, col='red')

hist(posts$f2samp, col='grey', breaks=seq(0,1,by=0.01), main='Final sample freq\nwith obs', xlab='Frequency', freq=FALSE)
abline(v=dat[myloc,Freq_2], lty=1, col='red')

hist(posts$f1, col='grey', breaks=seq(0,1,by=0.01), main='Initial pop freq\nwith obs', xlab='Frequency', freq=FALSE)
abline(v=dat[myloc,Freq_1], lty=1, col='red')

hist(posts$f2, col='grey', breaks=seq(0,1,by=0.01), main='Final pop freq\nwith obs', xlab='Frequency', freq=FALSE)
abline(v=dat[myloc,Freq_2], lty=1, col='red')

bks <- seq(xlimdiff[1],xlimdiff[2],by=0.01)
hst <- hist(posts$f2 - posts$f1, breaks=bks, plot=FALSE)
hst2 <- hist(sims$f2-sims$f1, plot=FALSE, breaks=bks) # histogram for the null model change in pop frequency
ylims <- c(0, max(c(hst$density, hst2$density)))
hist(posts$f2 - posts$f1, col='grey', breaks=bks, main='Change in pop freq\nwith obs & null', xlab='Frequency', freq=FALSE, xlim=xlimdiff, ylim=ylims)
abline(v=dat[myloc,Freq_2 - Freq_1], lty=1, col='red')

	# overplot the null model (s=0)
	lines(hst2$mids, hst2$density, type='l', pch=16, col='green') # change in pop freq

hst3 <- hist(posts$f2samp - posts$f1samp, breaks=bks, plot=FALSE)
hst4 <- hist(sims$f2samp-sims$f1samp, plot=FALSE, breaks=bks) # histogram for the null model change in sample frequency
ylims <- c(0, max(c(hst3$density, hst4$density)))
hist(posts$f2samp - posts$f1samp, col='grey', breaks=bks, main='Change in sample freq\nwith obs & null', xlab='Frequency', freq=FALSE, xlim=xlimdiff, ylim=ylims)
abline(v=dat[myloc,Freq_2 - Freq_1], lty=1, col='red')

	# overplot the null model (s=0)
	lines(hst4$mids, hst4$density, type='l', pch=16, col='green') # change in sample freq


hist(posts$s, col='grey', breaks=seq(-1,1,by=0.05), main='Inferred s\nwith min,max', xlab='s')
abline(v=range(posts$s), lty=1, col='red')

mtext(paste('Locus ', dat[myloc,CHROM], ' ', dat[myloc,POS], ' (#', myloc, '), ', nrow(posts), ' sims', sep=''), side=3, outer=TRUE, line=1)

dev.off()