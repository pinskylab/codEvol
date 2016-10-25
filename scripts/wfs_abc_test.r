# Testing the ABC procedue

# load functions
source('scripts/wfs.r')
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel)
	require(data.table)
	require(abc)
	require(boa)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(locfit, lib.loc="/projects/cees/lib/R_packages/") # needed by abc
	require(abc.data, lib.loc="/projects/cees/lib/R_packages/") # needed by abc
	require(abc, lib.loc="/projects/cees/lib/R_packages/")
}


# load data
load('analysis/wfs_sims_notrim_fullsampsize.rdata') # load the simulations: 'out' array, trimmed to only sims with full sample size

# turn into a data.frame for easy access
sims <- as.data.frame(t(out))
	rm(out)
	sims$diff <- sims$f2samp - sims$f1samp

# plots
ss <- sample(1:nrow(sims), 10000)
colrmp <- colorRamp(colors=brewer.pal(11, 'RdYlBu'))
getcol <- function(x, alpha=1){ rgbs <-colrmp(x); return(rgb(rgbs[,1], rgbs[,2], rgbs[,3], alpha*256, maxColorValue=256))}

plot(sims$f1samp[ss], sims$f1[ss], col=rgb(0,0,0,0.2))

plot(sims$f1samp[ss]-sims$f2samp[ss], sims$s[ss], col=rgb(0,0,0,0.2))

quartz(width=7, height=6)
layout(matrix(c(1,2),ncol=2),widths=c(4,1))
par(mai=c(1,1,0.2, 0.2), las=1, mgp=c(2, 0.6, 0), tcl=-0.2)
with(sims[ss,], plot(f1samp, f2samp-f1samp, xlab='Observed initial frequency', ylab='Observed change in frequency', col=getcol(s/2+0.5, 0.2))) # background
	abline(h=0, col='grey')
plot(rep(1,100),seq(-1,1,length.out=100), col=getcol(seq(-1,1,length.out=100)/2+0.5), bty='n', xaxt='n', ylab='Posterior mean s', xlab='', pch=15, cex=0.7) # legend

# sims with low initial frequency
inds <- sims$f1samp < 0.01
sum(inds)
with(sims[inds,], plot(s, f1samp-f2samp, main='Initial freq < 0.01', col=rgb(0,0,0,0.01)))

plot(density(sims$s[sims$f1samp < 0.05 & abs(sims$f2samp - sims$f1samp) < 0.05]))
boa.hpd(sims$s[sims$f1samp < 0.05 & abs(sims$f2samp - sims$f1samp) < 0.05], alpha=0.05)

# principal components analysis
pc <- prcomp(sims[1:1000,c('f1samp', 'f2samp', 'diff', 'fsiprime', 'fsdprime', 's', 'f1', 'ne')], scale.=TRUE, center=TRUE)
biplot(pc)