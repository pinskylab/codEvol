# some basic checking of the Wright-Fisher with selection ABC simulations
# to run on a cod node

##############
## data prep
##############

if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
	require(ff)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(bit, lib.loc="/projects/cees/lib/R_packages/") # for use with ff
	require(ff, lib.loc="/projects/cees/lib/R_packages/") # for big objects shared across processes
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}

# the abc simulations
sampsize <- '44,44'
ffnm <- paste('analysis/temp/wfs_sims_ff', sampsize, sep='')
ffload(ffnm, overwrite=TRUE) # takes 30 sec or so. loads thisout.ff

head(thisout.ff)

# meansds for untransforming the sumstats
load('analysis/wfs_sims_meansds.rdata')
meansds <- meansds[alcnt1==44 & alcnt2==44,]



### basic plots
i <- 1:2000 # select a few simulations to simplify plotting

# fsi and fsd vs. diff and f1samp
pdf('analysis/figures/wfs_abc_sims_fs_vs_f1samp_diff.pdf')
par(mfrow=c(2,2))
plot(thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean, thisout.ff['fsiprime', i]*meansds$fsi.sd + meansds$fsi.mean, xlab='f1samp', ylab='fsi')
plot(thisout.ff['f2samp',i]-thisout.ff['f1samp',i]*meansds$f1samp.sd - meansds$f1samp.mean, thisout.ff['fsiprime', i]*meansds$fsi.sd + meansds$fsi.mean, xlab='diff', ylab='fsi')
plot(thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean, thisout.ff['fsdprime', i]*meansds$fsd.sd + meansds$fsd.mean, xlab='f1samp', ylab='fsd')
plot(thisout.ff['f2samp',i]-thisout.ff['f1samp',i]*meansds$f1samp.sd - meansds$f1samp.mean, thisout.ff['fsdprime', i]*meansds$fsd.sd + meansds$fsd.mean, xlab='diff', ylab='fsd')
dev.off()


# s vs. fsi and fsd and f1samp
pdf('analysis/figures/wfs_abc_sims_s_vs_fsi_fsd_f1samp.pdf', width=9, height=4)
par(mfrow=c(1,3))
plot(thisout.ff['fsiprime',i]*meansds$fsi.sd + meansds$fsi.mean, thisout.ff['s', i], xlab='fsi', ylab='s')
plot(thisout.ff['fsdprime',i]*meansds$fsd.sd + meansds$fsd.mean, thisout.ff['s', i], xlab='fsd', ylab='s')
plot(thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean, thisout.ff['s', i], xlab='f1samp', ylab='s')
dev.off()

# s vs. diff and f1samp
pdf('analysis/figures/wfs_abc_sims_s_vs_diff_f1samp.pdf', width=9, height=4)
par(mfrow=c(1,2))
plot(thisout.ff['f2samp',i]-thisout.ff['f1samp',i]*meansds$f1samp.sd - meansds$f1samp.mean, thisout.ff['s', i], xlab='diff', ylab='s')
plot(thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean, thisout.ff['s', i], xlab='f1samp', ylab='s')
dev.off()

	# in slices of f1samp
i <- 1:10000 # selection a few simulations to simplify plotting

pdf('analysis/figures/wfs_abc_sims_s_vs_obsdiff_by_f1samp.pdf', width=9, height=4)
par(mfrow=c(2,5))
ints <- list(c(0,0.1), c(0.1, 0.2), c(0.2, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.6), c(0.6, 0.7), c(0.7, 0.8), c(0.8, 0.9), c(0.9, 1))
for(j in 1:length(ints)){
	k <- which((thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean) > ints[[j]][1] & (thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean) < ints[[j]][2])
	plot(thisout.ff['f2samp',i][k]-thisout.ff['f1samp',i][k]*meansds$f1samp.sd - meansds$f1samp.mean, thisout.ff['s', i][k], xlab='diff', ylab='s', main=paste('f1samp', paste(ints[[j]], collapse='-')), xlim=c(-1,1), ylim=c(-1,1), cex=0.5)
}
dev.off()


# s vs. diff and f1
	# in slices of f1
i <- 1:10000 # selection a few simulations to simplify plotting

pdf('analysis/figures/wfs_abc_sims_s_vs_diff_by_f1.pdf', width=9, height=4)
par(mfrow=c(2,5))
ints <- list(c(0,0.1), c(0.1, 0.2), c(0.2, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.6), c(0.6, 0.7), c(0.7, 0.8), c(0.8, 0.9), c(0.9, 1))
for(j in 1:length(ints)){
	k <- which(thisout.ff['f1',i] > ints[[j]][1] & thisout.ff['f1',i] < ints[[j]][2])
	plot(thisout.ff['f2',i][k]-thisout.ff['f1',i][k], thisout.ff['s', i][k], xlab='diff', ylab='s', main=paste('f1', paste(ints[[j]], collapse='-')), xlim=c(-1,1), ylim=c(-1,1), cex=0.5)
}
dev.off()

# s vs. f2 in slices of f1
i <- 1:10000 # selection a few simulations to simplify plotting

pdf('analysis/figures/wfs_abc_sims_s_vs_f2_by_f1.pdf', width=9, height=4)
par(mfrow=c(2,5))
ints <- list(c(0,0.1), c(0.1, 0.2), c(0.2, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.6), c(0.6, 0.7), c(0.7, 0.8), c(0.8, 0.9), c(0.9, 1))
for(j in 1:length(ints)){
	k <- which(thisout.ff['f1',i] > ints[[j]][1] & thisout.ff['f1',i] < ints[[j]][2])
	plot(thisout.ff['f2',i][k], thisout.ff['s', i][k], xlab='f2', ylab='s', main=paste('f1', paste(ints[[j]], collapse='-')), xlim=c(0,1), ylim=c(-1,1), cex=0.5, col=rgb(0.2, 0.2, 0.2, 0.2))
}
dev.off()

# s vs. f2samp in slices of f1samp
i <- 1:10000 # selection a few simulations to simplify plotting

pdf('analysis/figures/wfs_abc_sims_s_vs_f2samp_by_f1samp.pdf', width=9, height=4)
par(mfrow=c(2,5))
ints <- list(c(0,0.1), c(0.1, 0.2), c(0.2, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.6), c(0.6, 0.7), c(0.7, 0.8), c(0.8, 0.9), c(0.9, 1))
f1samp <- thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean
f2samp <- thisout.ff['f2samp',i]
s <- thisout.ff['s',i]
for(j in 1:length(ints)){
	k <- which(f1samp > ints[[j]][1] & f1samp < ints[[j]][2])
	plot(f2samp[k], s[k], xlab='f2samp', ylab='s', main=paste('f1samp', paste(ints[[j]], collapse='-')), xlim=c(0,1), ylim=c(-1,1), cex=0.5, col=rgb(0.2, 0.2, 0.2, 0.2))
}
dev.off()



# s vs. fsi in slices of f1samp
i <- 1:10000 # selection a few simulations to simplify plotting

pdf('analysis/figures/wfs_abc_sims_s_vs_fsi_by_f1samp.pdf', width=9, height=4)
par(mfrow=c(2,5))
ints <- list(c(0,0.1), c(0.1, 0.2), c(0.2, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.6), c(0.6, 0.7), c(0.7, 0.8), c(0.8, 0.9), c(0.9, 1))
f1samp <- thisout.ff['f1samp',i]*meansds$f1samp.sd + meansds$f1samp.mean
fsi <- thisout.ff['fsiprime',i]*meansds$fsi.sd + meansds$fsi.mean
s <- thisout.ff['s',i]
xlims <- range(fsi, na.rm=TRUE)
for(j in 1:length(ints)){
	k <- which(f1samp > ints[[j]][1] & f1samp < ints[[j]][2])
	plot(fsi[k], s[k], xlab='fsi', ylab='s', main=paste('f1samp', paste(ints[[j]], collapse='-')), xlim=xlims, ylim=c(-1,1), cex=0.5, col=rgb(0.2, 0.2, 0.2, 0.2))
}
dev.off()



## clean up
	delete(thisout.ff)
	rm(thisout.ff)
