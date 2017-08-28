# create an LD decay plot

require(data.table)

# read in data
dat <- fread('analysis/LOF_S_14.geno.ld')
setnames(dat, 5, 'r2')

# calculate distance
dat[,dist:=abs(POS2-POS1)]

# plot raw data
dat[sample(1:nrow(dat), 10000),plot(dist, r2)] # noisy

# calculate averages within bins (like Bario et al. 2016 eLife)
stp1 <- 5 # step size for small distances
thresh <- 50 # use stp2 above this distance
stp2 <- 500

dat[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

bins <- dat[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2n=.N), by=distclass] # average within distance classes
bins[,r2se:=r2sd/r2n]
setkey(bins, distclass)

# plot binned data
quartz(width=3, height=3)
# pdf(width=3, height=3, file='analysis/figures/ld_decay_LOF_S_14.pdf')
par(mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.5, 0), tcl=-0.2)
bins[,plot(distclass, r2ave, ylim=c(0,0.7), type='o', xlab='Distance (bp)', ylab='Average correlation (r2)', cex=0.7, main='LOF_S_14', log='x')]
bins[,lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se)), by=distclass] # 95% CI. too small to be visible

dev.off()