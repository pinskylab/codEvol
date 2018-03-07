# create an LD decay plot

require(data.table)
# require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

require(RColorBrewer)

######################
# calculate LD decay
######################

# read in data
dat14 <- fread('zcat analysis/LOF_S_14.geno.ld')
dat11 <- fread('zcat analysis/LOF_S_11.geno.ld')
dat07 <- fread('zcat analysis/LOF_07.geno.ld')
dat40 <- fread('zcat analysis/CAN40.geno.ld')
datMod <- fread('zcat analysis/CANMod.geno.ld')

# remove ^ from a column name
setnames(dat14, 5, 'r2')
setnames(dat11, 5, 'r2')
setnames(dat07, 5, 'r2')
setnames(dat40, 5, 'r2')
setnames(datMod, 5, 'r2')

# calculate distance
dat14[,dist:=abs(POS2-POS1)]
dat11[,dist:=abs(POS2-POS1)]
dat07[,dist:=abs(POS2-POS1)]
dat40[,dist:=abs(POS2-POS1)]
datMod[,dist:=abs(POS2-POS1)]

# plot raw data
#dat[sample(1:nrow(dat), 10000),plot(dist, r2)] # noisy

# calculate averages within bins (like Bario et al. 2016 eLife)
stp1 <- 5 # step size for small distances
thresh <- 50 # use stp2 above this distance
stp2 <- 500

dat14[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat14[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class
dat11[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat11[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class
dat07[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat07[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class
dat40[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
dat40[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class
datMod[,distclass:=floor(dist/stp1)*stp1+stp1/2] # calculate a distance class
datMod[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2] # calculate a distance class

bins14 <- dat14[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass] # average within distance classes
bins11 <- dat11[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]
bins07 <- dat07[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]
bins40 <- dat40[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]
binsMod <- datMod[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]

bins14[,r2se:=r2sd/r2n]
bins11[,r2se:=r2sd/r2n]
bins07[,r2se:=r2sd/r2n]
bins40[,r2se:=r2sd/r2n]
binsMod[,r2se:=r2sd/r2n]

setkey(bins14, distclass)
setkey(bins11, distclass)
setkey(bins07, distclass)
setkey(bins40, distclass)
setkey(binsMod, distclass)

# label by population
bins14[,pop:='LOF_S_14']
bins11[,pop:='LOF_S_11']
bins07[,pop:='LOF_07']
bins40[,pop:='CAN40']
binsMod[,pop:='CANMod']

# merge
bins <- rbind(bins07, bins11, bins14, bins40, binsMod)

# write out
write.csv(bins, file='analysis/ld_decay.csv')

###############
# plot
###############
bins <- fread('analysis/ld_decay.csv')
cols <- brewer.pal(5, 'Set1')
cex=0.5

# plot binned data
quartz(width=3, height=3)
# pdf(width=3, height=3, file='analysis/figures/ld_decay.pdf')
par(mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)
bins[pop=='LOF_07',plot(distclass, r2ave, ylim=c(0,0.7), type='o', xlab='Distance (bp)', ylab='Average correlation (r2)', cex=cex, main='LD decay', log='x', col=cols[1])]
bins[pop=='LOF_S_11',lines(distclass, r2ave, type='o', cex=cex, col=cols[2])]
bins[pop=='LOF_S_14',lines(distclass, r2ave, type='o', cex=cex, col=cols[3])]
bins[pop=='CAN40',lines(distclass, r2ave, type='o', cex=cex, col=cols[4])]
bins[pop=='CANMod',lines(distclass, r2ave, type='o', cex=cex, col=cols[5])]

bins[pop=='LOF_07',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[1]), by=distclass] # 95% CI. too small to be visible
bins[pop=='LOF_S_11',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[2]), by=distclass]
bins[pop=='LOF_S_14',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[3]), by=distclass]
bins[pop=='CAN40',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[4]), by=distclass]
bins[pop=='CANMod',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[5]), by=distclass]

legend('topright', legend=c('LOF_07', 'LOF_S_11', 'LOF_S_14', 'CAN40', 'CANMod'), col=cols, pch=1, lty=1, bty='n', cex=0.6)

dev.off()