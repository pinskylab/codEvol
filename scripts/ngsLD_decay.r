# create an LD decay plot from ngsLD output

require(data.table)
require(RColorBrewer)

######################
# calculate LD decay
######################

# read in data (gatk loci only)
datCan40 <- fread('analysis/ld.Can_40.gatk.gz')
datCan14 <- fread('analysis/ld.Can_14.gatk.gz')
dat07 <- fread('analysis/ld.Lof_07.gatk.gz')
dat11 <- fread('analysis/ld.Lof_11.gatk.gz')
dat14 <- fread('analysis/ld.Lof_14.gatk.gz')

# add column names
nms <- c('pos1', 'pos2', 'dist', 'r2', 'D', 'Dprime', 'r2em')
setnames(datCan14, nms)
setnames(datCan40, nms)
setnames(dat14, nms)
setnames(dat11, nms)
setnames(dat07, nms)


# calculate averages within bins (like Bario et al. 2016 eLife)
stp1 <- 5 # step size for small distances
thresh <- 50 # use stp2 above this distance
stp2 <- 500

dat14[, distclass := floor(dist/stp1)*stp1 + stp1/2] # calculate a distance class
dat14[dist > thresh, distclass := floor(dist/stp2)*stp2 + stp2/2] # calculate a distance class for distances above thresh
dat11[,distclass:=floor(dist/stp1)*stp1+stp1/2]
dat11[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2]
dat07[,distclass:=floor(dist/stp1)*stp1+stp1/2]
dat07[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2]
datCan40[,distclass:=floor(dist/stp1)*stp1+stp1/2]
datCan40[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2]
datCan14[,distclass:=floor(dist/stp1)*stp1+stp1/2]
datCan14[dist>thresh,distclass:=floor(dist/stp2)*stp2+stp2/2]

bins14 <- dat14[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass] # average within distance classes
bins11 <- dat11[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]
bins07 <- dat07[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]
binsCan40 <- datCan40[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]
binsCan14 <- datCan14[!is.na(r2),.(r2ave=mean(r2), r2sd=sd(r2), r2l95=quantile(r2,probs=0.025), r2u95=quantile(r2,probs=0.975), r2l75=quantile(r2,probs=0.125), r2u75=quantile(r2,probs=0.875), r2n=.N), by=distclass]

bins14[,r2se:=r2sd/r2n]
bins11[,r2se:=r2sd/r2n]
bins07[,r2se:=r2sd/r2n]
binsCan40[,r2se:=r2sd/r2n]
binsCan14[,r2se:=r2sd/r2n]

setkey(bins14, distclass)
setkey(bins11, distclass)
setkey(bins07, distclass)
setkey(binsCan40, distclass)
setkey(binsCan14, distclass)

# label by population
bins14[,pop:='LOF_S_14']
bins11[,pop:='LOF_S_11']
bins07[,pop:='LOF_07']
binsCan40[,pop:='CAN40']
binsCan14[,pop:='CAN14']

# merge
bins <- rbind(bins07, bins11, bins14, binsCan40, binsCan14)

# write out
write.csv(bins, file='analysis/ld_decay_ngsLD.gatk.csv')

###############
# plot
###############
bins <- fread('analysis/ld_decay_ngsLD.gatk.csv')
cols <- brewer.pal(5, 'Set1')
cex=0.5

# plot binned data
quartz(width=3, height=3)
# pdf(width=3, height=3, file='figures/ld_decay_ngsLD.pdf')
par(mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)
bins[pop=='LOF_07',plot(distclass, r2ave, ylim=c(0,0.7), type='o', xlab='Distance (bp)', ylab='Average correlation (r2)', cex=cex, main='LD decay', log='x', col=cols[1])]
bins[pop=='LOF_S_11',lines(distclass, r2ave, type='o', cex=cex, col=cols[2])]
bins[pop=='LOF_S_14',lines(distclass, r2ave, type='o', cex=cex, col=cols[3])]
bins[pop=='CAN40',lines(distclass, r2ave, type='o', cex=cex, col=cols[4])]
bins[pop=='CAN14',lines(distclass, r2ave, type='o', cex=cex, col=cols[5])]

bins[pop=='LOF_07',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[1]), by=distclass] # 95% CI. too small to be visible
bins[pop=='LOF_S_11',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[2]), by=distclass]
bins[pop=='LOF_S_14',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[3]), by=distclass]
bins[pop=='CAN40',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[4]), by=distclass]
bins[pop=='CAN14',lines(c(distclass, distclass), c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols[5]), by=distclass]

legend('topright', legend=c('LOF_07', 'LOF_S_11', 'LOF_S_14', 'CAN40', 'CAN14'), col=cols, pch=1, lty=1, bty='n', cex=0.6)

dev.off()