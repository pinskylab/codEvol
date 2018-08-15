#################
## Parameters and functions
#################
require(data.table)
require(RColorBrewer)
require(png)

# takes values 0-1 and turns them into colors. takes ColorBrewer pallete as argument.
colramp <- function(x, pal='RdBu', alpha=255){
	cr <- colorRamp(brewer.pal(9, pal))
	i <- !is.na(x)
	cols <- rgb(cr(x[i]), alpha=rep(alpha, sum(i)), maxColorValue=255)
	out <- rep(NA, length(x))
	out[i] <- cols
	return(out)
}


cols <- c('#a6cee333', '#1f78b433') # light blue, blue, partially transparent: for alternating LGs
cols2 <- c('#a6cee3', '#1f78b4') # light blue, blue: for alternating LGs

cols3 <- c('#F7941E', '#F2592A', '#BF1E2D', '#FFDE17', '#BEA512') # Lof07, Lof11, Lof14, Can40, CanMod

cols4 <- c('#b2df8a55', '#a6cee355') # green, blue light, partially transparent
cols5 <- c('#ff7f00', '#e31a1c') # orange/red
cols6 <- c('#ff7f0055', '#e31a1c55') # orange/red, partially transparent
lw <- 2 # line width


######################
## Fig. 3 Outliers
######################

# read in files
dat <- fread('gunzip -c analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz')
freqs <- fread('analysis/freq_init_final_outliers.csv')
ld <- fread('analysis/ld_outliers.csv')
ldcomb <- fread('analysis/ld_outliers_combinedpop.csv')
pi <- fread('analysis/pi_outliers.csv', )
picomb <- fread('analysis/pi_outliers_combinedpop.csv', )
d <- fread('analysis/TajimaD_outliers_1000bp.csv', )
dcomb <- fread('analysis/TajimaD_outliers_1000bp_combinedpop.csv', )

# trim loci
dat <- dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & kmer25==1]
	nrow(dat)

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- dat[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
chrmax[, mid := start + len/2]

setkey(chrmax, CHROM)
setkey(dat, CHROM)

dat <- dat[chrmax[,.(CHROM, start)], ]
dat[,POSgen:=POS+start]

# add a vector for color by LG
dat[,lgcol := cols[1]]
dat[CHROM %in% chrmax$CHROM[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]
dat[,lgcol2 := cols2[1]]
dat[CHROM %in% chrmax$CHROM[seq(2, nrow(chrmax),by=2)], lgcol2 := cols2[2]]





### set up plot
quartz(height=8, width=6)
# png(height=8, width=6, units='in', res=300, file='figures/figure3.png')
layout(matrix(c(1, 1:9), nrow=5, byrow=TRUE))
par(las=1, mai=c(0.6, 0.6, 0.1, 0.1))
adjlet <- -0.07 # adjustment for subplot letter A
adjlet2 <- -0.18 # adjustment for subplot letter B-I
cexlet <- 1
linelet <- -0.7
linelet2 <- 0

### A: Manhattan plot
ymax <- dat[dpCanFlag==TRUE | dpLofFlag==TRUE, max(-log10(c(pLof071114, pCan, p.comb071114Can)), na.rm=TRUE)]
xlims <- dat[dpCanFlag==TRUE | dpLofFlag==TRUE, range(POSgen/1e6, na.rm=TRUE)]
fdrln <- min(c(dat[outlierLof071114_Can_q3==1, -log10(p.comb071114Can)], dat[outlierCan_q3==1, -log10(pCan)], dat[outlierLof071114_q3==1, -log10(pLof071114)])) # find FDR cutoff across all 3 populations

plot(0, 0, type='n', xaxt='n', xlab='', ylab=expression(-log[10]*(italic(p))), bty='l', xlim=xlims, ylim=c(0,ymax), xaxs='i', cex.axis=1, cex.lab=1.5) # set up axes
opar <- par(mgp=c(2,0.5,0)); axis(side=1, at=chrmax$mid/1e6, labels=gsub('LG|LG0', '', chrmax$CHROM), tick=FALSE, las=0, hadj=0.5, cex.axis=1); par(opar) # plot x-axis

abline(h=fdrln, lty=2) # draw line for FDR cutoff

dat[dpCanFlag==TRUE  & dpLofFlag==TRUE & outlierLof071114_Can_q3==0,points(POSgen/1e6, -log10(p.comb071114Can), type='p', cex=0.2, col=lgcol)] # plot non-outliers points
dat[outlierLof071114_Can_q3==1, points(POSgen/1e6, -log10(p.comb071114Can), type='p', cex=1, col=lgcol2)] # plot outliers combined
dat[outlierLof071114_q3==1, points(POSgen/1e6, -log10(pLof071114), type='p', cex=1, pch= 4, col=lgcol2)] # plot outliers Lof
dat[outlierCan_q3==1, points(POSgen/1e6, -log10(pCan), type='p', cex=1, pch = 17, col=lgcol2)] # plot outliers Can

legend('bottomleft', pch=c(1,4,17), col=cols2[2], legend=c('Combined', 'Norway', 'Canada'), ncol=3, bty='o', cex=0.6, bg='white', box.lwd=0)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)


### B: Norway allele frequencies
ylims <- c(0,10.5)
lwd=1

plot(0,0, type='n', ylim=ylims, xlim=c(0,1), xlab='Allele frequency', ylab='Density', main='', bty='l', cex.lab=1.5)

	# initial
freqs[outlierLof071114_q3==1,lines(density(Freq_07low, from=0, to=1), col=cols3[1], lwd=lwd)]
freqs[outlierLof071114_Can_q3==1,lines(density(Freq_07low, from=0, to=1), col=cols3[1], lwd=lwd, lty=2)]
freqs[notoutlierLof==TRUE,lines(density(Freq_07low, from=0, to=1), col=cols3[1], lty=3)]

	# final
freqs[outlierLof071114_q3==1,lines(density(Freq_11low, from=0, to=1), col=cols3[2], lwd=lwd)]
freqs[outlierLof071114_q3==1,lines(density(Freq_14low, from=0, to=1), col=cols3[3], lwd=lwd)]

freqs[outlierLof071114_Can_q3==1,lines(density(Freq_14low, from=0, to=1), col=cols3[2], lwd=lwd, lty=2)]
freqs[outlierLof071114_Can_q3==1,lines(density(Freq_11low, from=0, to=1), col=cols3[3], lwd=lwd, lty=2)]

freqs[notoutlierLof==TRUE,lines(density(Freq_11low, from=0, to=1), col=cols3[2], lty=3)]
freqs[notoutlierLof==TRUE,lines(density(Freq_14low, from=0, to=1), col=cols3[3], lty=3)]

legend('topright', legend=c('Norway 1907', 'Norway 2011', 'Norway 2014', 'Norway outlier', 'Combined outlier', 'Not outlier'), col=c(cols3[1:3], 'black', 'black', 'black'), lty=c(rep(1,3), 1, 2, 3), bty='n', cex=0.8)
mtext(side=3, 'B', adj=adjlet2, line=linelet2, cex=cexlet)


# C: Canada allele frequencies
ylims <- c(0,11)
lwd<-1
plot(0,0, type='n', ylim=ylims, xlim=c(0,1), xlab='Allele frequency', ylab='Density', main='', bty='l', cex.lab=1.5)

	# initial
freqs[outlierCan_q3==1,lines(density(Freq_Can40low, from=0, to=1, adjust=0.1), col=cols3[4], lwd=lwd)]
freqs[outlierLof071114_Can_q3==1,lines(density(Freq_Can40low, from=0, to=1, adjust=1), col=cols3[4], lwd=lwd, lty=2)]
freqs[notoutlierCAN==TRUE,lines(density(Freq_Can40low, from=0, to=1), col=cols3[4], lty=3)]

	# final
freqs[outlierCan_q3==1,lines(density(Freq_CanModlow, from=0, to=1, adjust=50), col=cols3[5], lwd=lwd)]
freqs[outlierLof071114_Can_q3==1,lines(density(Freq_CanModlow, from=0, to=1, adjust=1), col=cols3[5], lwd=lwd, lty=2)]
freqs[notoutlierCAN==TRUE,lines(density(Freq_CanModlow, from=0, to=1), col=cols3[5], lty=3)]

legend('topright', legend=c('Canada 1940', 'Canada 2013', 'Canada outlier', 'Combined outlier', 'Not outlier'), col=c(cols3[4:5], 'black', 'black', 'black'), lty=c(rep(1,2), 1, 2, 3), bty='n', cex=0.8)
mtext(side=3, 'C', adj=adjlet2, line=linelet2, cex=cexlet)


### D: Norway LD
ofs=c(1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4) # small offsets so error bars don't overlap
lwd=1

plot(1, 1, xlim=c(5,5000), ylim=c(0,1), type='n', xlab='Distance (bp)', ylab=expression(r^2), log='x', bty='l', cex.lab=1.5)

ld[pop=='LOF_07' & outlier==TRUE,lines(distclass*ofs[1], r2ave, type='l', col=cols3[1], lwd=lwd)] # Norway outliers
ld[pop=='LOF_S_11' & outlier==TRUE,lines(distclass*ofs[2], r2ave, type='l', col=cols3[2], lwd=lwd)]
ld[pop=='LOF_S_14' & outlier==TRUE,lines(distclass*ofs[3], r2ave, col=cols3[3], lwd=lwd)]

ldcomb[pop=='LOF_07' & outlier==TRUE,lines(distclass*ofs[4], r2ave, type='l', col=cols3[1], lwd=lwd, lty=2)] # Combined outliers
ldcomb[pop=='LOF_S_11' & outlier==TRUE,lines(distclass*ofs[5], r2ave, type='l', col=cols3[2], lwd=lwd, lty=2)]
ldcomb[pop=='LOF_S_14' & outlier==TRUE,lines(distclass*ofs[6], r2ave, col=cols3[3], lwd=lwd, lty=2)]

ld[pop=='LOF_07' & outlier==FALSE,lines(distclass*ofs[7], r2ave, type='l', col=cols3[1], lty=3, lwd=1)] # not outliers
ld[pop=='LOF_S_11' & outlier==FALSE,lines(distclass*ofs[8], r2ave, type='l', col=cols3[2], lty=3, lwd=1)]
ld[pop=='LOF_S_14' & outlier==FALSE,lines(distclass*ofs[9], r2ave, type='l', col=cols3[3], lty=3, lwd=1)]

ld[pop=='LOF_07' & outlier==TRUE,lines(c(distclass, distclass)*ofs[1], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[1]), by=distclass] # 95% CI on the mean.
ld[pop=='LOF_S_11' & outlier==TRUE,lines(c(distclass, distclass)*ofs[2], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[2]), by=distclass]
ld[pop=='LOF_S_14' & outlier==TRUE,lines(c(distclass, distclass)*ofs[3], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[3]), by=distclass]

ldcomb[pop=='LOF_07' & outlier==TRUE,lines(c(distclass, distclass)*ofs[4], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[1]), by=distclass]
ldcomb[pop=='LOF_S_11' & outlier==TRUE,lines(c(distclass, distclass)*ofs[5], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[2]), by=distclass]
ldcomb[pop=='LOF_S_14' & outlier==TRUE,lines(c(distclass, distclass)*ofs[6], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[3]), by=distclass]

ld[pop=='LOF_07' & outlier==FALSE,lines(c(distclass, distclass)*ofs[7], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[1]), by=distclass]
ld[pop=='LOF_S_11' & outlier==FALSE,lines(c(distclass, distclass)*ofs[8], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[2]), by=distclass]
ld[pop=='LOF_S_14' & outlier==FALSE,lines(c(distclass, distclass)*ofs[9], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[3]), by=distclass]

mtext(side=3, 'D', adj=adjlet2, line=linelet2, cex=cexlet)


### E: Canada LD
lwd=1

plot(1, 1, xlim=c(5,5000), ylim=c(0,1), type='n', xlab='Distance (bp)', ylab=expression(r^2), log='x', bty='l', cex.lab=1.5)

ld[pop=='CAN40' & outlier==TRUE,lines(distclass*ofs[1], r2ave, col=cols3[4], lwd=lwd)]
ld[pop=='CANMod' & outlier==TRUE,lines(distclass*ofs[2], r2ave, col=cols3[5], lwd=lwd)]

ldcomb[pop=='CAN40' & outlier==TRUE,lines(distclass*ofs[3], r2ave, col=cols3[4], lwd=lwd, lty=2)]
ldcomb[pop=='CANMod' & outlier==TRUE,lines(distclass*ofs[4], r2ave, col=cols3[5], lwd=lwd, lty=2)]

ld[pop=='CAN40' & outlier==FALSE,lines(distclass*ofs[5], r2ave, type='l', col=cols3[4], lty=3, lwd=lwd)]
ld[pop=='CANMod' & outlier==FALSE,lines(distclass*ofs[6], r2ave, type='l', col=cols3[5], lty=3, lwd=lwd)] # not outliers

ld[pop=='CAN40' & outlier==TRUE,lines(c(distclass, distclass)*ofs[1], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[4]), by=distclass]
ld[pop=='CANMod' & outlier==TRUE,lines(c(distclass, distclass)*ofs[2], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[5]), by=distclass]
ldcomb[pop=='CAN40' & outlier==TRUE,lines(c(distclass, distclass)*ofs[3], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[4]), by=distclass]
ldcomb[pop=='CANMod' & outlier==TRUE,lines(c(distclass, distclass)*ofs[4], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[5]), by=distclass]
ld[pop=='CAN40' & outlier==FALSE,lines(c(distclass, distclass)*ofs[5], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[4]), by=distclass]
ld[pop=='CANMod' & outlier==FALSE,lines(c(distclass, distclass)*ofs[6], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[5]), by=distclass]

mtext(side=3, 'E', adj=adjlet2, line=linelet2, cex=cexlet)

# F: Norway pi
ylims <- pi[!is.na(disttype),range(c(piave-1.96*pise, piave+1.96*pise), na.rm=TRUE)]

plot(1, 1, xlim=c(5,5000), ylim=ylims, type='n', xlab='Distance (bp)', ylab='Pi per site', main='', log='x', bty='l', cex.lab=1.5)

pi[pop=='LOF_07' & disttype=='outlier',lines(distclass, piave, col=cols3[1], lwd=lwd)]
pi[pop=='LOF_S_11' & disttype=='outlier',lines(distclass, piave, col=cols3[2], lwd=lwd)]
pi[pop=='LOF_S_14' & disttype=='outlier',lines(distclass, piave, col=cols3[3], lwd=lwd)]

picomb[pop=='LOF_07' & disttype=='outlier',lines(distclass, piave, col=cols3[1], lwd=lwd, lty=2)]
picomb[pop=='LOF_S_11' & disttype=='outlier',lines(distclass, piave, col=cols3[2], lwd=lwd, lty=2)]
picomb[pop=='LOF_S_14' & disttype=='outlier',lines(distclass, piave, col=cols3[3], lwd=lwd, lty=2)]

pi[pop=='LOF_07' & disttype=='notoutlier',lines(distclass, piave, type='l', col=cols3[1], lty=3, lwd=lwd)]
pi[pop=='LOF_S_11' & disttype=='nooutlier',lines(distclass, piave, type='l', col=cols3[2], lty=3, lwd=lwd)]
pi[pop=='LOF_S_14' & disttype=='nooutlier',lines(distclass, piave, type='l', col=cols3[3], lty=3, lwd=lwd)]

pi[pop=='LOF_07' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[1]), by=distclass]
pi[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[2]), by=distclass]
pi[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[3]), by=distclass] # 95% CI on the mean.

picomb[pop=='LOF_07' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[1]), by=distclass]
picomb[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[2]), by=distclass]
picomb[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[3]), by=distclass]

pi[pop=='LOF_07' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[1]), by=distclass]
pi[pop=='LOF_S_11' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[2]), by=distclass]
pi[pop=='LOF_S_14' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[3]), by=distclass]

mtext(side=3, 'F', adj=adjlet2, line=linelet2, cex=cexlet)

# G: Canada pi
plot(1, 1, xlim=c(5,5000), ylim=ylims, type='n', xlab='Distance (bp)', ylab='Pi per site', main='', log='x', bty='l', cex.lab=1.5)

pi[pop=='CAN40' & disttype=='outlier',lines(distclass, piave, col=cols3[4], lwd=lwd)]
pi[pop=='CANMod' & disttype=='outlier',lines(distclass, piave, col=cols3[5], lwd=lwd)]

picomb[pop=='CAN40' & disttype=='outlier',lines(distclass, piave, col=cols3[4], lwd=lwd, lty=2)]
picomb[pop=='CANMod' & disttype=='outlier',lines(distclass, piave, col=cols3[5], lwd=lwd, lty=2)]

pi[pop=='CAN40' & disttype=='notoutlier',lines(distclass, piave, col=cols3[4], lty=3, lwd=lwd)]
pi[pop=='CANMod' & disttype=='notoutlier',lines(distclass, piave, col=cols3[5], lty=3, lwd=lwd)]

pi[pop=='CAN40' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[4]), by=distclass]
pi[pop=='CANMod' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[5]), by=distclass]
pi[pop=='CAN40' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[4]), by=distclass]
pi[pop=='CANMod' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[5]), by=distclass]

mtext(side=3, 'G', adj=adjlet2, line=linelet2, cex=cexlet)


# H. Norway Tajima’s D
ofs <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)
ylims <- d[!is.na(disttype),range(c(Dave-1.96*Dse, Dave+1.96*Dse), na.rm=TRUE)]

plot(1, 1, xlim=c(0,5000), ylim=ylims, type='n', xlab='Distance (bp)', ylab="Tajima's D", main='', log='', bty='l', cex.lab=1.5)

d[pop=='LOF_07' & disttype=='outlier',lines(distclass+ofs[1], Dave, col=cols3[1], lwd=lwd)]
d[pop=='LOF_S_11' & disttype=='outlier',lines(distclass+ofs[2], Dave, col=cols3[2], lwd=lwd)]
d[pop=='LOF_S_14' & disttype=='outlier',lines(distclass+ofs[3], Dave, col=cols3[3], lwd=lwd)]

dcomb[pop=='LOF_07' & disttype=='outlier',lines(distclass+ofs[4], Dave, col=cols3[1], lwd=lwd, lty=2)]
dcomb[pop=='LOF_S_11' & disttype=='outlier',lines(distclass+ofs[5], Dave, col=cols3[2], lwd=lwd, lty=2)]
dcomb[pop=='LOF_S_14' & disttype=='outlier',lines(distclass+ofs[6], Dave, col=cols3[3], lwd=lwd, lty=2)]

d[pop=='LOF_07' & disttype=='notoutlier',lines(distclass+ofs[7], Dave, col=cols3[1], lty=3)]
d[pop=='LOF_S_11' & disttype=='nooutlier',lines(distclass+ofs[8], Dave, col=cols3[2], lty=3)]
d[pop=='LOF_S_14' & disttype=='nooutlier',lines(distclass+ofs[9], Dave, col=cols3[3], lty=3)]

d[pop=='LOF_07' & disttype=='outlier',lines(c(distclass, distclass)+ofs[1], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[1]), by=distclass] # 95% CI on the mean.
d[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass, distclass)+ofs[2], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[2]), by=distclass]
d[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass)+ofs[3], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[3]), by=distclass]

dcomb[pop=='LOF_07' & disttype=='outlier',lines(c(distclass, distclass)+ofs[4], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[1]), by=distclass]
dcomb[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass, distclass)+ofs[5], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[2]), by=distclass]
dcomb[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass)+ofs[6], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[3]), by=distclass]

d[pop=='LOF_07' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[7], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[1]), by=distclass]
d[pop=='LOF_S_11' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[8], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[2]), by=distclass]
d[pop=='LOF_S_14' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[9], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[3]), by=distclass]

mtext(side=3, 'H', adj=adjlet2, line=linelet2, cex=cexlet)


# I. Canada Tajima’s D
plot(1, 1, xlim=c(0,5000), ylim=ylims, type='n', xlab='Distance (bp)', ylab="Tajima's D", main='', log='', bty='l', cex.lab=1.5)

d[pop=='CAN40' & disttype=='outlier',lines(distclass+ofs[1], Dave, type='l', col=cols3[4], lwd=lwd)]
d[pop=='CANMod' & disttype=='outlier',lines(distclass+ofs[2], Dave, ylim=ylims, col=cols3[5], lwd=lwd)]

dcomb[pop=='CAN40' & disttype=='outlier',lines(distclass+ofs[3], Dave, type='l', col=cols3[4], lwd=lwd, lty=2)]
dcomb[pop=='CANMod' & disttype=='outlier',lines(distclass+ofs[4], Dave, ylim=ylims, col=cols3[5], lwd=lwd, lty=2)]

d[pop=='CAN40' & disttype=='notoutlier',lines(distclass+ofs[5], Dave, type='l', col=cols3[4], lty=3, lwd=lwd)]
d[pop=='CANMod' & disttype=='notoutlier',lines(distclass+ofs[6], Dave, type='l', col=cols3[5], lty=3, lwd=lwd)]

d[pop=='CAN40' & disttype=='outlier',lines(c(distclass, distclass)+ofs[1], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[4]), by=distclass]
d[pop=='CANMod' & disttype=='outlier',lines(c(distclass, distclass)+ofs[2], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[5]), by=distclass]
dcomb[pop=='CAN40' & disttype=='outlier',lines(c(distclass, distclass)+ofs[3], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[4]), by=distclass]
dcomb[pop=='CANMod' & disttype=='outlier',lines(c(distclass, distclass)+ofs[4], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[5]), by=distclass]
d[pop=='CAN40' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[5], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[4]), by=distclass]
d[pop=='CANMod' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[6], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[5]), by=distclass]

mtext(side=3, 'I', adj=adjlet2, line=linelet2, cex=cexlet)



dev.off()