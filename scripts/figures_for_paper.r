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
cols2 <- c('#a6cee3', '#1f78b4') # light blue, blue: for alternating LGs outliers

#cols3 <- c('#F7941E', '#F2592A', '#BF1E2D', '#FFDE17', '#BEA512') # Lof07, Lof11, Lof14, Can40, CanMod (reds, yellows, from Bastiaan)
cols3 <- c('#e7d4e8', '#af8dc3', '#762a83', '#d9f0d3', '#1b7837') # Lof07, Lof11, Lof14, Can40, CanMod (PRGn colorbrewer)


######################
## Fig. 3 Outliers
######################

# read in files
dat <- fread('gunzip -c analysis/wfs_nullmodel_outliers_07-11-14_Can.csv.gz')
freqs <- fread('analysis/freq_init_final_outliers.csv')
ld <- fread('analysis/ld_outliers_union.csv')
pi <- fread('analysis/pi_outliers_union.csv', )
d <- fread('analysis/TajimaD_outliers_1000bp_union.csv', )

# trim loci
dat <- dat[!(CHROM %in% c('Unplaced'))]
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
ymax <- dat[, max(-log10(c(pLof0711, pLof0714, pCan, p.comb0711Can, p.comb0714Can)), na.rm=TRUE)]
xlims <- dat[, range(POSgen/1e6, na.rm=TRUE)]
fdrln <- min(c(dat[q3.comb0711Can<0.2, -log10(p.comb0711Can)], dat[q3.comb0714Can<0.2, -log10(p.comb0714Can)], dat[q3.Can<0.2, -log10(pCan)], dat[q3.Lof0711<0.2, -log10(pLof0711)], dat[q3.Lof0714<0.2, -log10(pLof0714)])) # find FDR cutoff across all 3 populations

plot(0, 0, type='n', xaxt='n', xlab='', ylab=expression(-log[10]*(italic(p))), bty='l', xlim=xlims, ylim=c(0,ymax), xaxs='i', cex.axis=1, cex.lab=1.5) # set up axes
opar <- par(mgp=c(2,0.5,0)); axis(side=1, at=chrmax$mid/1e6, labels=gsub('LG|LG0', '', chrmax$CHROM), tick=FALSE, las=0, hadj=0.5, cex.axis=1); par(opar) # plot x-axis

abline(h=fdrln, lty=2) # draw line for FDR cutoff

dat[!(q3.comb0711Can<0.2 | q3.comb0714Can<0.2 | q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2),points(POSgen/1e6, -log10(p.comb0711Can), type='p', cex=0.2, col=lgcol)] # plot non-outliers points
dat[q3.comb0711Can<0.2, points(POSgen/1e6, -log10(p.comb0711Can), type='p', cex=1, col=lgcol2)] # plot outliers combined
dat[q3.comb0714Can<0.2, points(POSgen/1e6, -log10(p.comb0714Can), type='p', cex=1, col=lgcol2)] # plot outliers combined
dat[q3.Lof0711<0.2, points(POSgen/1e6, -log10(pLof0711), type='p', cex=1, pch= 4, col=lgcol2)] # plot outliers Lof
dat[q3.Lof0714<0.2, points(POSgen/1e6, -log10(pLof0714), type='p', cex=1, pch= 4, col=lgcol2)] # plot outliers Lof
dat[q3.Can<0.2, points(POSgen/1e6, -log10(pCan), type='p', cex=1, pch = 17, col=lgcol2)] # plot outliers Can

legend('bottomright', pch=c(1,4,17), col=cols2[2], legend=c('Combined', 'Norway', 'Canada'), ncol=3, bty='o', cex=0.6, bg='white', box.lwd=0)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)


### B: Norway allele frequencies
ylims <- c(0,10.5)
lwd=1

plot(0,0, type='n', ylim=ylims, xlim=c(0,1), xlab='Allele frequency', ylab='Density', main='', bty='l', cex.lab=1.5)

	# initial
freqs[outlierLof071114_q3==1 | outlierLof071114_Can_q3==1,lines(density(Freq_07low, from=0, to=1), col=cols3[1], lwd=lwd)] # union
freqs[notoutlierLof==TRUE,lines(density(Freq_07low, from=0, to=1), col=cols3[1], lty=3)]

	# final
freqs[outlierLof071114_q3==1 | outlierLof071114_Can_q3==1,lines(density(Freq_11low, from=0, to=1), col=cols3[2], lwd=lwd)]
freqs[outlierLof071114_q3==1 | outlierLof071114_Can_q3==1,lines(density(Freq_14low, from=0, to=1), col=cols3[3], lwd=lwd)]

freqs[notoutlierLof==TRUE,lines(density(Freq_11low, from=0, to=1), col=cols3[2], lty=3)]
freqs[notoutlierLof==TRUE,lines(density(Freq_14low, from=0, to=1), col=cols3[3], lty=3)]

legend('topright', legend=c('Norway 1907', 'Norway 2011', 'Norway 2014', 'Norway outlier', 'Not outlier'), col=c(cols3[1:3], 'black', 'black'), lty=c(rep(1,3), 1, 3), bty='n', cex=0.8)
mtext(side=3, 'B', adj=adjlet2, line=linelet2, cex=cexlet)


# C: Canada allele frequencies
ylims <- c(0,11)
lwd<-1
plot(0,0, type='n', ylim=ylims, xlim=c(0,1), xlab='Allele frequency', ylab='Density', main='', bty='l', cex.lab=1.5)

	# initial
freqs[outlierCan_q3==1 | outlierLof071114_Can_q3==1,lines(density(Freq_Can40low, from=0, to=1), col=cols3[4], lwd=lwd)]
freqs[notoutlierCAN==TRUE,lines(density(Freq_Can40low, from=0, to=1), col=cols3[4], lty=3)]

	# final
freqs[outlierCan_q3==1 | outlierLof071114_Can_q3==1,lines(density(Freq_CanModlow, from=0, to=1), col=cols3[5], lwd=lwd)]
freqs[notoutlierCAN==TRUE,lines(density(Freq_CanModlow, from=0, to=1), col=cols3[5], lty=3)]

legend('topright', legend=c('Canada 1940', 'Canada 2013', 'Canada outlier', 'Not outlier'), col=c(cols3[4:5], 'black', 'black'), lty=c(rep(1,2), 1, 3), bty='n', cex=0.8)
mtext(side=3, 'C', adj=adjlet2, line=linelet2, cex=cexlet)


### D: Norway LD
ofs=c(1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4) # small offsets so error bars don't overlap
lwd=1

plot(1, 1, xlim=c(15,5000), ylim=c(0,1), type='n', xlab='Distance (bp)', ylab=expression(r^2), log='x', bty='l', cex.lab=1.5)

ld[pop=='LOF_07' & outlierclass=='outlier',lines(distclass*ofs[1], r2ave, type='l', col=cols3[1], lwd=lwd)] # Norway outliers
ld[pop=='LOF_S_11' & outlierclass=='outlier',lines(distclass*ofs[2], r2ave, type='l', col=cols3[2], lwd=lwd)]
ld[pop=='LOF_S_14' & outlierclass=='outlier',lines(distclass*ofs[3], r2ave, col=cols3[3], lwd=lwd)]

ld[pop=='LOF_07' & outlierclass=='notoutlier',lines(distclass*ofs[4], r2ave, type='l', col=cols3[1], lty=3, lwd=1)] # not outliers
ld[pop=='LOF_S_11' & outlierclass=='notoutlier',lines(distclass*ofs[5], r2ave, type='l', col=cols3[2], lty=3, lwd=1)]
ld[pop=='LOF_S_14' & outlierclass=='notoutlier',lines(distclass*ofs[6], r2ave, type='l', col=cols3[3], lty=3, lwd=1)]

ld[pop=='LOF_07' & outlierclass=='outlier',lines(c(distclass, distclass)*ofs[1], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[1]), by=distclass] # 95% CI on the mean.
ld[pop=='LOF_S_11' & outlierclass=='outlier',lines(c(distclass, distclass)*ofs[2], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[2]), by=distclass]
ld[pop=='LOF_S_14' & outlierclass=='outlier',lines(c(distclass, distclass)*ofs[3], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[3]), by=distclass]

ld[pop=='LOF_07' & outlierclass=='notoutlier',lines(c(distclass, distclass)*ofs[4], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[1]), by=distclass]
ld[pop=='LOF_S_11' & outlierclass=='notoutlier',lines(c(distclass, distclass)*ofs[5], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[2]), by=distclass]
ld[pop=='LOF_S_14' & outlierclass=='notoutlier',lines(c(distclass, distclass)*ofs[6], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[3]), by=distclass]

mtext(side=3, 'D', adj=adjlet2, line=linelet2, cex=cexlet)


### E: Canada LD
lwd=1

plot(1, 1, xlim=c(15,5000), ylim=c(0,1), type='n', xlab='Distance (bp)', ylab=expression(r^2), log='x', bty='l', cex.lab=1.5)

ld[pop=='CAN40' & outlierclass=='outlier',lines(distclass*ofs[1], r2ave, col=cols3[4], lwd=lwd)]
ld[pop=='CANMod' & outlierclass=='outlier',lines(distclass*ofs[2], r2ave, col=cols3[5], lwd=lwd)]

ld[pop=='CAN40' & outlierclass=='notoutlier',lines(distclass*ofs[3], r2ave, type='l', col=cols3[4], lty=3, lwd=lwd)]
ld[pop=='CANMod' & outlierclass=='notoutlier',lines(distclass*ofs[4], r2ave, type='l', col=cols3[5], lty=3, lwd=lwd)] # not outliers

ld[pop=='CAN40' & outlierclass=='outlier',lines(c(distclass, distclass)*ofs[1], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[4]), by=distclass]
ld[pop=='CANMod' & outlierclass=='outlier',lines(c(distclass, distclass)*ofs[2], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[5]), by=distclass]
ld[pop=='CAN40' & outlierclass=='notoutlier',lines(c(distclass, distclass)*ofs[3], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[4]), by=distclass]
ld[pop=='CANMod' & outlierclass=='notoutlier',lines(c(distclass, distclass)*ofs[4], c(r2ave-1.96*r2se, r2ave+1.96*r2se), col=cols3[5]), by=distclass]

mtext(side=3, 'E', adj=adjlet2, line=linelet2, cex=cexlet)


# F: Norway pi
ylims <- pi[!is.na(disttype),range(c(piave-1.96*pise, piave+1.96*pise), na.rm=TRUE)]

plot(1, 1, xlim=c(5,5000), ylim=ylims, type='n', xlab='Distance (bp)', ylab='Pi per site', main='', log='x', bty='l', cex.lab=1.5)

pi[pop=='LOF_07' & disttype=='outlier',lines(distclass, piave, col=cols3[1], lwd=lwd)]
pi[pop=='LOF_S_11' & disttype=='outlier',lines(distclass, piave, col=cols3[2], lwd=lwd)]
pi[pop=='LOF_S_14' & disttype=='outlier',lines(distclass, piave, col=cols3[3], lwd=lwd)]

pi[pop=='LOF_07' & disttype=='notoutlier',lines(distclass, piave, type='l', col=cols3[1], lty=3, lwd=lwd)]
pi[pop=='LOF_S_11' & disttype=='nooutlier',lines(distclass, piave, type='l', col=cols3[2], lty=3, lwd=lwd)]
pi[pop=='LOF_S_14' & disttype=='nooutlier',lines(distclass, piave, type='l', col=cols3[3], lty=3, lwd=lwd)]

pi[pop=='LOF_07' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[1]), by=distclass]
pi[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[2]), by=distclass]
pi[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[3]), by=distclass] # 95% CI on the mean.

pi[pop=='LOF_07' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[1]), by=distclass]
pi[pop=='LOF_S_11' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[2]), by=distclass]
pi[pop=='LOF_S_14' & disttype=='notoutlier',lines(c(distclass, distclass), c(piave-1.96*pise, piave+1.96*pise), col=cols3[3]), by=distclass]

mtext(side=3, 'F', adj=adjlet2, line=linelet2, cex=cexlet)

# G: Canada pi
plot(1, 1, xlim=c(5,5000), ylim=ylims, type='n', xlab='Distance (bp)', ylab='Pi per site', main='', log='x', bty='l', cex.lab=1.5)

pi[pop=='CAN40' & disttype=='outlier',lines(distclass, piave, col=cols3[4], lwd=lwd)]
pi[pop=='CANMod' & disttype=='outlier',lines(distclass, piave, col=cols3[5], lwd=lwd)]

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

d[pop=='LOF_07' & disttype=='notoutlier',lines(distclass+ofs[4], Dave, col=cols3[1], lty=3)]
d[pop=='LOF_S_11' & disttype=='nooutlier',lines(distclass+ofs[5], Dave, col=cols3[2], lty=3)]
d[pop=='LOF_S_14' & disttype=='nooutlier',lines(distclass+ofs[6], Dave, col=cols3[3], lty=3)]

d[pop=='LOF_07' & disttype=='outlier',lines(c(distclass, distclass)+ofs[1], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[1]), by=distclass] # 95% CI on the mean.
d[pop=='LOF_S_11' & disttype=='outlier',lines(c(distclass, distclass)+ofs[2], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[2]), by=distclass]
d[pop=='LOF_S_14' & disttype=='outlier',lines(c(distclass, distclass)+ofs[3], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[3]), by=distclass]

d[pop=='LOF_07' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[4], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[1]), by=distclass]
d[pop=='LOF_S_11' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[5], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[2]), by=distclass]
d[pop=='LOF_S_14' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[6], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[3]), by=distclass]

mtext(side=3, 'H', adj=adjlet2, line=linelet2, cex=cexlet)


# I. Canada Tajima’s D
plot(1, 1, xlim=c(0,5000), ylim=ylims, type='n', xlab='Distance (bp)', ylab="Tajima's D", main='', log='', bty='l', cex.lab=1.5)

d[pop=='CAN40' & disttype=='outlier',lines(distclass+ofs[1], Dave, type='l', col=cols3[4], lwd=lwd)]
d[pop=='CANMod' & disttype=='outlier',lines(distclass+ofs[2], Dave, ylim=ylims, col=cols3[5], lwd=lwd)]

d[pop=='CAN40' & disttype=='notoutlier',lines(distclass+ofs[3], Dave, type='l', col=cols3[4], lty=3, lwd=lwd)]
d[pop=='CANMod' & disttype=='notoutlier',lines(distclass+ofs[4], Dave, type='l', col=cols3[5], lty=3, lwd=lwd)]

d[pop=='CAN40' & disttype=='outlier',lines(c(distclass, distclass)+ofs[1], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[4]), by=distclass]
d[pop=='CANMod' & disttype=='outlier',lines(c(distclass, distclass)+ofs[2], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[5]), by=distclass]
d[pop=='CAN40' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[3], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[4]), by=distclass]
d[pop=='CANMod' & disttype=='notoutlier',lines(c(distclass, distclass)+ofs[4], c(Dave-1.96*Dse, Dave+1.96*Dse), col=cols3[5]), by=distclass]

mtext(side=3, 'I', adj=adjlet2, line=linelet2, cex=cexlet)



dev.off()



######################
## Table S5 Outliers
######################
dat <- fread('analysis/outlier_annotation.csv')
hpdsLof <- fread('analysis/wfs_abc_hpds_Lof.csv')
hpdsCan <- fread('analysis/wfs_abc_hpds_Can.csv')

# force POS to numeric
set(dat, j='POS', value=as.numeric(dat$POS))

# add column for functional location
dat[,annotation:=as.character(NA)]
dat[NearGene != '', annotation:='<25kb from gene']
dat[grepl('CDS', feature), annotation:='coding']
dat[grepl('three_prime_UTR', feature), annotation:="3' UTR"]
dat[grepl('five_prime_UTR', feature), annotation:="5' UTR"]
dat[grepl('mRNA', feature) & !grepl('CDS|three_prime_UTR|five_prime_UTR', feature), annotation:="transcript"]

# polarize s posterior for the allele that increases in frequency
hpdsLof[f1samp < f3samp, ':=' (slowLof.mean=s.mean, slowLof.l95=s.l95, slowLof.u95=s.u95)]
hpdsLof[f1samp >= f3samp, ':=' (slowLof.mean=-s.mean, slowLof.l95=-s.l95, slowLof.u95=-s.u95)]
hpdsCan[f1samp < f2samp, ':=' (slowCan.mean=s.mean, slowCan.l95=s.l95, slowCan.u95=s.u95)]
hpdsCan[f1samp >= f2samp, ':=' (slowCan.mean=-s.mean, slowCan.l95=-s.l95, slowCan.u95=-s.u95)]

# combine s estimates
hpdsLof[,sLof := paste(round(slowLof.mean,2), ' [', round(slowLof.l95,2), '-', round(slowLof.u95,2), ']', sep='')]
hpdsCan[,sCan := paste(round(slowCan.mean,2), ' [', round(slowCan.l95,2), '-', round(slowCan.u95,2), ']', sep='')]
hpdsLof[slowLof.l95>slowLof.u95, sLof := paste(round(slowLof.mean,2), ' [', round(slowLof.u95,2), '-', round(slowLof.l95,2), ']', sep='')]
hpdsCan[slowCan.l95>slowCan.u95, sCan := paste(round(slowCan.mean,2), ' [', round(slowCan.u95,2), '-', round(slowCan.l95,2), ']', sep='')]

# merge
setkey(dat, CHROM, POS)
setkey(hpdsLof, CHROM, POS)
setkey(hpdsCan, CHROM, POS)
dat2 <- merge(dat, hpdsLof, all.x=TRUE)
dat2 <- merge( dat2, hpdsCan, all.x=TRUE)
	dim(dat)
	dim(dat2)

# range of s estimates for outlier loci
dat2[q3.Lof071114 < 0.2, summary(slowLof.mean)]
dat2[q3.comb071114Can < 0.2, summary(slowLof.mean)]
dat2[q3.Can < 0.2, summary(slowCan.mean)]
dat2[q3.comb071114Can < 0.2, summary(slowCan.mean)]

mean(c(dat2[q3.Lof071114 < 0.2, slowLof.mean], dat2[q3.Can < 0.2, slowCan.mean])) # mean across Lof and Can outliers
sd(c(dat2[q3.Lof071114 < 0.2, slowLof.mean], dat2[q3.Can < 0.2, slowCan.mean])) # SE across Lof and Can outliers


# make Table S5
out <- dat2[,.(linkagegroup=CHROM, position=POS, norwayFDR=round(q3.Lof071114,2), canadaFDR=round(q3.Can,2), combinedFDR=round(q3.comb071114Can,2), sLof, sCan, annotation,
	gene_name=lapply(strsplit(paste(WithinAnno, NearAnno), split=':'), FUN=function(x) return(x[1])))]


# convert to text and turn NA to ''
for (j in names(out)) set(out, j = j, value = as.character(out[[j]]))
for (j in names(out)) set(out, i=which(is.na(out[[j]])), j = j, value = '')

# write out
write.csv(out, file='figures/tableS5.csv', row.names=FALSE)