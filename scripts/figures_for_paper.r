#################
## Parameters and functions
#################
source('scripts/error.bar.R')
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


cols <- c('#a6cee333', '#1f78b433') # light blue, blue, partially transparent
cols2 <- c('#a6cee3', '#1f78b4') # light blue, blue
cols3 <- c('#b2df8a', '#a6cee3') # green, blue light
cols4 <- c('#b2df8a55', '#a6cee355') # green, blue light, partially transparent
cols5 <- c('#ff7f00', '#e31a1c') # orange/red
cols6 <- c('#ff7f0055', '#e31a1c55') # orange/red, partially transparent
lw <- 2 # line width


######################
## Fig. 3 Outliers
######################

dat <- fread('gunzip -c analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz')

# trim
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

# set up plot
quartz(height=10, width=13)
# png(height=3, width=6.5, units='in', res=300, file='figures/figure3.png')
layout(matrix(1:3, nrow=3))
par(las=1, mai=c(0.25, 0.8, 0.1, 0.1))

# A: Manhattan plot
ymax <- dat[dpCanFlag==TRUE | dpLofFlag==TRUE, max(-log10(c(pLof071114, pCan, p.comb071114Can)), na.rm=TRUE)]
xlims <- dat[dpCanFlag==TRUE | dpLofFlag==TRUE, range(POSgen/1e6, na.rm=TRUE)]
fdrln <- min(c(dat[outlierLof071114_Can_q3==1, -log10(p.comb071114Can)], dat[outlierCan_q3==1, -log10(pCan)], dat[outlierLof071114_q3==1, -log10(pLof071114)])) # find FDR cutoff across all 3 populations

dat[dpCanFlag==TRUE  & dpLofFlag==TRUE,plot(0, 0, type='n', xaxt='n', xlab='', ylab='-log10(p)', bty='l', xlim=xlims, ylim=c(0,ymax), xaxs='i', cex.axis=1.2, cex.lab=1.5)] # set up axes
opar <- par(mgp=c(2,0.5,0)); axis(side=1, at=chrmax$mid/1e6, labels=chrmax$CHROM, tick=FALSE, las=3, hadj=1, cex.axis=1.2); par(opar) # plot x-axis
abline(h=fdrln, lty=2) # draw line for FDR cutoff
dat[dpCanFlag==TRUE  & dpLofFlag==TRUE & outlierLof071114_Can_q3==0,points(POSgen/1e6, -log10(p.comb071114Can), type='p', cex=0.2, col=lgcol)] # plot non-outliers points
dat[outlierLof071114_Can_q3==1, points(POSgen/1e6, -log10(p.comb071114Can), type='p', cex=1, col=lgcol2)] # plot outliers combined
dat[outlierLof071114_q3==1, points(POSgen/1e6, -log10(pLof071114), type='p', cex=1, pch= 4, col=lgcol2)] # plot outliers Lof
dat[outlierCan_q3==1, points(POSgen/1e6, -log10(pCan), type='p', cex=1, pch = 17, col=lgcol2)] # plot outliers Can

legend('topleft', pch=c(1,4,17), col=cols2[2], legend=c('Combined', 'Norway', 'Canada'), bty='n')


