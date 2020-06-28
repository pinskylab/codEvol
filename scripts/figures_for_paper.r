#################
## Parameters and functions
#################
require(data.table)
require(RColorBrewer)

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
colout <- '#b2182b' # red, part of RdGy colorbrewer, for outlier loci
#cols3 <- c('#F7941E', '#F2592A', '#BF1E2D', '#FFDE17', '#BEA512') # Lof07, Lof11, Lof14, Can40, CanMod (reds, yellows, from Bastiaan)
cols3 <- c('#e7d4e8', '#af8dc3', '#762a83', '#d9f0d3', '#1b7837') # Lof07, Lof11, Lof14, Can40, CanMod (PRGn colorbrewer)

# takes p-values and returns a color
# uses colout for p<thresh, otherwise a colorbrewer palette
# alternates color by chromosome
lgcolsramp <- function(x, lg = 1, thresh = 0.05){
  if(lg == 1) rmp <- colorRamp(brewer.pal(11, 'BrBG')[5:1])
  if(lg == 2) rmp <- colorRamp(brewer.pal(11, 'BrBG')[7:11])
  if(!(lg %in% c(1,2))) error('lg must be 1 or 2') 
  out <- rep(colout, length(x))
  out[x <= thresh] <- colout
  rgbs <- rmp(-log10(x[x > thresh])/-log10(thresh))
  out[x > thresh] <- rgb(rgbs[,1], rgbs[,2], rgbs[,3], maxColorValue = 256)
  return(out)
}

# draw a vertical colorbar at xposmin to xposmax, from yposmin to yposmax
# cols in hex, x are values to go with each col
color.bar <- function(cols, x, axis = TRUE, cex = 1, nticks=11, 
                      xposmin, xposmax, yposmin, yposmax, title = '', titley = yposmax) {
  scale = (length(cols)-1)/(yposmax - yposmin)
  
  for (i in 1:(length(cols)-1)) {
    y = (i-1)/scale + yposmin
    rect(xposmin,y,xposmax,y+1/scale, col=cols[i], border=NA)
  }
  
  if(axis){
    line(c(xposmax, xposmax), c(yposmin, yposmax))
    line(c(xposmax, xposmax + (xposmax-xposmin)/2), c(yposmin, yposmin))
    line(c(xposmax, xposmax + (xposmax-xposmin)/2), c(yposmax, yposmax))
    at = seq(yposmin, yposmax, length.out = nticks)
    labs <- signif(x[round(seq(1, length(x), length.out = nticks))], 2)
    text(rep(xposmax + (xposmax-xposmin)/2, nticks), y = at, labels = labs, adj = 0, cex = cex)
  }
  
 text(xposmin, titley, title, adj = 0.5, cex = cex)
}

###############################
## Fig. 2 Manhattan plot FSTs
###############################

# read in data: sliding window fst and site-shuffle p-values from ANGSD (GATK nodam2 unlinked sites)
dat <- fread('output/fst_siteshuffle.angsd.gatk.csv.gz') # output by angsd_fst_siteshuffle_null_stats.r

# LG mid-points and genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

setkey(dat, CHROM)
dat <- dat[chrmax[, .(CHROM = chr, start)], ]
dat[, posgen := midPos + start]
dat[,start := NULL]

# trim out Unplaced or windows with <2 loci
dat <- dat[!(CHROM %in% c('Unplaced')) & nloci > 1, ]
	nrow(dat)

# report stats
dat[pop == 'can', min(p)]
dat[pop == 'lof0711', min(p)]
dat[pop == 'lof0714', min(p)]
dat[pop == 'lof1114', min(p)]

# add a vector for color by LG
lgs <- dat[, sort(unique(CHROM))]
dat[,lgcol := lgcolsramp(p, lg = 1, thresh = 0.1)]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := lgcolsramp(p, lg = 2, thresh = 0.1)]

# and for plotting a colorbar
cbar <- data.frame(logx = seq(0, 1.30103, length.out=50))
cbar$x <- 10^(-cbar$logx)
cbar$col1 <- lgcolsramp(cbar$x, lg = 1, thresh = 0.1)
cbar$col2 <- lgcolsramp(cbar$x, lg = 2, thresh = 0.1)

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

# quartz(height=8, width=6)
png(height=5, width=6, units='in', res=300, file='figures/figure2.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.6, 0.1, 0.1))

ymax <- dat[, max(fst, na.rm=TRUE)]
xlims <- dat[, range(posgen, na.rm=TRUE)]

dat[pop == 'can', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)

legend('topleft', legend = c(2, 5, 10, 100), pch = 1, pt.cex = log(c(2,5,10,100))*cexsc, title = '# of SNPs', bty = 'n', cex = 0.5)
color.bar(cbar$col1, cbar$x, axis = FALSE, nticks = 5, xposmin =52e6, xposmax = 56e6, yposmin = 0.2, yposmax = 0.35)
color.bar(cbar$col2, cbar$x, cex = 0.5, axis = TRUE, nticks = 5, 
          xposmin =56e6, xposmax = 60e6, yposmin = 0.2, yposmax = 0.35, title = 'p-value', titley = 0.375)


dat[pop == 'lof0711', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)


dev.off()


####################################################
## Fig. S5 Manhattan plot FST site-shuffle p-value
####################################################

# read in data: sliding window fst and site-shuffle p-values from ANGSD (GATK nodam2 unlinked sites)
dat <- fread('output/fst_siteshuffle.angsd.gatk.csv.gz') # output by angsd_fst_siteshuffle_null_stats.r

# LG mid-points and genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

setkey(dat, CHROM)
dat <- dat[chrmax[, .(CHROM = chr, start)], ]
dat[, posgen := midPos + start]
dat[,start := NULL]

# trim out Unplaced or windows with <2 loci
dat <- dat[!(CHROM %in% c('Unplaced')) & nloci > 1, ]
nrow(dat)

# add a vector for color by LG
lgs <- dat[, sort(unique(CHROM))]
dat[,lgcol := lgcolsramp(0.2, lg = 1, thresh = 0.1)]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := lgcolsramp(0.2, lg = 2, thresh = 0.1)]

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS5.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.6, 0.1, 0.1))

ymax <- max(c(dat[, max(-log10(p), na.rm=TRUE)], -log10(0.05)))
xlims <- dat[, range(posgen, na.rm=TRUE)]

dat[pop == 'can', plot(posgen, -log10(p), type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                       ylim = c(0,ymax), ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = -log10(0.05), lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)

legend('topleft', legend = c(2, 5, 10, 100), pch = 1, pt.cex = log(c(2,5,10,100))*cexsc, title = '# of SNPs', bty = 'n', cex = 0.5)


dat[pop == 'lof0711', plot(posgen, -log10(p), type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                           ylim = c(0,ymax), ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = -log10(0.05), lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(posgen, -log10(p), type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                           ylim = c(0,ymax), ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = -log10(0.05), lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(posgen, -log10(p), type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                           ylim = c(0,ymax), ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = -log10(0.05), lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)


dev.off()


####################################################
## Fig. S6 Manhattan plot change in pi by region
####################################################
winsz = 5e4 # for scaling the windowed pi and thetaW values

# read in data: sliding window pi from ANGSD (GATK nodam2 unlinked sites)
dat <- fread('analysis/theta_siteshuffle.angsd.gatk.csv.gz') # output by angsd_theta_siteshuffle_null_stats.r

# LG mid-points
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

# add a vector for color by LG
lgs <- dat[, sort(unique(Chromo))]
dat[,lgcol := lgcolsramp(0.2, lg = 1, thresh = 0.1)]
dat[Chromo %in% lgs[seq(2, length(lgs),by=2)], lgcol := lgcolsramp(0.2, lg = 2, thresh = 0.1)]

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS6.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1), mgp = c(4, 1, 0))

ylims <- dat[, range(tPd/winsz, na.rm=TRUE)]
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, tPd/winsz, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression("Change in" ~ pi), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)


dat[pop == 'lof0711', plot(POSgen, tPd/winsz, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression("Change in" ~ pi), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(POSgen, tPd/winsz, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression("Change in" ~ pi), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(POSgen, tPd/winsz, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression("Change in" ~ pi), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)


dev.off()



######################
## Table S5 Outliers
######################
library(data.table)

dat <- fread('tables/outlier_annotation.csv')

# add column for functional location
dat[,annotation:=as.character(NA)]
dat[NearGene != '', annotation:='<25kb from gene']
dat[grepl('CDS', feature), annotation:='coding']
dat[grepl('three_prime_UTR', feature), annotation:="3' UTR"]
dat[grepl('five_prime_UTR', feature), annotation:="5' UTR"]
dat[grepl('mRNA', feature) & !grepl('CDS|three_prime_UTR|five_prime_UTR', feature), annotation:="transcript"]




	
# make Table S5
out <- dat[,.(LG=CHROM, Position=POS, Pop=comp, Test=paste0(test, ' ', gsub(' p-value', '', testvaltype)), p = signif(testval,2), annotation,
	gene_name=lapply(strsplit(paste(WithinAnno, NearAnno), split=':'), FUN=function(x) return(x[1])))]


# convert to text and turn NA to ''
for (j in names(out)) set(out, j = j, value = as.character(out[[j]]))
for (j in names(out)) set(out, i=which(is.na(out[[j]])), j = j, value = '')

# write out
write.csv(out, file='figures/tableS5.csv', row.names=FALSE)
