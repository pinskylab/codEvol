#################
## Parameters and functions
#################
require(data.table)
require(RColorBrewer)
require(ggplot2) # for Fig. 1 hexplots. also requires package hexbin
require(gridExtra) # for Fig. 1

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

#####################################################
# Part of Fig. 1: allele frequency change heatmap
#####################################################
colLof <- '#F26523' # matching Bastiaan
colCan <- '#BEA512'

# read in and merge data
datCan40 <- fread('data_31_01_20/Can_40_freq.mafs.gz')
datCan14 <- fread('data_31_01_20/Can_14_freq.mafs.gz')
datLof07 <- fread('data_31_01_20/Lof_07_freq.mafs.gz')
datLof11 <- fread('data_31_01_20/Lof_11_freq.mafs.gz')
datLof14 <- fread('data_31_01_20/Lof_14_freq.mafs.gz')

datCan <- merge(datCan40[, .(chromo, position, freq1 = knownEM)], datCan14[, .(chromo, position, freq2 = knownEM)])
datLof0711 <- merge(datLof07[, .(chromo, position, freq1 = knownEM)], datLof11[, .(chromo, position, freq2 = knownEM)])
datLof0714 <- merge(datLof07[, .(chromo, position, freq1 = knownEM)], datLof14[, .(chromo, position, freq2 = knownEM)])
datLof1114 <- merge(datLof11[, .(chromo, position, freq1 = knownEM)], datLof14[, .(chromo, position, freq2 = knownEM)])

# trim to unlinked loci
ulnkCan <- fread('analysis/ld.unlinked.Can.gatk.nodam.csv.gz')
ulnkLof <- fread('analysis/ld.unlinked.Lof.gatk.nodam.csv.gz')

datCan <- merge(datCan, ulnkCan[, .(chromo = CHROM, position = POS)])
datLof0711 <- merge(datLof0711, ulnkLof[, .(chromo = CHROM, position = POS)])
datLof0714 <- merge(datLof0714, ulnkLof[, .(chromo = CHROM, position = POS)])
datLof1114 <- merge(datLof1114, ulnkLof[, .(chromo = CHROM, position = POS)])

# r2 values
datCan[, cor(freq1, freq2)]^2
datLof0711[, cor(freq1, freq2)]^2
datLof0714[, cor(freq1, freq2)]^2
datLof1114[, cor(freq1, freq2)]^2

# plots
bks <- c(0, 1, 10, 100, 1000, 10000)
p1 <- ggplot(datCan, aes(freq1, freq2)) +
  geom_hex() +
  scale_fill_gradient(low="grey", high=colCan, trans = 'log', breaks = bks, labels = bks) +
  labs(title = 'Canada', x = 'Historical', y = 'Modern') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        legend.position = 'none')
p2 <- ggplot(datLof0711, aes(freq1, freq2)) +
  geom_hex() +
  scale_fill_gradient(low="grey", high=colLof, trans = 'log', breaks = bks, labels = bks) +
  labs(title = 'Norway 07-11', x = 'Historical', y = 'Modern') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        legend.position = 'none')
p3 <- ggplot(datLof0714, aes(freq1, freq2)) +
  geom_hex() +
  scale_fill_gradient(low="grey", high=colLof, trans = 'log', breaks = bks, labels = bks) +
  labs(title = 'Norway 07-14', x = 'Historical', y = 'Modern') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        legend.position = 'none')
p4 <- ggplot(datLof1114, aes(freq1, freq2)) +
  geom_hex() +
  scale_fill_gradient(low="grey", high=colLof, trans = 'log', breaks = bks, labels = bks) +
  labs(title = 'Norway 11-14', x = 'Historical', y = 'Modern') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        legend.position = 'none')

p5 <- grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave(p5, filename = 'figures/figure1.png')

########################################
## Fig. 2 Manhattan plot FSTs by region
########################################

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


####################################################
## Fig. S7 Manhattan plot change in D by region
####################################################

# read in data: sliding window D from ANGSD (GATK nodam2 unlinked sites)
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

png(height=5, width=6, units='in', res=300, file='figures/figureS7.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1), mgp = c(4, 1, 0))

ylims <- dat[, range(tDd, na.rm=TRUE)]
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, tDd, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = "Change in D", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)


dat[pop == 'lof0711', plot(POSgen, tDd, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = "Change in D", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(POSgen, tDd, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = "Change in D", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(POSgen, tDd, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = "Change in D", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)


dev.off()



####################################################
## Fig. S8 Manhattan plot change in LD by region
####################################################

# read in data: sliding window D from ngsLD (GATK nodam2 sites)
dat <- fread('analysis/ld_change_region_5e4_ngsLD.gatk.csv.gz'); width='5e4' # output by ngsLD_region_change.r

# Remove unplaced
dat <- dat[chr != 'Unplaced', ]

# LG mid-points and genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]
dat <- merge(dat, chrmax[, .(chr, start)], by = 'chr')
dat[, POSgen := BIN_START + start + as.numeric(width)/2]

# long format for plotting
datl <- melt(dat[, .(chr, POSgen, ld_diff_Can, ld_diff_Lof0711, ld_diff_Lof0714, ld_diff_Lof1114)], id.vars = c('chr', 'POSgen'), variable.name = 'pop', value.name = 'ldchange')
datl[, pop := gsub('ld_diff_', '', pop)]

# add a vector for color by LG
lgs <- datl[, sort(unique(chr))]
datl[,lgcol := lgcolsramp(0.2, lg = 1, thresh = 0.1)]
datl[chr %in% lgs[seq(2, length(lgs),by=2)], lgcol := lgcolsramp(0.2, lg = 2, thresh = 0.1)]

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS8.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1), mgp = c(4, 1, 0))

ylims <- datl[, range(ldchange, na.rm=TRUE)]
xlims <- datl[, range(POSgen, na.rm=TRUE)]

datl[pop == 'Can', plot(POSgen, ldchange, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = "Change in LD", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.75)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)


datl[pop == 'Lof0711', plot(POSgen, ldchange, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = "Change in LD", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.75)
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

datl[pop == 'Lof0714', plot(POSgen, ldchange, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = "Change in LD", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

datl[pop == 'Lof1114', plot(POSgen, ldchange, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = "Change in LD", bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)


dev.off()



####################################################
## Fig. S9 Manhattan plot p-value for change in pi by region
####################################################

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

png(height=5, width=6, units='in', res=300, file='figures/figureS9.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1), mgp = c(4, 1, 0))

ylims <- c(0, dat[, max(c(-log10(tPd.p), -log10(0.05)), na.rm=TRUE)])
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, -log10(tPd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)


dat[pop == 'lof0711', plot(POSgen, -log10(tPd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(POSgen, -log10(tPd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(POSgen, -log10(tPd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)


dev.off()



####################################################
## Fig. S10 Manhattan plot p-value for change in D by region
####################################################

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

png(height=5, width=6, units='in', res=300, file='figures/figureS10.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1), mgp = c(4, 1, 0))

ylims <- c(0, dat[, max(c(-log10(tDd.p), -log10(0.05)), na.rm=TRUE)])
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, -log10(tDd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)


dat[pop == 'lof0711', plot(POSgen, -log10(tDd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(POSgen, -log10(tDd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(POSgen, -log10(tDd.p), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)


dev.off()



#########################################
## Fig. S11 Manhattan plot FSTs by SNP
#########################################

# WFS results, from wfs_nullmodel_analysis.r
wfs <- fread('analysis/wfs_nullmodel_padj.csv.gz')

# 1907-2011-2014
fst <- fread('analysis/gatk.lof07-11-14.weir.fst', header=TRUE)
datLof <- merge(fst, wfs[pop == 'Lof', .(CHROM, POS, pop, p.adj)], by = c('CHROM', 'POS'))

# Canada
fstCan <- fread('analysis/gatk.can.weir.fst', header=TRUE)
datCan <- merge(fstCan, wfs[pop == 'Can', .(CHROM, POS, pop, p.adj)], by = c('CHROM', 'POS'))
nrow(datCan)

# combine
dat <- rbind(datLof, datCan)

# add genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$length)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

dat <- merge(dat, chrmax[, .(CHROM = chr, start)], by = c('CHROM'))
dat[, POSgen := POS + start]
dat[,start := NULL]

# trim out NAs
dat <- dat[!is.na(WEIR_AND_COCKERHAM_FST), ]


# add a vector for color by LG
lgs <- dat[, sort(unique(CHROM))]
dat[,lgcol := lgcolsramp(p.adj, lg = 1, thresh = 0.1)]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := lgcolsramp(p.adj, lg = 2, thresh = 0.1)]

# and for plotting a colorbar
cbar <- data.frame(logx = seq(0, 1.30103, length.out=50))
cbar$x <- 10^(-cbar$logx)
cbar$col1 <- lgcolsramp(cbar$x, lg = 1, thresh = 0.1)
cbar$col2 <- lgcolsramp(cbar$x, lg = 2, thresh = 0.1)

# order
setorder(dat, pop, -p.adj)

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

ylims <- dat[, range(WEIR_AND_COCKERHAM_FST, na.rm=TRUE)]
xlims <- dat[, range(POSgen, na.rm=TRUE)]


png(height=3.5, width=6, units='in', res=300, file='figures/figureS11.png')
par(mfrow = c(2,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1))

dat[pop == 'Can', plot(POSgen, WEIR_AND_COCKERHAM_FST, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression(F[ST]), bty = 'l', cex.lab = 1, cex.axis = 0.8, cex = 0.3, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.5, line = -0.6) # plot x-axis
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)

color.bar(cbar$col1, cbar$x, axis = FALSE, nticks = 5, xposmin =12e6, xposmax = 16e6, yposmin = -0.1, yposmax = 0.3)
color.bar(cbar$col2, cbar$x, cex = 0.5, axis = TRUE, nticks = 5, 
          xposmin =16e6, xposmax = 20e6, yposmin = -0.1, yposmax = 0.3, title = 'p-value', titley = 0.4)


dat[pop == 'Lof', plot(POSgen, WEIR_AND_COCKERHAM_FST, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(F[ST]), bty = 'l', cex.lab = 1, cex.axis = 0.8, cex = 0.3, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.5, line = -0.6) # plot x-axis
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)


dev.off()




####################################################
## Fig. S12 Manhattan plot p-value for pcangsd
####################################################

# read in data: outlier test from pcangsd (GATK nodam2 unlinked sites)
dat <- fread('analysis/pcangsd_outlier.gatk.nodam.unlinked.csv.gz') # output by angsd_pcangsd_plot_selection.r

# LG mid-points and genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$length)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

dat <- merge(dat, chrmax[, .(CHROM = chr, start)], by = c('CHROM'))
dat[, POSgen := POS + start]
dat[,start := NULL]

# add a vector for color by LG
lgs <- dat[, sort(unique(CHROM))]
dat[,lgcol := lgcolsramp(0.2, lg = 1, thresh = 0.1)]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := lgcolsramp(0.2, lg = 2, thresh = 0.1)]

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS12.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1), mgp = c(4, 1, 0))

ylims <- c(0, dat[, max(c(-log10(pfdr), -log10(0.05)), na.rm=TRUE)])
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')


dat[pop == 'lof0711', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0714', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof1114', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')


dev.off()


############################################################
## Fig. S13 Power to detect selection from FST site-shuffle
############################################################
require(ggplot2)

# read in data: outlier test from pcangsd (GATK nodam2 unlinked sites)
dat <- fread('analysis/slim_fst_siteshuffle.summary.csv.gz') # output by slim_fst_siteshuffle_plot.R

# plot
fs13 <- ggplot(dat[comb == 1, ], aes(s, fpl05, group = f, color = as.factor(f))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(~ ne, labeller = label_both) +
  coord_cartesian(ylim = c(0,1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(1, "lines")) +
  labs(y = 'Proportion that\ndetect outlier loci', color = "Initial frequency")
ggsave(plot = fs13, filename = 'figures/figureS13.png', width = 7, height = 2, dpi = 300)



############################################################
## Fig. S14 Power to detect selection from FST site-shuffle
############################################################
require(ggplot2)

# read in data: outlier test from pcangsd (GATK nodam2 unlinked sites)
dat <- fread('analysis/slim_pcangsd.summary.csv.gz') # output by slim_fst_siteshuffle_plot.R

# plot
fs14 <- ggplot(dat[comb == 1, ], aes(s, prop, group = f, color = as.factor(f))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(~ ne, labeller = label_both) +
  coord_cartesian(ylim = c(0,1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(1, "lines")) +
  labs(y = 'Proportion that\ndetect outlier loci', color = "Initial frequency")
ggsave(plot = fs14, filename = 'figures/figureS14.png', width = 7, height = 2, dpi = 300)



######################
## Table S5 Outliers
######################
library(data.table)

# read in data
dat <- fread('tables/outlier_annotation.csv') # from annotate_outliers.r
freqCan40 <- fread('data_31_01_20/Can_40_freq.mafs.gz')
freqCan14 <- fread('data_31_01_20/Can_14_freq.mafs.gz')
freqLof07 <- fread('data_31_01_20/Lof_07_freq.mafs.gz')
freqLof11 <- fread('data_31_01_20/Lof_11_freq.mafs.gz')
freqLof14 <- fread('data_31_01_20/Lof_14_freq.mafs.gz')

# add column for functional location
dat[,annotation:=as.character(NA)]
dat[NearGene != '', annotation:='<25kb from gene']
dat[grepl('CDS', feature), annotation:='coding']
dat[grepl('three_prime_UTR', feature), annotation:="3' UTR"]
dat[grepl('five_prime_UTR', feature), annotation:="5' UTR"]
dat[grepl('mRNA', feature) & !grepl('CDS|three_prime_UTR|five_prime_UTR', feature), annotation:="transcript"]


# Add allele frequencies
dat[, freq := NA_character_]
dat[, sort(unique(comp))] # make sure all cases are handled below
for(i in 1:nrow(dat)){
  if(dat$comp[i] %in% c('can', 'Canada')){
    freq1 <- freqCan40[chromo == dat$CHROM[i] & position == dat$POS[i], round(knownEM, 2)]
    freq2 <- freqCan14[chromo == dat$CHROM[i] & position == dat$POS[i], round(knownEM, 2)]
    dat[i, freq := paste0(freq1, ', ', freq2)]
  }
  if(dat$comp[i] == 'Norway 1907-2011-2014'){
    freq1 <- freqLof07[chromo == dat$CHROM[i] & position == dat$POS[i], round(knownEM, 2)]
    freq2 <- freqLof11[chromo == dat$CHROM[i] & position == dat$POS[i], round(knownEM, 2)]
    freq3 <- freqLof14[chromo == dat$CHROM[i] & position == dat$POS[i], round(knownEM, 2)]
    dat[i, freq := paste0(freq1, ', ', freq2, ', ', freq3)]
  }
}

	
# make Table S5
out <- dat[,.(LG=CHROM, Position=POS, Pop=gsub('can', 'Canada', comp), Freq = freq, Test=paste0(test, ' ', gsub(' p-value', '', testvaltype)), p = signif(testval,2), annotation,
	gene_name=lapply(strsplit(paste(WithinAnno, NearAnno), split=':'), FUN=function(x) return(x[1])))]


# convert to text and turn NA to ''
for (j in names(out)) set(out, j = j, value = as.character(out[[j]]))
for (j in names(out)) set(out, i=which(is.na(out[[j]])), j = j, value = '')

# write out
write.csv(out, file='figures/tableS5.csv', row.names=FALSE)
