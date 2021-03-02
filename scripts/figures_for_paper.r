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
# Part of Fig. 1 and S6: allele frequency change heatmap
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
        legend.position = 'none') +
  coord_cartesian(xlim =c(0, 1), ylim = c(0, 1))
p2 <- ggplot(datLof0711, aes(freq1, freq2)) +
  geom_hex() +
  scale_fill_gradient(low="grey", high=colLof, trans = 'log', breaks = bks, labels = bks) +
  labs(title = 'Norway 07-11', x = 'Historical', y = 'Modern') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        legend.position = 'none') +
  coord_cartesian(xlim =c(0, 1), ylim = c(0, 1))
p3 <- ggplot(datLof0714, aes(freq1, freq2)) +
  geom_hex() +
  scale_fill_gradient(low="grey", high=colLof, trans = 'log', breaks = bks, labels = bks) +
  labs(title = 'Norway 07-14', x = 'Historical', y = 'Modern') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        legend.position = 'none') +
  coord_cartesian(xlim =c(0, 1), ylim = c(0, 1))
p4 <- ggplot(datLof1114, aes(freq1, freq2)) +
  geom_hex() +
  scale_fill_gradient(low="grey", high=colLof, trans = 'log', breaks = bks, labels = bks) +
  labs(title = 'Norway 11-14', x = 'Historical', y = 'Modern') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        legend.position = 'none') +
  coord_cartesian(xlim =c(0, 1), ylim = c(0, 1))

p5 <- grid.arrange(p1, p2, p3, p4, ncol = 4)

ggsave(p5, filename = 'figures/figure1.png', width = 7.5, height = 2)
ggsave(p5, filename = 'figures/figure1.pdf', width = 7.5, height = 2)



########################################
## Fig. 2 Manhattan plot FSTs by region
########################################

# read in data: sliding window fst from ANGSD (GATK nodam2 sites)
dat1 <- fread('analysis/Can_40.Can_14.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1) # output by angsd_fst.sh. skip headers.
dat3 <- fread('analysis/Lof_07.Lof_14.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1)

dat1[, pop := 'can']
dat3[, pop := 'lof0714']

dat <- rbind(dat1, dat3)

# trim to non-overlapping windows
dat <- dat[(midPos - 1) %% 25000 == 0, ]

# read in shared outlier regions
outl <- fread('tables/outlier_annotation.csv') # annotated list of outliers, from annotate_outliers.r
outl[, midPos := as.numeric(midPos) + 1] # add one to match dat centers
outl <- outl[region == 1 & test == '99th percentile across pops', .(CHROM, midPos, outl = 1)] # trim to just regions
dat <- merge(dat, outl, all.x = TRUE)
dat[is.na(outl), outl := 0]

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
dat[,lgcol := cols[1]]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := cols[2]]


### set up plot
adjlet <- -0.14 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- 0.35
cexsc <- 1/5
xlims <- dat[, range(posgen, na.rm=TRUE)]
ylims <- c(0,0.4)

# quartz(height=8, width=6)
#png(height=3.3, width=6, units='in', res=300, file='figures/figure2.png')
pdf(height=3.3, width=6, file='figures/figure2.pdf')
par(mfrow = c(2,1), las=1, mai=c(0.3, 0.6, 0.25, 0.1), omi = c(0.3, 0, 0, 0))


dat[pop == 'can', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, 
                       ylim = ylims, xlab = '', ylab = expression(F[ST]), bty = 'l', 
                       cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
dat[pop == 'can' & outl == 1, points(posgen, fst, col = colout, cex=log(nloci)*cexsc)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8, mgp = c(1, 0, 0)) # plot x-axis
mtext(side=3, 'A. Canada 1940-2013', adj=adjlet, line=linelet, cex=cexlet)

legend('topleft', legend = c(2, 5, 10, 100), pch = 1, pt.cex = log(c(2,5,10,100))*cexsc, title = '# of SNPs', bty = 'n', cex = 0.7)

dat[pop == 'lof0714', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, 
                           ylim = ylims, xlab = '', ylab = expression(F[ST]), bty = 'l', 
                           cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
dat[pop == 'lof0714' & outl == 1, points(posgen, fst, col = colout, cex=log(nloci)*cexsc)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8, mgp = c(1, 0, 0))
mtext(side=3, 'B. Norway 1907-2014', adj=adjlet, line=linelet, cex=cexlet)
mtext(side=1, 'Linkage group', outer = TRUE)

dev.off()

###################################################
## Fig. S5: read balance at heterozygote genotypes
## (check for reference bias)
###################################################
# read in
can13 <- fread('analysis/allele_balance_popCan13.csv.gz', header= FALSE)
can40 <- fread('analysis/allele_balance_popCan40.csv.gz', header= FALSE)
lof07 <- fread('analysis/allele_balance_popLof07.csv.gz', header= FALSE)
lof11 <- fread('analysis/allele_balance_popLof11.csv.gz', header= FALSE)
lof14 <- fread('analysis/allele_balance_popLof14.csv.gz', header= FALSE)

# plot
linex <- 2
cexx <- 0.6

png(height=4, width=6.5, units='in', res=300, file='figures/figureS5.png')
par(mfrow = c(2,3), las=1, mai=c(0.5, 0.7, 0.25, 0.1), mgp = c(4,0.8,0))

hist(lof07$V1, main = 'Norway 1907')
mtext('REF allele read proportion', 1, cex = cexx, line = linex)

hist(lof11$V1, main = 'Norway 2011')
mtext('REF allele read proportion', 1, cex = cexx, line = linex)

hist(lof14$V1, main = 'Norway 2014')
mtext('REF allele read proportion', 1, cex = cexx, line = linex)

hist(can40$V1, main = 'Canada 1940')
mtext('REF allele read proportion', 1, cex = cexx, line = linex)

hist(can13$V1, main = 'Canada 2013')
mtext('REF allele read proportion', 1, cex = cexx, line = linex)

dev.off()


####################################################
## Fig. S7 Manhattan plot FSTs 1907-2011 and 2011-2014
####################################################

# read in data: sliding window fst from ANGSD (GATK nodam2 sites)
dat2 <- fread('analysis/Lof_07.Lof_11.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1)
dat4 <- fread('analysis/Lof_11.Lof_14.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1)

dat2[, pop := 'lof0711']
dat4[, pop := 'lof1114']

dat <- rbind(dat2, dat4)

# trim to non-overlapping windows
dat <- dat[(midPos - 1) %% 25000 == 0, ]

# read in outlier regions
outl <- fread('tables/outlier_annotation.csv') # annotated list of outliers, from annotate_outliers.r
outl[, midPos := as.numeric(midPos) + 1] # add one to match dat centers
outl <- outl[region == 1, .(CHROM, midPos, outl = 1)] # trim to just regions
dat <- merge(dat, outl, all.x = TRUE)
dat[is.na(outl), outl := 0]

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
dat[,lgcol := cols[1]]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := cols[2]]


### set up plot
adjlet <- -0.14 # horizontal adjustment for subplot letter
cexlet <- 0.7
linelet <- 0.2
cexsc <- 1/5

# quartz(height=8, width=6)
png(height=3, width=6, units='in', res=300, file='figures/figureS7.png')
par(mfrow = c(2,1), las=1, mai=c(0.3, 0.6, 0.3, 0.1))

xlims <- dat[, range(posgen, na.rm=TRUE)]



dat[pop == 'lof0711', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
dat[pop == 'lof0711' & outl == 1, points(posgen, fst, col = colout, cex=log(nloci)*cexsc)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8, mgp = c(1, 0, 0))
mtext(side=3, 'A. Norway 1907-2011', adj=adjlet, line=linelet, cex=cexlet)

legend('topleft', legend = c(2, 5, 10, 100), pch = 1, pt.cex = log(c(2,5,10,100))*cexsc, title = '# of SNPs', bty = 'n', cex = 0.5)

dat[pop == 'lof1114', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8, mgp = c(1, 0, 0))
mtext(side=3, 'B. Norway 2011-2014', adj=adjlet, line=linelet, cex=cexlet)


dev.off()


####################################################
## Fig. S8 Manhattan plot FST site-shuffle p-value
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
adjlet <- -0.14 # horizontal adjustment for subplot letter
cexlet <- 0.8
linelet <- 0.2
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS8.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.6, 0.2, 0.1))

ymin <- min(c(dat[, min(p, na.rm=TRUE)], 0.05))
xlims <- dat[, range(posgen, na.rm=TRUE)]

dat[pop == 'can', plot(posgen, p, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                       ylim = c(1, ymin), ylab = 'p', bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'A. Canada 1940-2013', adj=adjlet, line=linelet, cex=cexlet)

legend('topleft', legend = c(2, 5, 10, 100), pch = 1, pt.cex = log(c(2,5,10,100))*cexsc, title = '# of SNPs', bty = 'n', cex = 0.5)


dat[pop == 'lof0711', plot(posgen, p, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                           ylim = c(1, ymin), ylab = 'p', bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'B. Norway 1907-2011', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(posgen, p, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                           ylim = c(1, ymin), ylab = 'p', bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'C. Norway 1907-2014', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(posgen, p, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', 
                           ylim = c(1, ymin), ylab = 'p', bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8) # plot x-axis
mtext(side=3, 'D. Norway 2011-2014', adj=adjlet, line=linelet, cex=cexlet)


dev.off()


####################################################
## Fig. S9 Manhattan plot change in pi by region
####################################################
winsz = 5e4 # for scaling the windowed pi and thetaW values

# read in data: sliding window pi from ANGSD (GATK nodam2 unlinked sites)
dat1 <- fread('analysis/theta_change_region_50000.Can_40.Can_14.gatk.csv.gz', drop = 1) # from angsd_theta_siteshuffle_null.r
dat2 <- fread('analysis/theta_change_region_50000.Lof_07.Lof_11.gatk.csv.gz', drop = 1)
dat3 <- fread('analysis/theta_change_region_50000.Lof_07.Lof_14.gatk.csv.gz', drop = 1)
dat4 <- fread('analysis/theta_change_region_50000.Lof_11.Lof_14.gatk.csv.gz', drop = 1)

dat1[, pop := 'can']
dat2[, pop := 'lof0711']
dat3[, pop := 'lof0714']
dat4[, pop := 'lof1114']

dat <- rbind(dat1, dat2, dat3, dat4)

# LG mid-points and genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

dat <- merge(dat[, .(pop, CHROM = Chromo, midPos = WinCenter, tPd)], chrmax[, .(CHROM = chr, start)])
dat[, POSgen := midPos + start]
dat[, start := NULL]

# add a vector for color by LG
lgs <- dat[, sort(unique(CHROM))]
dat[,lgcol := cols[1]]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := cols[2]]

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS9.png')
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
## Fig. S10 Manhattan plot change in D by region
####################################################

# read in data: sliding window D from ANGSD (GATK nodam2 unlinked sites)
dat1 <- fread('analysis/theta_change_region_50000.Can_40.Can_14.gatk.csv.gz', drop = 1) # from angsd_theta_siteshuffle_null.r
dat2 <- fread('analysis/theta_change_region_50000.Lof_07.Lof_11.gatk.csv.gz', drop = 1)
dat3 <- fread('analysis/theta_change_region_50000.Lof_07.Lof_14.gatk.csv.gz', drop = 1)
dat4 <- fread('analysis/theta_change_region_50000.Lof_11.Lof_14.gatk.csv.gz', drop = 1)

dat1[, pop := 'can']
dat2[, pop := 'lof0711']
dat3[, pop := 'lof0714']
dat4[, pop := 'lof1114']

dat <- rbind(dat1, dat2, dat3, dat4)

# LG mid-points and genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

dat <- merge(dat[, .(pop, CHROM = Chromo, midPos = WinCenter, tDd)], chrmax[, .(CHROM = chr, start)])
dat[, POSgen := midPos + start]
dat[, start := NULL]

# add a vector for color by LG
lgs <- dat[, sort(unique(CHROM))]
dat[,lgcol := cols[1]]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := cols[2]]

### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS10.png')
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
## Fig. S11 Manhattan plot change in LD by region
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

png(height=5, width=6, units='in', res=300, file='figures/figureS11.png')
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
## Fig. S12 Manhattan plot p-value for change in pi by region
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
adjlet <- -0.14 # horizontal adjustment for subplot letter
cexlet <- 0.8
linelet <- 0.2
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS12.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.6, 0.2, 0.1))

ylims <- c(1, dat[, min(c(tPd.p, 0.05), na.rm=TRUE)])
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, tPd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A. Canada 1940-2013', adj=adjlet, line=linelet, cex=cexlet)


dat[pop == 'lof0711', plot(POSgen, tPd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B. Norway 1907-2011', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(POSgen, tPd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C. Norway 1907-2014', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(POSgen, tPd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D. Norway 2011-2014', adj=adjlet, line=linelet, cex=cexlet)


dev.off()



####################################################
## Fig. S13 Manhattan plot p-value for change in D by region
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
adjlet <- -0.14 # horizontal adjustment for subplot letter
cexlet <- 0.8
linelet <- 0.2
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS13.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.6, 0.2, 0.1))

ylims <- c(1, dat[, min(c(tDd.p, 0.05), na.rm=TRUE)])
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, tDd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A. Canada 1940-2013', adj=adjlet, line=linelet, cex=cexlet)


dat[pop == 'lof0711', plot(POSgen, tDd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B. Norway 1907-2011', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(POSgen, tDd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C. Norway 1907-2014', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(POSgen, tDd.p, type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = 'p', bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
abline(h = 0.05, lty = 2, col = 'grey')
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D. Norway 2011-2014', adj=adjlet, line=linelet, cex=cexlet)


dev.off()


############################################################
## Fig. S14 Power to detect selection from FST site-shuffle
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
ggsave(plot = fs13, filename = 'figures/figureS14.png', width = 7, height = 2, dpi = 300)


#########################################
# Fig. S15 Fst around putative outliers
#########################################
library(data.table)
ncol = 3 # number of columns in graph
rowinch <- 2 # inches per row in graph
rng <- 100000 # how many bp to go in either direction from SNP

# read in outliers
dat <- fread('tables/outlier_annotation.csv') # annotated list of outliers, from annotate_outliers.r
dat[, midPos := as.numeric(dat$midPos)]

# read in Fsts
fstLof1 <- fread('analysis/Lof_07.Lof_11.fst.AB.gz', col.names = c('CHROM', 'POS', 'A', 'B')) # from angsd_fst.sh
fstLof2 <- fread('analysis/Lof_07.Lof_14.fst.AB.gz', col.names = c('CHROM', 'POS', 'A', 'B'))
fstCan <- fread('analysis/Can_40.Can_14.fst.AB.gz', col.names = c('CHROM', 'POS', 'A', 'B'))

# trim fsts to nodam2 loci
gatk <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab', col.names = c('CHROM', 'POS', 'REF', 'ALT'))
fstLof1 <- merge(fstLof1, gatk[, .(CHROM, POS)])
fstLof2 <- merge(fstLof2, gatk[, .(CHROM, POS)])
fstCan <- merge(fstCan, gatk[, .(CHROM, POS)])

# calc and merge fsts
fstLof1[, fst := A/B]
fstLof2[, fst := A/B]
fstLof <- merge(fstLof1[, .(CHROM, POS, fst1 = A/B)], fstLof2[, .(CHROM, POS, fst2 = A/B)])
fstLof[, fst := rowMeans(cbind(fst1, fst2))]
fstCan[, fst := A/B]

# read in and merge sample size data
datCan40 <- fread('data_31_01_20/Can_40_freq.mafs.gz')
datCan14 <- fread('data_31_01_20/Can_14_freq.mafs.gz')
datLof07 <- fread('data_31_01_20/Lof_07_freq.mafs.gz')
datLof11 <- fread('data_31_01_20/Lof_11_freq.mafs.gz')
datLof14 <- fread('data_31_01_20/Lof_14_freq.mafs.gz')

datCan <- merge(datCan40[, .(CHROM = chromo, POS = position, n1 = nInd)], datCan14[, .(CHROM = chromo, POS = position, n2 = nInd)])
datLof <- merge(datLof07[, .(CHROM = chromo, POS = position, n1 = nInd)], datLof11[, .(CHROM = chromo, POS = position, n2 = nInd)])
datLof <- merge(datLof, datLof14[, .(CHROM = chromo, POS = position, n3 = nInd)])

# harmonic mean sample size
datCan[, n := 2/(1/n1 + 1/n2)]
datLof[, n12 := 2/(1/n1 + 1/n2)]
datLof[, n13 := 2/(1/n1 + 1/n3)]
datLof[, n123 := 3/(1/n1 + 1/n2 + 1/n3)]

datCan[, n.sc := (n - min(n))/(max(n) - min(n))] # scale 0-1
datLof[, n12.sc := (n12 - min(n12))/(max(n12) - min(n12))] # scale 0-1
datLof[, n13.sc := (n13 - min(n13))/(max(n13) - min(n13))] # scale 0-1
datLof[, n123.sc := (n123 - min(n123))/(max(n123) - min(n123))] # scale 0-1

fstCan <- merge(fstCan, datCan[, .(CHROM, POS, n, n.sc)], by = c('CHROM', 'POS'))
fstLof1 <- merge(fstLof1, datLof[, .(CHROM, POS, n12, n12.sc)], by = c('CHROM', 'POS'))
fstLof2 <- merge(fstLof2, datLof[, .(CHROM, POS, n13, n13.sc)], by = c('CHROM', 'POS'))
fstLof <- merge(fstLof, datLof[, .(CHROM, POS, n12, n13, n123, n12.sc, n13.sc, n123.sc)], by = c('CHROM', 'POS'))

# merge adjacent regions of the same comparison & test
dat2 <- dat # the version to modify
dat2[, todelete := 0]
for(i in 1:nrow(dat)){
  j <- dat[, which(CHROM == dat$CHROM[i] & comp == dat$comp[i] & test == dat$test[i] & abs(midPos - dat$midPos[i]) < 100000)]
  if(length(j) > 1){
    print(i)
    dat2$POS[j[1]] <- paste0(dat$POS[j], collapse = ',')
    dat2$midPos[j[1]] <- mean(dat$midPos[j])
    dat2$todelete[j[-1]] <- 1
  }
}
dat2[, .(CHROM, POS, midPos, comp, test, todelete)]
dat2 <- dat2[todelete == 0, ]

# make plot
nrow = ceiling(nrow(dat)/ncol)
bty <- 'l'
outlcol <- 'red'

png(height= 9, width=6.5, units='in', res=300, file='figures/figureS15.png')
par(mfrow = c(nrow, ncol), mai = c(0.3, 0.3, 0.4, 0.1), omi = c(0.3, 0.3, 0, 0))
for(i in 1:nrow(dat2)){
  if(dat2$comp[i] %in% c('can', 'Canada')){
    thisdat <- fstCan[CHROM == dat2$CHROM[i] & abs(POS - dat2$midPos[i]) < rng, ]
    outl <- fstCan[CHROM == dat2$CHROM[i] & POS == dat2$midPos[i], ]
    if(grepl(',', dat2$POS[i])) outl <- fstCan[CHROM == dat2$CHROM[i] & POS %in% unlist(strsplit(dat2$POS[i], split = ',')), ] # if a joined region
    pop <- 'Canada'
    thisdat[, plot(POS/1e6, fst, cex = n.sc, xlab = '', ylab = '',
                   main = paste0(dat2$CHROM[i], ' ', pop), bty = bty)]
    outl[, points(POS/1e6, fst, col = outlcol)]
  }
  if(dat2$comp[i] == 'Norway 1907-2011'){
    thisdat <- fstLof1[CHROM == dat2$CHROM[i] & abs(POS - dat2$midPos[i]) < rng, ]
    outl <- fstLof1[CHROM == dat2$CHROM[i] & POS == dat2$midPos[i], ]
    if(grepl(',', dat2$POS[i])) outl <- fstCan[CHROM == dat2$CHROM[i] & POS %in% unlist(strsplit(dat2$POS[i], split = ',')), ] # if a joined region
    pop <- 'Norway 1907-2011'
    thisdat[, plot(POS/1e6, fst, cex = n12.sc, xlab = '', ylab = '',
                   main = paste0(dat2$CHROM[i], ' ', pop), bty = bty)]
    outl[, points(POS/1e6, fst, col = outlcol)]
  }
  if(dat2$comp[i] == 'Norway 1907-2014'){
    thisdat <- fstLof2[CHROM == dat2$CHROM[i] & abs(POS - dat2$midPos[i]) < rng, ]
    outl <- fstLof2[CHROM == dat2$CHROM[i] & POS == dat2$midPos[i], ]
    if(grepl(',', dat2$POS[i])) outl <- fstCan[CHROM == dat2$CHROM[i] & POS %in% unlist(strsplit(dat2$POS[i], split = ',')), ] # if a joined region
    pop <- 'Norway 1907-2014'
    thisdat[, plot(POS/1e6, fst, cex = n13.sc, xlab = '', ylab = '',
                   main = paste0(dat2$CHROM[i], ' ', pop), bty = bty)]
    outl[, points(POS/1e6, fst, col = outlcol)]
  }
  if(dat2$comp[i] == 'Norway 1907-2011-2014'){
    thisdat <- fstLof[CHROM == dat2$CHROM[i] & abs(POS - dat2$midPos[i]) < rng, ]
    outl <- fstLof[CHROM == dat2$CHROM[i] & POS == dat2$midPos[i], ]
    if(grepl(',', dat2$POS[i])) outl <- fstCan[CHROM == dat2$CHROM[i] & POS %in% unlist(strsplit(dat2$POS[i], split = ',')), ] # if a joined region
    pop <- 'Norway 1907-2011-2014'
    thisdat[, plot(POS/1e6, fst, cex = n123.sc, xlab = '', ylab = '',
                   main = paste0(dat2$CHROM[i], ' ', pop), cex.main = 0.95, bty = bty)]
    outl[, points(POS/1e6, fst, col = outlcol)]
  }
  if(dat2$comp[i] == 'Canada-Norway'){
    thisdat1 <- fstLof[CHROM == dat2$CHROM[i] & abs(POS - dat2$midPos[i]) < rng, ]
    thisdat2 <- fstCan[CHROM == dat2$CHROM[i] & abs(POS - dat2$midPos[i]) < rng, ]
    pop <- 'All'
    thisdat1[, plot(POS/1e6, fst, cex = n123.sc, xlab = '', ylab = '',
                    main = paste0(dat2$CHROM[i], ' ', pop), bty = bty)]
    thisdat2[, points(POS/1e6, fst, cex = n.sc, col = 'grey')]
  }
  title(main = dat2$test[i], line = 0.4, cex.main = 0.7)
}

mtext('Position (Mb)', side = 1, outer = TRUE)
mtext('Fst', side = 2, outer = TRUE)

dev.off()



#########################################
## Fig. S16 Manhattan plot FSTs by SNP
#########################################

# WFS results, from wfs_nullmodel_analysis.r
wfs <- fread('analysis/wfs_nullmodel_padj.csv.gz')

# 1907-2011-2014 fsts
fst <- fread('analysis/gatk.lof07-11-14.weir.fst', header=TRUE) # from vcftools --gzvcf data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz --weir-fst-pop data_2020.05.07/popLof07.txt --weir-fst-pop data_2020.05.07/popLof11.txt --weir-fst-pop data_2020.05.07/popLof14.txt --out analysis/gatk.lof07-11-14
datLof <- merge(fst, wfs[pop == 'Lof', .(CHROM, POS, pop, p.adj)], by = c('CHROM', 'POS'))

# Canada fsts
fstCan <- fread('analysis/gatk.can.weir.fst', header=TRUE) # from vcftools --gzvcf data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz --weir-fst-pop data_2020.05.07/popCan40.txt --weir-fst-pop data_2020.05.07/popCan13.txt --out analysis/gatk.can
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


png(height=3.5, width=6, units='in', res=300, file='figures/figureS16.png')
par(mfrow = c(2,1), las=1, mai=c(0.3, 0.8, 0.1, 0.1))

dat[pop == 'Can', plot(POSgen, WEIR_AND_COCKERHAM_FST, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression(F[ST]), bty = 'l', cex.lab = 1, cex.axis = 0.8, cex = 0.3, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.5, line = -0.6) # plot x-axis
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)

color.bar(cbar$col1, cbar$x, axis = FALSE, nticks = 5, xposmin =20e6, xposmax = 26e6, yposmin = -0.95, yposmax = -0.2)
color.bar(cbar$col2, cbar$x, cex = 0.7, axis = TRUE, nticks = 5, 
          xposmin =26e6, xposmax = 32e6, yposmin = -0.95, yposmax = -0.2, title = 'p-value', titley = -0.05)


dat[pop == 'Lof', plot(POSgen, WEIR_AND_COCKERHAM_FST, type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression(F[ST]), bty = 'l', cex.lab = 1, cex.axis = 0.8, cex = 0.3, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.5, line = -0.6) # plot x-axis
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)


dev.off()



####################################################
## Fig. S17 Manhattan plot p-value for pcangsd
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
adjlet <- -0.14 # horizontal adjustment for subplot letter
cexlet <- 0.8
linelet <- 0.2
cexsc <- 1/5

png(height=5, width=6, units='in', res=300, file='figures/figureS17.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.6, 0.2, 0.1))

ylims <- c(0, dat[, max(c(-log10(pfdr), -log10(0.05)), na.rm=TRUE)])
xlims <- dat[, range(POSgen, na.rm=TRUE)]

dat[pop == 'can', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                       ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'A. Canada 1940-2013', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')


dat[pop == 'lof0711', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'B. Norway 1907-2011', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0714', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'C. Norway 1907-2014', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof1114', plot(POSgen, -log10(pfdr), type='p', col=lgcol, xlim = xlims, xlab = '', 
                           ylim = ylims, ylab = expression(-log[10](p)), bty = 'l', cex.lab = 1, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE, cex.axis = 0.8)
mtext(side=3, 'D. Norway 2011-2014', adj=adjlet, line=linelet, cex=cexlet)
abline(h = -log10(0.05), lty = 2, col = 'grey')


dev.off()


############################################################
## Fig. S18 Power to detect selection from PCAngsd
############################################################
require(ggplot2)

# read in data: outlier test from pcangsd (GATK nodam2 unlinked sites)
dat <- fread('analysis/slim_pcangsd.summary.csv.gz') # output by slim_pcangsdoutlier_plot.R

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
ggsave(plot = fs14, filename = 'figures/figureS18.png', width = 7, height = 2, dpi = 300)



######################
## Table S5 Outliers
######################
library(data.table)

# read in data
dat <- fread('tables/outlier_annotation.csv') # annotated list of outliers, from annotate_outliers.r
dat[, midPos := as.integer(midPos)]
freqCan40 <- fread('data_31_01_20/Can_40_freq.mafs.gz') # allele frequencies
freqCan14 <- fread('data_31_01_20/Can_14_freq.mafs.gz') # actually 2013, not 2014
freqLof07 <- fread('data_31_01_20/Lof_07_freq.mafs.gz')
freqLof11 <- fread('data_31_01_20/Lof_11_freq.mafs.gz')
freqLof14 <- fread('data_31_01_20/Lof_14_freq.mafs.gz')

setkey(freqCan40, chromo, position)
setkey(freqCan14, chromo, position)
setkey(freqLof07, chromo, position)
setkey(freqLof11, chromo, position)
setkey(freqLof14, chromo, position)

# p-value for 99th percentile is an fst, not a p-value, so erase it
dat[testvaltype == 'Average Fst', testval := NA_real_]
dat[testvaltype == 'Average Fst', testvaltype := '']

# adjust text slightly
dat[test == '99th percentile across pops', test := 'Shared 99th percentile']
dat[, test := gsub('pcangsd', 'PCAngsd', test)]
dat[testvaltype == 'Genome-wide p-value for the region', testvaltype := 'Genome-wide p']

# add column for functional location
dat[,annotation:='']
dat[NearGene != '', annotation:= paste(annotation, '<25kb from gene', sep = ', ')]
dat[grepl('CDS', feature), annotation:='coding']
dat[grepl('three_prime_UTR', feature), annotation:= paste(annotation, "3' UTR", sep = ', ')]
dat[grepl('five_prime_UTR', feature), annotation:= paste(annotation, "5' UTR", sep = ', ')]
dat[grepl('mRNA', feature) & !grepl('CDS|three_prime_UTR|five_prime_UTR', feature), annotation:= paste(annotation, "transcript", sep = ', ')]
dat[, annotation := gsub('^, ', '', annotation)]

# Add allele frequencies
dat[, Canada1940 := NA_character_]
dat[, Canada2013 := NA_character_]
dat[, Norway1907 := NA_character_]
dat[, Norway2011 := NA_character_]
dat[, Norway2014 := NA_character_]
dat[, sort(unique(comp))] # make sure all cases are handled below
for(i in 1:nrow(dat)){
  if(dat$region[i] == 0){ # can only do this for snps
    if(dat$comp[i] %in% c('can', 'Canada', 'Canada-Norway')){
      dat$Canada1940[i] <- freqCan40[.(dat$CHROM[i], dat$midPos[i]), round(knownEM, 2)]
      dat$Canada2013[i] <- freqCan14[.(dat$CHROM[i], dat$midPos[i]), round(knownEM, 2)]
    }
    if(dat$comp[i] %in% c('Norway 1907-2011-2014', 'Norway 1907-2011', 'Canada-Norway', 'Norway 1907-2014')){
      dat$Norway1907[i] <- freqLof07[.(dat$CHROM[i], dat$midPos[i]), round(knownEM, 2)]
    }
    if(dat$comp[i] %in% c('Norway 1907-2011-2014', 'Norway 1907-2011', 'Canada-Norway')){
      dat$Norway2011[i] <- freqLof11[.(dat$CHROM[i], dat$midPos[i]), round(knownEM, 2)]
    }
    if(dat$comp[i] %in% c('Norway 1907-2011-2014', 'Norway 1907-2014', 'Canada-Norway')){
      dat$Norway2014[i] <- freqLof14[.(dat$CHROM[i], dat$midPos[i]), round(knownEM, 2)]
    }
  }
}

	
# make Table S5
out <- dat[,.(LG=CHROM, Position=POS, Pop=gsub('can', 'Canada', comp), Canada1940, Canada2013, Norway1907, Norway2011, Norway2014, Test=paste0(test, ' ', gsub(' p-value', '', testvaltype)), p = signif(testval,2), annotation,
	Gene=lapply(strsplit(paste(WithinAnno, NearAnno), split=':'), FUN=function(x) return(x[1])))]


# convert to text and turn NA to ''
for (j in names(out)) set(out, j = j, value = as.character(out[[j]]))
for (j in names(out)) set(out, i=which(is.na(out[[j]])), j = j, value = '')

# write out
write.csv(out, file='tables/tableS5.csv', row.names=FALSE)
