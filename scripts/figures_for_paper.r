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
colout <- '#33a02c' # green, part of Paired colorbrewer, for outlier loci
#cols3 <- c('#F7941E', '#F2592A', '#BF1E2D', '#FFDE17', '#BEA512') # Lof07, Lof11, Lof14, Can40, CanMod (reds, yellows, from Bastiaan)
cols3 <- c('#e7d4e8', '#af8dc3', '#762a83', '#d9f0d3', '#1b7837') # Lof07, Lof11, Lof14, Can40, CanMod (PRGn colorbrewer)

##################################
# Fig. 2 Fishing and phenotypes
##################################

# read in data
canfish <- fread('data/phenotypes/Brattey_etal_2018_CSAS_TableA2-5_NorthernCod_fishingmortality.csv')
loffish <- fread('data/phenotypes/AFWG_2019_3_Northeast Arctic_Cod_Table3.18.csv')
canmat <- fread('output/age_50percmature_can.csv')
lofmat <- fread('data/phenotypes/Eikeset_Age_and_length_at_maturation.csv')

# summarize for plotting
canfishsum <- canfish[, .(Year, F = rowMeans(cbind(Age5, Age6, Age7, Age8, Age9, Age10)))]

# plot params
adjlet <- -0.2
cexlet <- 1
linelet <- 0.5

# plot
png(height=6, width=6, units='in', res=300, file='figures/figure2.png')

par(mfrow=c(2,2), mai = c(0.6, 0.7, 0.3, 0.1), las = 1, tcl = -0.3, mgp = c(2.5, 0.7, 0), bty = 'l')
canfishsum[, plot(Year, F, type = 'l', ylab = 'Fishing mortality rate')]
mtext(side=3, 'A)', adj=adjlet, line=linelet, cex=cexlet)

loffish[, plot(Year, `FBAR 5-10`, type = 'l', ylab = 'Fishing mortality rate')]
mtext(side=3, 'B)', adj=adjlet, line=linelet, cex=cexlet)

canmat[, plot(Year, age50, type = 'l', ylab = 'Age at 50% maturity')]
mtext(side=3, 'C)', adj=adjlet, line=linelet, cex=cexlet)

lofmat[, plot(yearData, DataAgeMaturation, type = 'l', xlab = 'Year', ylab = 'Age at 50% maturity')]
mtext(side=3, 'D)', adj=adjlet, line=linelet, cex=cexlet)

dev.off()

###############################
## Fig. 3 Manhattan plot FSTs
###############################

# read in data: sliding window fst and site-shuffle p-values from ANGSD (GATK sites)
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
dat[,lgcol := cols[1]]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := cols[2]]
dat[,lgcol2 := cols2[1]]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol2 := cols2[2]]


### set up plot
adjlet <- -0.11 # horizontal adjustment for subplot letter
cexlet <- 1
linelet <- -0.5
cexsc <- 1/5

# quartz(height=8, width=6)
png(height=5, width=6, units='in', res=300, file='figures/figure3.png')
par(mfrow = c(4,1), las=1, mai=c(0.3, 0.6, 0.1, 0.1))

ymax <- dat[, max(fst, na.rm=TRUE)]
xlims <- dat[, range(posgen, na.rm=TRUE)]

dat[pop == 'can', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
dat[pop == 'can' & p < 0.05, points(posgen, fst, type='p', cex=1, pch = 16, col=colout)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'A', adj=adjlet, line=linelet, cex=cexlet)
legend('topright', legend = c(2, 5, 10, 100), pch = 1, pt.cex = log(c(2,5,10,100))*cexsc, title = 'Number of loci', bty = 'n', cex = 0.5)

dat[pop == 'lof0711', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
dat[pop == 'lof0711' & p < 0.05, points(posgen, fst, type='p', cex=1, pch = 16, col=colout)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'B', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof0714', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
dat[pop == 'lof0714' & p < 0.05, points(posgen, fst, type='p', cex=1, pch = 16, col=colout)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)

dat[pop == 'lof1114', plot(posgen, fst, type='p', cex=log(nloci)*cexsc, col=lgcol, xlim = xlims, xlab = '', ylab = expression(F[ST]), bty = 'l', cex.lab = 1.5, xaxt = 'n', xaxs = 'i', tcl = -0.3)]
dat[pop == 'lof1114' & p < 0.05, points(posgen, fst, type='p', cex=1, pch = 16, col=colout)]
axis(side=1, at = chrmax$mid, labels = gsub('LG|LG0', '', chrmax$chr), tick = FALSE) # plot x-axis
mtext(side=3, 'C', adj=adjlet, line=linelet, cex=cexlet)


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