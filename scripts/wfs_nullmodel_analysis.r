# Examine the probabilities of null model producing results as extreme as our observations
# run after wfs_nullmodel_function.r and wfs_nullmodel_combine.r/wfs_nullmodel_lowp_rerun.r
# Script looks at the output
require(data.table)
require(ggplot2)

#####################
## Functions
#####################

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


##########
## Prep 
##########

# 1907-2011-2014
fst <- fread('analysis/gatk.lof07-11-14.weir.fst', header=TRUE)
dat <- readRDS(file=paste('analysis/wfs_nullmodel_pos&pvals_07-11-14.rds', sep='')) # has p-values
nrow(dat)
dat <- merge(fst, dat[, .(CHROM, POS, p)], by = c('CHROM', 'POS'))
nrow(dat)

# Canada
fstCan <- fread('analysis/gatk.can.weir.fst', header=TRUE)
datCan <- readRDS(file=paste('analysis/wfs_nullmodel_pos&pvals_Can.rds', sep='')) # has p-values
nrow(datCan)
datCan <- merge(fstCan, datCan[, .(CHROM, POS, p)], by = c('CHROM', 'POS'))
nrow(datCan)

# combine
dat[, pop := 'Lof']
datCan[, pop := 'Can']
dat <- rbind(dat, datCan)

# add genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]

dat <- merge(dat, chrmax[, .(CHROM = chr, start)], by = c('CHROM'))
dat[, POSgen := POS + start]
dat[,start := NULL]

# trim to nodam2 locus set
nodam2 <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab')
setnames(nodam2, c('CHROM', 'POS', 'REF', 'ALT'))
dim(dat)
dat <- merge(dat, nodam2[, .(CHROM, POS)], by = c('CHROM', 'POS'))
dim(nodam2)
dim(dat)

# re-calculate FDR
dat[, p.adj := p.adjust(p), by = pop]



##################
## basic analysis
##################
dat[,min(p, na.rm=TRUE), by = pop]
dat[,min(p.adj, na.rm=TRUE), by = pop]

# lowest FDR-correct p-values
dat[pop == 'Can' & p.adj < 0.65, ]
dat[pop == 'Lof' & p.adj < 0.63, ]

#####################
## plots
#####################
# add a vector for color by LG
lgs <- dat[, sort(unique(CHROM))]
dat[,lgcol := lgcolsramp(p.adj, lg = 1, thresh = 0.1)]
dat[CHROM %in% lgs[seq(2, length(lgs),by=2)], lgcol := lgcolsramp(p.adj, lg = 2, thresh = 0.1)]

# order
setorder(dat, pop, -p.adj)

# histogram of p-values
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_hist_pvals.png', units='in', res=300)
ggplot(dat, aes(p, group = pop)) +
  geom_histogram() +
  facet_grid(~pop)

	dev.off()

# histogram of FDR-adjusted p-values
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_hist_padjvals.png', units='in', res=300)
	ggplot(dat, aes(p.adj, group = pop)) +
	  geom_histogram() +
	  facet_grid(~pop)
	
	dev.off()



# Per locus FST vs. genome position, colored by -log10(p)
# for plotting a colorbar
cbar <- data.frame(logx = seq(0, 1.30103, length.out=50))
cbar$x <- 10^(-cbar$logx)
cbar$col1 <- lgcolsramp(cbar$x, lg = 1, thresh = 0.1)
cbar$col2 <- lgcolsramp(cbar$x, lg = 2, thresh = 0.1)
	
	
png(filename = paste0('figures/wfs_nullmodel_fst_vs_position_by_p.png'), width = 20, height = 12, units = "in", res = 150)
par(mfrow = c(2,1))
dat[pop == 'Can' & !is.na(WEIR_AND_COCKERHAM_FST), plot(POSgen/1e6, WEIR_AND_COCKERHAM_FST, xlab = 'Position (Mb)', ylab = "FST",
                       main = 'Canada', col = lgcol, cex = 0.5)]

dat[pop == 'Lof' & !is.na(WEIR_AND_COCKERHAM_FST), plot(POSgen/1e6, WEIR_AND_COCKERHAM_FST, xlab = 'Position (Mb)', ylab = "FST",
                                                        main = 'Lof', col = lgcol, cex = 0.5)]

color.bar(cbar$col1, cbar$x, axis = FALSE, xposmin =640, xposmax = 645, 
          yposmin = 0.2, yposmax = 0.5)
color.bar(cbar$col2, cbar$x, cex = 0.5, axis = TRUE, nticks = 5, 
          xposmin =645, xposmax = 650, yposmin = 0.2, yposmax = 0.5, title = 'FDR p-value', titley = 0.55)

dev.off()
	


# -log10(p.adj) vs. genome position
	col='red'
	quartz(width=8, height=6)
	# png(width=8, height=6, file=paste('figures/wfs_nullmodel_padj_vs_position', suffix, '.png', sep=''), res=300, units='in')
	plot(dat$POSgen, -log10(dat$p.adj),type='p', cex=0.2, xlab='Genome position', 
	     ylab='-log10(FDR-adjusted p)', col=rgb(0.1,0.1,0.1, 0.2), ylim = c(0,1.35))
		abline(h=-log10(0.05), col='red', lty=2)

		# add LG labels
	lgs <- sort(unique(dat[,CHROM]))
	for(j in 1:length(lgs)){
		rng <- range(dat[CHROM==lgs[j], POSgen])
		if(j %% 2 == 0) lines(x=rng, y=c(0,0), col=col, lwd=2)
		if(j %% 2 == 1) lines(x=rng, y=c(0.015,0.015), col=col, lwd=2)
		text(x=mean(rng), y=0.03, labels=lgs[j], col=col, cex=0.5)
	
	}

		# add mean p
	lines(meanp$mids, meanp$p.adj, col='orange')

	dev.off()

