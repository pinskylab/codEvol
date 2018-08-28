# analyze h and s posteriors from ABC

require(data.table)
require(MASS)
require(png)

# read in data
posts <- readRDS('analysis/wfs_abc_posts_Lof_07-1114.rds'); pop<-'Lof'

# sort
setkey(posts, CHROM, POS)

# plot single loci
posts[CHROM=='LG03' & POS=='8094943', plot(h, s, col='#00000011', type='p', cex=0.3)]

posts[CHROM=='LG03' & POS=='21120013', plot(h, s, col='#00000011', type='p', cex=0.3)]

posts[CHROM=='LG04' & POS=='2924609', plot(h, s, col='#00000011', type='p', cex=0.3)]

posts[CHROM=='LG09' & POS=='3062018', plot(h, s, col='#00000011', type='p', cex=0.3)] # different pattern

posts[CHROM=='LG22' & POS=='10459232', plot(h, s, col='#00000011', type='p', cex=0.3)] # different pattern again


# plot all loci as a panel of plots
n <- sum(!duplicated(posts[,.(CHROM, POS)]))
nc <- ceiling(sqrt(n))
nr <- ceiling(n/nc)

	# as transparent dots
	quartz(width=8, height=6.5)
	par(mfrow=c(nr, nc), mai=c(0.1, 0.15, 0.1, 0.05), tcl=-0.15, mgp=c(2,0.15,0), las=1, omi=c(0.2, 0.2, 0, 0))
	posts[,plot(h, s, col='#00000005', type='p', cex=0.3, xlab='', ylab='', xlim=c(0,1.5), ylim=c(-1,1), main=paste(CHROM, POS), cex.main=0.5, cex.axis=0.5), by=c('CHROM', 'POS')]
	mtext('h', side=1, outer=TRUE)
	mtext('s', side=2, outer=TRUE)
	
	# as 2D posterior
	quartz(width=8, height=6.5)
	# png(width=8, height=6.5, units='in', res=300, file=paste('figures/wfs_abc_posterior_s_vs_h_kde2d_', pop, '.png', sep=''))
	par(mfrow=c(nr, nc), mai=c(0.1, 0.15, 0.1, 0.05), tcl=-0.15, mgp=c(2,0.15,0), las=1, omi=c(0.2, 0.2, 0, 0))
	posts[,image(kde2d(h,s,n=300, lims=c(0,1.5,-1,1)), xlab='', ylab='', main=paste(CHROM, POS), cex.main=0.5, cex.axis=0.5), by=c('CHROM', 'POS')]

	mtext('h', side=1, outer=TRUE)
	mtext('s', side=2, outer=TRUE)

	dev.off()