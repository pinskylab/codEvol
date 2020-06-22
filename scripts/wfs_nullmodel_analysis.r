# Calculate probability of null model producing results as extreme as our observations
# run after wfs_nullmodel_function.r and wfs_nullmodel_combine.r
# Script analyzes the output

##########
## Prep 
##########
require(data.table)

# Choose one set to read in
	# 1907-2011-2014
	locnms <- fread('data_2018.09.05/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	locnms2 <- fread('data_2018.09.05/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(locnms2, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_3', 'Freq_3', 'ABS_DIFF2'))
	setkey(locnms, CHROM, POS, N_CHR_1, Freq_1)
	setkey(locnms2, CHROM, POS, N_CHR_1, Freq_1)
	locnms <- locnms[locnms2,]
	suffix <- '_07-11-14'

	# Canada
	locnms <- fread('data_2019_03_18/Frequency_table_CAN_40_TGA.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	suffix <- '_Can'

# Then continue here:

# Read in p-values by locus
infile <- paste('analysis/wfs_nullmodel_pos&pvals', suffix, '.rds', sep='')
infile
dat <- readRDS(file=infile) # has p-values

# add genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
chrmax[, mid := rowSums(cbind(start, length/2))]

setkey(dat, CHROM)
dat <- dat[chrmax[, .(CHROM = chr, start)], ]
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
dat[, p.adj := p.adjust(p)]




##################
## basic analysis
##################
# all have kmer25 and depth trimming
dat[,min(p, na.rm=TRUE)]
dat[,min(p.adj, na.rm=TRUE)]



#####################
## plots
#####################

# histogram of p-values
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_hist_pvals.png', units='in', res=300)
	hist(dat$p, col='grey', breaks=100)

	dev.off()

	# zoom in on low p-values
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_hist_pvals_zoom.png', units='in', res=300)
	hist(dat$p, col='grey', breaks=8000, xlim=c(0,0.01))

	dev.off()
	
# histogram of FDR-adjusted p-values
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_hist_padjvals.png', units='in', res=300)
	hist(dat$p.adj, col='grey', breaks=100)

	dev.off()

	# zoom in on low p-values
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_hist_padjvals_zoom.png', units='in', res=300)
	hist(dat$p.adj, col='grey', breaks=2000, xlim=c(0,0.1), ylim=c(0,600), main='FDR-corrected p-values under null model')

	dev.off()
	

# -log10(p) vs. genome position
	quartz(width=8, height=6)
	# png(width=8, height=6, file=paste('figures/wfs_nullmodel_p_vs_position', suffix, '.png', sep=''), res=300, units='in')
	plot(dat$POSgen, -log10(dat$p),type='p', cex=0.2, xlab='Genome position', ylab='-log10(p)', col=rgb(0.1,0.1,0.1, 0.2), ylim=c(-0.2,6))
#		abline(h=-log10(0.0001), col='red', lty=2)

		# add 25kmer loci outliers
		dat[p.adj3<0.2,points(POSgen, -log10(p), cex=0.4, col='red')]

		# add LG labels
	col = 'red'
	lgs <- sort(unique(dat[,CHROM]))
	for(j in 1:length(lgs)){
		rng <- range(dat[CHROM==lgs[j], POSgen])
		if(j %% 2 == 0) lines(x=rng, y=c(-0.1,-0.1), col=col, lwd=2)
		if(j %% 2 == 1) lines(x=rng, y=c(-0.05,-0.05), col=col, lwd=2)
		text(x=mean(rng), y=-0.18, labels=lgs[j], col=col, cex=0.5)
	
	}

		# add mean p
#	lines(meanp$mids, meanp$p, col='orange')

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

