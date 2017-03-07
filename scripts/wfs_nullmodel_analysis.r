# Calculate probability of null model producing results as extreme as our observations
# run after wfs_nullmodel_function.r and wfs_nullmodel_combine.r
# Script analyzes the output

##########
## Prep 
##########
require(data.table)
locnms <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals.rdata') # dat

# make a nucleotide position for the whole genome
chrmax <- locnms[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])

setkey(locnms, CHROM)
setkey(chrmax, CHROM)
locnms <- locnms[chrmax[,.(CHROM, start)], ]
locnms[,POSgen:=POS+start]

# trim locus names to match rest of data
locnms <- locnms[!(locnms$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')),] # trim out inversions
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting

	# make sure it matches
	nrow(dat)
	nrow(locnms)


# merge
dat <- merge(dat, locnms[,.(locusnum, CHROM, POS, POSgen, Freq_1, Freq_2, ABS_DIFF)], by='locusnum')

# calculate a running mean -log10(p-value)
stp = 1e5
meanp <- data.frame(mids <- seq(stp/2, max(dat$POSgen), by=stp), mean <- rep(NA, length(mids))) # for mean of p-adj
nrow(meanp)
for(j in 1:nrow(meanp)){
	if(j %% 100 == 0) cat(j)
	inds <- dat$POSgen >= (meanp$mids[j]-stp/2) & dat$POSgen < (meanp$mids[j]+stp/2)
	meanp$mean[j] <- mean(-log10(dat$p.adj[inds]))
}


##################
## basic analysis
##################
# most diverged loci?
sum(selinds <- dat$p.adj <0.025 & dat$CHROM != 'Unplaced')

dat[selinds,]
summary(dat$ABS_DIFF[selinds])

	# find distance among loci
	cands <- as.data.table(dat[selinds,])
	setkey(cands, CHROM, POS) # sort by position
	cands[,ndist := NA] # nearest neighbor at an earlier position (measure only to left)
	for(i in 1:nrow(cands)){
		j <- which(cands$CHROM == cands$CHROM[i]) # other loci on same chromosome
		j <- j[j<i] # remove focal locus and any later loci
		if(length(j)>0) cands$ndist[i] <- min(abs(cands$POS[i] - cands$POS[j]))
	} 

	selinds2 <- which(cands$ndist < 100000 & !is.na(cands$ndist)) # loci close to a previous locus
	selinds2 <- sort(c(selinds2, selinds2-1)) # add locus before

	cands[selinds2,]


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
	hist(dat$p.adj, col='grey', breaks=2000, xlim=c(0,0.05), ylim=c(0,80), main='FDR-corrected p-values under null model')

	dev.off()
	
# -log10(p) vs. genome position
	col='red'
	quartz(width=8, height=6)
	#png(width=8, height=6, file='analysis/figures/wfs_nullmodel_p_vs_position.png', res=300, units='in')
	plot(dat$POSgen, -log10(dat$p.adj),type='p', cex=0.2, xlab='Genome position', ylab='-log10(FDR-adjusted p)', col=rgb(0.1,0.1,0.1, 0.2))
		abline(h=-log10(0.05), col='red', lty=2)
		abline(h=-log10(0.03), col='red', lty=3)

		# add LG labels
	lgs <- sort(unique(locnms[i,CHROM]))
	for(j in 1:length(lgs)){
		rng <- range(locnms[CHROM==lgs[j], POSgen])
		if(j %% 2 == 0) lines(x=rng, y=c(0,0), col=col, lwd=2)
		if(j %% 2 == 1) lines(x=rng, y=c(0.015,0.015), col=col, lwd=2)
		text(x=mean(rng), y=0.03, labels=lgs[j], col=col, cex=0.5)
	
	}

		# add mean p
	lines(meanp$mids, meanp$mean, col='orange')

	dev.off()

