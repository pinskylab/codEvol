# Calculate probability of null model producing results as extreme as our observations
# run after wfs_nullmodel_function.r and wfs_nullmodel_combine.r
# Script analyzes the output

##########
## Prep 
##########
require(data.table)

# Choose one set to read in
	# 1907-2014
	locnms <- fread('data_2018.09.05/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	locnms11 <- fread('data_2018.09.05/Frequency_table_Lof07_Lof11.txt', header=TRUE); 
	setnames(locnms11, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # for 1907-2011 (for comparison)
	suffix <- '_07-14'

	# 1907-2011
	locnms <- fread('data_2018.09.05/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	suffix <- '_07-11'

	# 2011-2014
	locnms <- fread('data_2018.09.05/Frequency_table_Lof11_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	suffix <- '_11-14'

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

# calculate a running mean -log10(p-value) (FDR-adjusted)
stp = 1e5
mids <- seq(stp/2, max(dat$POSgen), by=stp)
meanp <- data.frame(mids=mids, p=rep(NA, length(mids)), p.adj=rep(NA, length(mids)), p.adj2=rep(NA, length(mids))) # for mean of p and p-adj
nrow(meanp)
for(j in 1:nrow(meanp)){ # takes a couple minutes
	if(j %% 100 == 0) cat(j)
	inds <- dat$POSgen >= (meanp$mids[j]-stp/2) & dat$POSgen < (meanp$mids[j]+stp/2)
	meanp$p[j] <- mean(-log10(dat$p[inds]))
	meanp$p.adj[j] <- mean(-log10(dat$p.adj[inds]))
	meanp$p.adj2[j] <- mean(-log10(dat$p.adj2[inds]))
}

######################
# basic sanity checks
######################

# number of unique p-values for a single sample size
# Need a different sample size if using Can
if(suffix != '_07-11-14') numps <- dat[cnt1==42 & cnt2==48 & !is.na(p),length(unique(p)), by=.(Freq_1, Freq_2)]
if(suffix == '_07-11-14') numps <- dat[cnt1==42 & cnt2==48 & cnt3==40 & !is.na(p),length(unique(p)), by=.(Freq_1, Freq_2, Freq_3)]
	summary(numps) # should all be 1

# number of unique p-values for all sample sizes and starting and ending frequencies
if(suffix != '_07-11-14') cnts <- unique(dat[,.(cnt1, cnt2)], by=c('cnt1', 'cnt2'))
if(suffix == '_07-11-14') cnts <- unique(dat[,.(cnt1, cnt2, cnt3)], by=c('cnt1', 'cnt2', 'cnt3'))

`if(suffix != '_07-11-14'){
	numps <- dat[!is.na(p) & cnt1==cnts[1,cnt1] & cnt2==cnts[1,cnt2], .(nump=length(unique(p))), by=.(Freq_1, Freq_2)]
	numps[,':=' (cnt1=cnts[1,cnt1], cnt2=cnts[1,cnt2])]
	for(i in 2:nrow(cnts)){
		temp <- dat[!is.na(p) & cnt1==cnts[i,cnt1] & cnt2==cnts[i,cnt2], .(nump=length(unique(p))), by=.(Freq_1, Freq_2)]
		temp[,':=' (cnt1=cnts[i,cnt1], cnt2=cnts[i,cnt2])]
		numps <- rbind(numps, temp)
	}
} 
if(suffix == '_07-11-14'){
	numps <- dat[!is.na(p) & cnt1==cnts[1,cnt1] & cnt2==cnts[1,cnt2] & cnt3==cnts[1,cnt3], .(nump=length(unique(p))), by=.(Freq_1, Freq_2, Freq_3)]
	numps[,':=' (cnt1=cnts[1,cnt1], cnt2=cnts[1,cnt2], cnt3=cnts[1,cnt3])]
	for(i in 2:nrow(cnts)){
		temp <- dat[!is.na(p) & cnt1==cnts[i,cnt1] & cnt2==cnts[i,cnt2] & cnt3==cnts[i,cnt3], .(nump=length(unique(p))), by=.(Freq_1, Freq_2, Freq_3)]
		temp[,':=' (cnt1=cnts[i,cnt1], cnt2=cnts[i,cnt2], cnt3=cnts[i,cnt3])]
		numps <- rbind(numps, temp)
	}
} 
`

unique(numps$nump) # should all be 1

##################
## basic analysis
##################
# all have kmer25 and depth trimming
dat[,min(p, na.rm=TRUE)]
dat[,min(p.adj, na.rm=TRUE)]
#dat[,min(p.adj2, na.rm=TRUE)]
dat[,min(p.adj3, na.rm=TRUE)]


# most diverged loci? using p.adj
#pthresh <- 0.1
#pthresh <- 0.01
#sum(selinds <- dat$p.adj <pthresh & !(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), na.rm=TRUE) # not in Inversions or on Unplaced (the ones we want)
#sum(dat$p.adj <pthresh & dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')) # within inversions
#sum(dat$p.adj <pthresh & dat$CHROM %in% c('Unplaced')) # on Unplaced
#
#	# with kmer25 trimming
#sum(selinds <- dat$p.adj < pthresh & !(dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & dat$kmer25==1, na.rm=TRUE) # not in Inversions or on Unplaced (the ones we want)
#sum(dat$p.adj <pthresh & dat$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12') & dat$kmer25==1) # within inversions
#sum(dat$p.adj <pthresh & dat$CHROM %in% c('Unplaced') & dat$kmer25==1) # on Unplaced

# most diverged loci? using p.adj3 (after inversions and unplaced removed, all with kmer25 and depth trimming)
pthresh <- 0.2
#pthresh <- 0.1
sum(selinds <- dat$p.adj3 < pthresh, na.rm=TRUE) # not in Inversions or on Unplaced (the ones we want)



dat[selinds,]
summary(dat$ABS_DIFF[selinds])
dat[selinds, table(CHROM)]

	# find distance among outlier loci
	cands <- as.data.table(dat[selinds,])
	setkey(cands, CHROM, POS) # sort by position
	cands[,ndist := NA] # nearest neighbor at an earlier position (measure only to left)
	for(i in 1:nrow(cands)){
		j <- which(cands$CHROM == cands$CHROM[i]) # other loci on same chromosome
		j <- j[j<i] # remove focal locus and any later loci
		if(length(j)>0) cands$ndist[i] <- min(abs(cands$POS[i] - cands$POS[j]))
	} 

	selinds2 <- which(cands$ndist < 10000 & !is.na(cands$ndist)) # loci close to a previous locus
	selinds2 <- unique(sort(c(selinds2, selinds2-1))) # add locus before
	length(selinds2) # how many outlier loci that are in clusters

	# label each cluster
	indx <- 1
	cands[,cluster:=as.numeric(NA)]
	for(i in selinds2){
		cands[i,cluster:=indx]
		if(i < nrow(cands)){
			if(cands[i+1,ndist>=100000 | is.na(ndist)]) indx <- indx+1
		}
	}

	cands[!is.na(cluster),length(unique(cluster))] # how many clusters?
	cands[,table(cluster)]

	cands[selinds2,] # the outlier loci nearby to other outlier loci
	cands[selinds2, .(CHROM, POS, cnt1, cnt2, cnt3, Freq_1, Freq_2, Freq_3, p, p.adj, p.adj2, ndist, cluster)][!(Freq_1 %in% c(0,1)),]

	# examine all loci in a region
		# LG08
	dat[CHROM=='LG08' & POS >= 16344554 & POS <=16347801,]
	locnms11[CHROM=='LG08' & POS >= 16344554 & POS <=16347801,] # same region in 1907-2011

		# LG06
	dat[CHROM=='LG06' & POSgen >= 143400000 & POSgen <=143600000,]
	locnms11[CHROM=='LG06' & POS >= 2780283 & POS <=2866691,] # same region in 1907-2011

# write out
if(suffix != '_07-11-14') write.csv(cands[selinds2,.(CHROM, POS, POSgen, locusnum, cluster, ndist, cnt1, cnt2, Freq_1, Freq_2, ABS_DIFF, p, n, p.adj, p.adj3)], file=paste('analysis/wfs_nullmodel_candidates', suffix, '.csv', sep=''))
if(suffix == '_07-11-14') write.csv(cands[selinds2,.(CHROM, POS, POSgen, locusnum, cluster, ndist, cnt1, cnt2, cnt3, Freq_1, Freq_2, Freq_3, ABS_DIFF, p, n, p.adj, p.adj3)], file=paste('analysis/wfs_nullmodel_candidates', suffix, '.csv', sep=''))

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
	
# p-values vs. abs_diff
	require(RColorBrewer)
	colrmp <- colorRamp(brewer.pal(9, 'RdYlBu'))
	colrmpfun <- function(x){
		rgbs <- colrmp(x)
		return(rgb(rgbs, maxColorValue=256))
	}
	dat[,cntmn:=2/((1/cnt1)+(1/cnt2))] # harmonic mean sample size
	dat[,cntmin:=pmin(cnt1,cnt2)] # min sample size

	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_pvals_vs_abs_diff.png', units='in', res=300)
#	plot(dat$ABS_DIFF, -log10(dat$pmax),type='p', cex=0.2, xlab='Frequency change', ylab='-log10(p)', col=rgb(0.1,0.1,0.1, 0.2)) # grey
	plot(dat$ABS_DIFF, -log10(dat$p),type='p', cex=0.2, xlab='Frequency change', ylab='-log10(p)', col=colrmpfun(dat[,cntmn/49])) # colored by harmonic mean sample size
#	plot(dat$ABS_DIFF, -log10(dat$pmax),type='p', cex=0.2, xlab='Frequency change', ylab='-log10(p)', col=colrmpfun(dat[,cntmin/48])) # colored by min sample size
	plot(dat$Freq_2 - dat$Freq_1, -log10(dat$p),type='p', cex=0.2, xlab='Frequency change', ylab='-log10(p)', col=colrmpfun(dat[,Freq_1])) # colored Freq_1
	
	dev.off()

	# just a single sample size
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_pvals_vs_abs_diff_cnt46,44.png', units='in', res=300)
	dat[cnt1==46 & cnt2==44, plot(ABS_DIFF, -log10(p), type='p', cex=0.2, xlab='Frequency change', ylab='-log10(p)', col=rgb(0.1,0.1,0.1,0.1), main='null model analysis, ctn1=46, cnt2=44')] # not colored
	dat[cnt1==42 & cnt2==48, plot(Freq_2-Freq_1, -log10(p), type='p', cex=0.2, xlab='Frequency change', ylab='-log10(pmax)', col=colrmpfun(Freq_1), main='null model analysis, ctn1=46, cnt2=44')] # colored by Freq_1
	dev.off()

	# just a single sample size and starting frequency
	dat[cnt1==42 & cnt2==48 & abs((Freq_2-Freq_1)-(7/46-9/44))<0.01, plot(Freq_2-Freq_1, -log10(p), type='p', cex=0.2, xlab='Frequency change', ylab='-log10(p)', col='red', main='null model analysis, ctn1=46, cnt2=44', ylim=c(0,3))]
	dat[cnt1==42 & cnt2==48 & abs(Freq_1-6/46) < 0.01 & abs(Freq_2-2/44)<40/44, plot(Freq_2-Freq_1, -log10(pmax), type='p', cex=0.2, xlab='ABS_DIFF', ylab='-log10(p)', col='blue', main='null model analysis, ctn1=46, cnt2=44')]
	dat[cnt1==46 & cnt2==44 & abs(Freq_1-0.3) < 0.01, points(ABS_DIFF, -log10(pmax), type='p', cex=0.2, xlab='Frequency change', ylab='-log10(p)', col='green', main='null model analysis, ctn1=46, cnt2=44')]

# FDR p-values vs. abs_diff
	# png(width=5, height=4, filename='analysis/figures/wfs_nullmodel_padjvals_vs_abs_diff.png', units='in', res=300)
#	plot(dat$ABS_DIFF, -log10(dat$p.adj),type='p', cex=0.2, xlab='Frequency change', ylab='-log10(FDR-adjusted p)', col=rgb(0.1,0.1,0.1, 0.2))
	plot(dat$ABS_DIFF, -log10(dat$p.adj2),type='p', cex=0.2, xlab='Frequency change', ylab='-log10(FDR-adjusted p)', col=colrmpfun(dat[,cntmn/49]))

	dev.off()
	
# p-value vs. Freq_1

	# just a single sample size
	quartz(width=6,height=6)
	# png(width=5, height=4, filename=paste('analysis/figures/wfs_nullmodel_pmax_vs_Freq_1_w_Freq_2_cnt46,44', suffix, '.png', sep=''), units='in', res=300)
	dat[cnt1==46 & cnt2==44, plot(Freq_1, -log10(pmax), type='p', cex=0.4, xlab='Initial frequency', ylab='-log10(pmax)', col=colrmpfun(Freq_2), main='null model analysis, ctn1=46, cnt2=44')]
		# add a legend
		points(rep(0,10), seq(3,4,length.out=10), col=colrmpfun(seq(0,1,length.out=10)), pch=16) 
		text(rep(0.05,5), seq(3,4,length.out=5), labels=seq(0,1,length.out=5), cex=0.5)
		text(0.025, 4.2, 'Freq_2', cex=0.8)


# Freq_1 vs. Freq_2 with p
	#cnt1=46;cnt2=44 # for Lof
	ct1=32;ct2=40 # for Can
	inds = dat$kmer25==1; suff2='_25kmer' # for trimming to 25-kmer
#	inds = TRUE; suff2='' # for not trimming
	
	require(RColorBrewer)
	colrmp <- colorRamp(brewer.pal(9, 'OrRd'))
	colrmpfun <- function(x){
		rgbs <- colrmp(x)
		return(rgb(rgbs, maxColorValue=256))
	}

	# just a single sample size
	quartz(width=6,height=6)
	# png(width=5, height=5, filename=paste('figures/wfs_nullmodel_Freq_1_vs_Freq_2_w_pmax_cnt', ct1, ',', ct2, suffix, suff2, '.png', sep=''), units='in', res=300)
	dat[inds & cnt1==ct1 & cnt2==ct2, plot(Freq_1, Freq_2, type='p', cex=0.8, pch=16, xlab='Freq_1', ylab='Freq_2', col=colrmpfun(-log10(p)/6), main=paste('null model analysis, cnt1=', ct1, 'cnt2=', ct2, '\n', suffix, suff2), xlim=c(0,1), ylim=c(0,1))]
		# add a legend
		points(rep(0,10), seq(0.6,0.8,length.out=10), col=colrmpfun(seq(0,1,length.out=10)), pch=16) 
		text(rep(0.05,5), seq(0.6,0.8,length.out=5), labels=seq(0,6,length.out=5), cex=0.5)
		text(-0.05, 0.85, '-log10(pmax)', cex=0.8, pos=4)
		
	dev.off()


# Freq_1 vs. Freq_2 with p.adj
	#cnt1=46;cnt2=44 # for Lof
	ct1=32;ct2=40 # for Can

	require(RColorBrewer)
	colrmp <- colorRamp(brewer.pal(9, 'OrRd'))
	colrmpfun <- function(x){
		rgbs <- colrmp(x)
		return(rgb(rgbs, maxColorValue=256))
	}

	# just a single sample size
	quartz(width=6,height=6)
	# png(width=5, height=4, filename=paste('analysis/figures/wfs_nullmodel_Freq_1_vs_Freq_2_w_p.adj_cnt46,44', suffix, '.png', sep=''), units='in', res=300)
	dat[cnt1==ct1 & cnt2==ct2, plot(Freq_1, Freq_2, type='p', cex=0.8, pch=16, xlab='Freq_1', ylab='Freq_2', col=colrmpfun(-log10(p.adj)/3), main=paste('null model analysis, cnt1=', ct1, 'cnt2=', ct2, '\n', suffix))]
		# add a legend
		points(rep(0,10), seq(0.6,0.8,length.out=10), col=colrmpfun(seq(0,1,length.out=10)), pch=16) 
		text(rep(0.05,5), seq(0.6,0.8,length.out=5), labels=seq(0,2,length.out=5), cex=0.5)
		text(-0.05, 0.85, '-log10(p.adj)', cex=0.8, pos=4)
		
	dev.off()

	
# -log10(p) vs. genome position
	quartz(width=8, height=6)
	# png(width=8, height=6, file=paste('figures/wfs_nullmodel_p_vs_position', suffix, '.png', sep=''), res=300, units='in')
	plot(dat$POSgen, -log10(dat$p),type='p', cex=0.2, xlab='Genome position', ylab='-log10(p)', col=rgb(0.1,0.1,0.1, 0.2), ylim=c(-0.2,6))
#		abline(h=-log10(0.0001), col='red', lty=2)

		# add 25kmer loci outliers
		dat[kmer25==1 & p.adj3<0.2,points(POSgen, -log10(p), cex=0.4, col='red')]

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


# -log10(p) vs. genome position (zoom)
	quartz(width=8, height=6)
	xlims=c(1.432e8, 1.439e8)
	# png(width=8, height=6, file=paste('figures/wfs_nullmodel_p_vs_position_zoom', suffix, '.png', sep=''), res=300, units='in')
	plot(dat$POSgen, -log10(dat$p),type='p', cex=0.2, xlab='Genome position', ylab='-log10(p)', col=rgb(0.1,0.1,0.1, 0.2), xlim=xlims)
#		abline(h=-log10(0.0001), col='red', lty=2)

		# add LG labels
	lgs <- sort(unique(dat[,CHROM]))
	for(j in 1:length(lgs)){
		rng <- range(dat[CHROM==lgs[j], POSgen])
		if(j %% 2 == 0) lines(x=rng, y=c(0,0), col=col, lwd=2)
		if(j %% 2 == 1) lines(x=rng, y=c(0.05,0.05), col=col, lwd=2)
		text(x=mean(rng), y=0.1, labels=lgs[j], col=col, cex=0.5)
	
	}

		# add mean p
	lines(meanp$mids, meanp$p, col='orange')

	dev.off()
	
# -log10(p.adj) vs. genome position
	col='red'
	quartz(width=8, height=6)
	# png(width=8, height=6, file=paste('figures/wfs_nullmodel_padj_vs_position', suffix, '.png', sep=''), res=300, units='in')
	plot(dat$POSgen, -log10(dat$p.adj),type='p', cex=0.2, xlab='Genome position', ylab='-log10(FDR-adjusted p)', col=rgb(0.1,0.1,0.1, 0.2))
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



# -log10(p.adj) vs. genome position for one chromosome (zoom)
	col='red'
	quartz(width=8, height=6)
	xlims <- c(00000, 22000000); chrom='LG19' # LG19
	# png(width=8, height=6, file=paste('figures/wfs_nullmodel_padj_vs_position_zoomLG19', suffix, '.png', sep=''), res=300, units='in')
	dat[CHROM==chrom, plot(POS,-log10(p.adj),type='p', cex=0.2, xlab='Genome position', ylab='-log10(FDR-adjusted p)', col=rgb(0.1,0.1,0.1, 0.8), xlim=xlims, main=paste('Zoom in on outliers on', chrom))]
		abline(h=-log10(0.05), col='red', lty=2)

	# harder to add meanp since based on POS rather than POSgen

	dev.off()

