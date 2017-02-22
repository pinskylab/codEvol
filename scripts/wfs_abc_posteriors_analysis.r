# run after wfs_abc.r (possibly sbatch version) and wfs_abc_hpds_calc.r

################################
# load functions and prep data
################################
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(RColorBrewer)
	require(data.table)
	require(boa)
	require(hexbin)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(RColorBrewer, lib.loc="/projects/cees/lib/R_packages/")
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

# load data
dat <- fread('analysis/LOF_07_to_LOF_S_14.wfabc', skip=2) # the data on samples sizes and observed allele frequencies, from WFABC input file
locnms <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_abc_hpds.rdata') # the HPD calculations


# prep data
	# extract sample frequencies
sampsize <- dat[seq(1,nrow(dat),by=2),] # sample sizes in # of chromosomes
obsfrqs <- dat[seq(2,nrow(dat),by=2),] # allele counts in # chromosomes
rm(dat)
setnames(obsfrqs, 1:2, c('count1', 'count2'))
setnames(sampsize, 1:2, c('size1', 'size2'))
obsfrqs <- cbind(obsfrqs, sampsize)
obsfrqs[,f1:=count1/size1]
obsfrqs[,f2:=count2/size2]
obsfrqs[,diff:=f2-f1]
obsfrqs[,locusnum:=1:nrow(obsfrqs)] # add locusnumber

	# trim locus names to match rest of data
#locnms <- locnms[CHROM=='LG03',] # if only looking at LG03
locnms <- locnms[!(locnms$CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12')),] # trim out inversions
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting

# make a nucleotide position for the whole genome
chrmax <- locnms[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])

setkey(locnms, CHROM)
setkey(chrmax, CHROM)
locnms <- locnms[chrmax[,.(CHROM, start)], ]
locnms[,POSgen:=POS+start]


#########################################################
## Initial exploration
#########################################################

# typical sample sizes
obsfrqs[,summary(size1)]
obsfrqs[,summary(size2)]


# find loci potentially under selection
# don't differentiate pos vs. neg: just differentiates which allele
selinds <- (hpds$s$p < 0.0001 & hpds$s$mean > 0 & !is.na(hpds$s$mean)) | (hpds$s$p < 0.0001 & hpds$s$mean < 0 & !is.na(hpds$s$mean))

sum(selinds) # number of loci

selinds2 <- locnms$CHROM!='Unplaced' & hpds$s$p < 0.001 & abs(hpds$s$u95) >0.5 & abs(hpds$s$l95) >0.5 # strongest evidence
	sum(selinds2)

# Examine candidates
print(locnms[selinds,], nrow=50)
print(obsfrqs[selinds,], nrow=50)

hpds$s[posinds,]

print(locnms[selinds2,], nrow=sum(selinds2))
print(obsfrqs[selinds,], nrow=50)

table(locnms[selinds2,CHROM])


# find distance among loci
cands <- locnms[selinds2,]
setkey(cands, CHROM, POS) # sort by position
cands[,ndist := NA] # nearest neighbor at an earlier position (measure only to left)
for(i in 1:nrow(cands)){
	j <- which(cands$CHROM == cands$CHROM[i]) # other loci on same chromosome
	j <- j[j<i] # remove focal locus and any later loci
	if(length(j)>0) cands$ndist[i] <- min(abs(cands$POS[i] - cands$POS[j]))
} 

selinds3 <- which(cands$ndist < 100000 & !is.na(cands$ndist)) # loci close to a previous locus
selinds3 <- sort(c(selinds3, selinds3-1)) # add locus before

cands[selinds3,]

################################
# plots of the results
################################
load('analysis/wfs_abc_hpds.rdata')

# plotting functions
colrmp <- colorRamp(colors=brewer.pal(11, 'RdYlBu'))
getcol <- function(x, alpha=1){ rgbs <-colrmp(x); return(rgb(rgbs[,1], rgbs[,2], rgbs[,3], alpha*256, maxColorValue=256))}

# histogram of s
hist(abs(hpds$s$mean), breaks=seq(0,1,by=0.05), col='grey')
hist(abs(hpds$s$mean), breaks=seq(0,1,by=0.05), col='grey', xlim=c(0.5,1), ylim=c(0,10000))


# Initial frequency vs. posterior s
	plot(hpds$f1$mean, hpds$s$mean, xlab='Posterior mean initial frequency', ylab='Posterior mean s', col=getcol(hpds$s$mean/2+0.5))
	points(hpds$f1$mean[selinds], hpds$s$mean[selinds], pch=16)


# Obs init freq vs. posterior true init freq
	png('analysis/figures/wfs_abc_postf1_vs_obsf1_with_s.png')
	plot(obsfrqs$f1, hpds$f1$mean, col=getcol(hpds$s$mean/2+0.5, 0.2), xlab='Observed initial frequency', ylab='Posterior mean initial frequency')
	points(obsfrqs$f1[selinds], hpds$f1$mean[selinds], pch=16)

	dev.off()

# Obs final freq vs. posterior true final freq
	png('analysis/figures/wfs_abc_postf2_vs_obsf2_with_s.png')
	plot(obsfrqs$f2, hpds$f2$mean, col=getcol(hpds$s$mean/2+0.5, 0.2), xlab='Observed final frequency', ylab='Posterior mean final frequency')
	points(obsfrqs$f2[selinds], hpds$f2$mean[selinds], pch=16)

	dev.off()

# Obs init freq vs. posterior init freq sampled
	png('analysis/figures/wfs_abc_postf1samp_vs_obsf1_with_s.png')
	plot(obsfrqs$f1, hpds$f1samp$mean, col=getcol(hpds$s$mean/2+0.5, 0.2), xlab='Observed initial frequency', ylab='Posterior mean observed initial frequency')
	points(obsfrqs$f1[selinds], hpds$f1samp$mean[selinds], pch=16)

	dev.off()

# Obs final freq vs. posterior final freq sampled
	png('analysis/figures/wfs_abc_postf2samp_vs_obsf2_with_s.png')
	plot(obsfrqs$f2, hpds$f2samp$mean, col=getcol(hpds$s$mean/2+0.5, 0.2), xlab='Observed final frequency', ylab='Posterior mean observed final frequency')
	points(obsfrqs$f2[selinds], hpds$f2samp$mean[selinds], pch=16)

	dev.off()


# Initial frequency vs. change in frequency vs. posterior s
# overplot loci "under selection"
	quartz(width=7, height=6)
	#png(width=7, height=6, filename='analysis/figures/wfs_abc_obsdiff_vs_obsf1_with_s.png', units='in', res=300)
	layout(matrix(c(1,2),ncol=2),widths=c(4,1))
	par(mai=c(1,1,0.2, 0.2), las=1, mgp=c(2.4, 0.6, 0), tcl=-0.2)
	plot(obsfrqs$f1, obsfrqs$diff, xlab='Observed initial frequency', ylab='Observed change in frequency', col=getcol(hpds$s$mean/2+0.5, 0.1), ylim=c(-1,1)) # background
	points(obsfrqs$f1[selinds], obsfrqs$diff[selinds], col='black', pch=16, cex=0.8)
	points(obsfrqs$f1[selinds2], obsfrqs$diff[selinds2], col='green', pch=16, cex=0.4)
	plot(rep(1,100),seq(-1,1,length.out=100), col=getcol(seq(-1,1,length.out=100)/2+0.5), bty='n', xaxt='n', ylab='Posterior mean s', xlab='', pch=15, cex=0.7) # legend

	dev.off()


# Initial frequency vs. change in frequency vs. AVERAGE posterior s (as a color)
# overplot loci "under selection"
# calculate average s within a box
	# prep data
	aves <- unique(obsfrqs[,.(f1=floor(f1/0.01)*0.01+0.005, diff=floor(diff/0.01)*0.01+0.005)])
	aves[,s:=NA]
	for(i in 1:nrow(aves)){
	#	if(i %% 100 == 0) cat(i)
		inds <- abs(obsfrqs$f1 - aves$f1[i]) < 0.005 & abs(obsfrqs$diff - aves$diff[i]) < 0.005
		aves$s[i] <- mean(hpds$s$mean[inds])
	}
	aves <- aves[!is.na(s),] # not sure why some are NA

	quartz(width=7, height=6)
	#png(width=7, height=6, filename='analysis/figures/wfs_abc_obsdiff_vs_obsf1_with_s.png', units='in', res=300)
	layout(matrix(c(1,2),ncol=2),widths=c(4,1))
	par(mai=c(1,1,0.2, 0.2), las=1, mgp=c(2.4, 0.6, 0), tcl=-0.2)
	plot(aves$f1, aves$diff, xlab='Observed initial frequency', ylab='Observed change in frequency', col=getcol(aves$s/2+0.5, 0.5), ylim=c(-1,1), pch=16, cex=0.4) # background
	points(obsfrqs$f1[selinds], obsfrqs$diff[selinds], col='black', pch=16, cex=0.8)
	points(obsfrqs$f1[selinds2], obsfrqs$diff[selinds2], col='green', pch=16, cex=0.4)
	plot(rep(1,100),seq(-1,1,length.out=100), col=getcol(seq(-1,1,length.out=100)/2+0.5), bty='n', xaxt='n', ylab='Posterior mean s', xlab='', pch=15, cex=0.7) # legend

	dev.off()

# Plot posterior mean s vs. allele frequency change, by initial frequency
	quartz(width=7, height=6)
	# png(width=7, height=7, filename='analysis/figures/wfs_abc_s_vs_obsdiff_by_obsf1.png', units='in', res=300)
	f1s <- seq(0,1.1,by=0.1)
	par(mfrow=c(3,4), las=1, tcl=-0.3, mgp=c(2.4,0.6,0), mai=c(0.5, 0.5, 0.3, 0.05))
	for(i in 2:length(f1s)){
		inds <- (obsfrqs$f1 >= f1s[i-1]) & (obsfrqs$f1 < f1s[i])
		plot(obsfrqs$diff[inds], hpds$s$l95[inds], xlab='Obs freq change', ylab='Posterior mean s', cex=0.05, main=paste('initial ', f1s[i-1], '-', f1s[i], sep=''), xlim=c(-0.52, 0.52), ylim=c(-1, 1), col=rgb(0.8,0.8,0.8,0.1))
		points(obsfrqs$diff[inds], hpds$s$u95[inds], col=rgb(0.8,0.8,0.8,0.1), cex=0.05)
		points(obsfrqs$diff[inds], hpds$s$mean[inds], col=rgb(0,0,0,0.1), cex=0.05)
	}

	dev.off()


# Plot of s vs. genome position
col='red'
quartz(width=8, height=6)
#png(width=8, height=6, file='analysis/figures/wfs_abc_s_vs_position.png', res=300, units='in')
i <- locnms$CHROM!='Unplaced'
plot(locnms$POSgen[i], abs(hpds$s$mean[i]),type='p', cex=0.2, ylim=c(0,1), xlab='Position', ylab='Posterior mean abs(s)', col=rgb(0.1,0.1,0.1, 0.2))
i2 <- locnms$CHROM!='Unplaced' & hpds$s$p < 0.001 & abs(hpds$s$u95) >0.5 & abs(hpds$s$l95) >0.5
	sum(i2)
points(locnms$POSgen[i2], abs(hpds$s$mean[i2]),cex=0.5, col=col)
for(j in which(i2)) lines(x=rep(locnms$POSgen[j],2), y=abs(c(hpds$s$l95[j], hpds$s$u95[j])), col=col) # too messy for this plot to show the CIs

	# add LG labels
lgs <- sort(unique(locnms[i,CHROM]))
for(j in 1:length(lgs)){
	rng <- range(locnms[CHROM==lgs[j], POSgen])
	if(j %% 2 == 0) lines(x=rng, y=c(0,0), col=col, lwd=2)
	if(j %% 2 == 1) lines(x=rng, y=c(0.015,0.015), col=col, lwd=2)
	text(x=mean(rng), y=0.03, labels=lgs[j], col=col, cex=0.5)
	
}

	# add mean s
mids <- seq(5e5, locnms[CHROM!='Unplaced', max(POSgen)], by=1e6)
means <- rep(NA, length(mids))
for(j in 1:length(mids)){
	inds <- locnms$POSgen > (mids[j]-5e5) & locnms$POSgen < (mids[j]+5e5)
	means[j] <- mean(hpds$s$mean[inds])
}
lines(mids, means, col=rgb(0,0,1, 0.5))

dev.off()



##############################
## explore a single locus
##############################

inds <- hpds$s$p.adj < 0.05 & abs(obsfrqs$diff) < 0.05

# do observed and posterior allele freq match?
plot(obsfrqs$f1[inds], hpds$f1$mean[inds], col=getcol(hpds$s$mean[inds]/2+0.5, 0.2), xlab='Observed initial frequency', ylab='Posterior mean initial frequency')

# pick one
head(which(inds))
i <- 36

obsfrqs[i,]
locnms[i,]

# look at target sumstats and posterior sumstats for this locus
obsstats <- fread('analysis/LOF_07_to_LOF_S_14.w_obs_stats.txt') # the observations of Fsi and Fsd
load('analysis/wfs_targ.rdata')
targ[i,]
obsstats[i,]
hpds$f1samp[i,] # matches the normalized f1samp.n in targ
hpds$f1[i,] # close to observed init freq
hpds$fsd[i,] # matches normalized fsd.n in targ
hpds$fsi[i,] # matches normalized fsi.n in targ
hpds$s[i,]

	# calculate fs' by hand for this locus
t <- 11 # generations
x <- obsfrqs[i,f1]; if(x>0.5) x <- 1-x # maf 1
y <- obsfrqs[i,f2]; if(y>0.5) y <- 1-y # maf 2
z <- (x+y)/2
fs <- (x-y)^2/(z*(1-z))
nx <- obsfrqs[i,size1]
ny <- obsfrqs[i,size2]
nstar <- 2/(1/nx + 1/ny)
fsprime <- 1/t * (fs*(1-1/(2*nstar)) - 2/nstar)/((1+fs/4)*(1-1/ny))
fsprime
obsstats[i,] # they match: great

# look at posterior samples
locset <- '000013-002257' # could find this from a script, but fast for just one locus to do this by hand
infile <- paste('analysis/temp/wfs_abc_sampsize', locnms[i,N_CHR_1], ',', locnms[i,N_CHR_2], '_locus', locset, '.csv.gz', sep='')

results <- read.csv(file=gzfile(infile))
dim(results)
results <- results[results$locus==i,] # trim to focal locus
dim(results)

hist(results$s)



# Examine MCMC sample distributions for individual loci
# WOULD NEED TO REWRITE THIS TO READ IN THE APPROPRIATE FILE
i <- which(obsfrqs$f1 > 0.8 & obsfrqs$f1 < 0.90 & obsfrqs$diff < -0.25) # initially high frequency loci under negative selection
j <- 2 # pick which one
plot(density(posts$s[i[j],]))
plot(density(posts$f1[i[j],])); abline

hexbinplot(posts$s[i[j],] ~ posts$f1[i[j],], colramp=rf, main=paste('locus', i[j]), xlab='posterior mean initial frequency', ylab='posterior mean s')

hexbinplot(posts$s[i[j],] ~ posts$ne[i[j],], colramp=rf, main=paste('locus', i[j]), xlab='posterior mean Ne', ylab='posterior mean s')
