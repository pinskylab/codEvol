# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
	require(plyr)
	require(parallel)
	ncores=3
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	ncores=20
}

# read in data (choose one)
dat <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=56; c2=48; nm='1907-2014'; gen=11 # for 1907 vs. 2014. sample sizes

# make a nucleotide position for the whole genome
chrmax <- dat[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])

setkey(dat, CHROM)
setkey(chrmax, CHROM)
dat <- dat[chrmax[,.(CHROM, start)], ]
dat[,POSgen:=POS+start]

# running mean of allele frequency difference
options(warn=1)
f <- function(i) {
	is.near <- abs(as.numeric(dat$POSgen[i] - dat$POSgen)) < 5e5 & dat$CHROM[i] == dat$CHROM
	mean(dat$ABS_DIFF[is.near])
}
f2 <- function(i) {
	is.near <- abs(as.numeric(dat$POSgen[i] - dat$POSgen)) < 5e4 & dat$CHROM[i] == dat$CHROM
	mean(dat$ABS_DIFF[is.near])
}
dat$absdiff_mean <- simplify2array(mclapply(1:nrow(dat), f, mc.cores=ncores)) # use multicores to go faster. about 10 min on 20 cores. a long time on 3 cores... didn't finish on 3 cores
dat$absdiff_mean1e5 <- simplify2array(mclapply(1:nrow(dat), f2, mc.cores=ncores))

save(dat, file='analysis/Frequency_table_Lof07_Lof14_runmean.rdata')

# load back in
load('analysis/Frequency_table_Lof07_Lof14_runmean.rdata')


# plot frequency difference
#dat[,plot(POSgen, ABS_DIFF, pch=16, cex=0.3)]
quartz(height=4, width=8)
#png(height=4, width=8, units='in', res=300, file='analysis/figures/abs_diff_vs_pos.png')
dat[,plot(POSgen/1e6, absdiff_mean, type='l', lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,0.1))] # for 1e6 window
#dat[,plot(POSgen/1e6, absdiff_mean1e5, type='l', lwd=0.3, xlab='Position (Mb)', ylab='Allele frequency change', ylim=c(0,0.2))] # for 1e5 window

col='red'
lgs <- sort(unique(dat[,CHROM]))
for(j in 1:length(lgs)){
	rng <- range(dat[CHROM==lgs[j], POSgen])
	if(j %% 2 == 0) lines(x=rng/1e6, y=c(0,0), col=col, lwd=2, lend=2)
	if(j %% 2 == 1) lines(x=rng/1e6, y=c(0.005,0.005), col=col, lwd=2, lend=2)
	text(x=mean(rng/1e6), y=0.01, labels=lgs[j], col=col, cex=0.5)
	
}

dev.off()


