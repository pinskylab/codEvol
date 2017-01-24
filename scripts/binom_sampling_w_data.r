

# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
	require(plyr)
	require(parallel)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
}

# set up parameters
limstep <- 0.05 # the size for the Freq_1 bins
alimstep <- 0.05 # the size for the ABS_DIFF bins

# read in data
dat <- fread('analysis/Frequency_table'); c1=56; c2=48; nm='1907-2011'; gen=20 # for 1907 vs. 2011. sample sizes
dat <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=56; c2=48; nm='1907-2014'; gen=20 # for 1907 vs. 2014. sample sizes
dat <- fread('analysis/Frequency_table_Lof11_Lof14.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); c1=48; c2=48; nm='2011-2014'; gen=1 # for 2011 vs. 2014

# analyze the data by frequency bin
dat[,f1r:=floor(Freq_1/limstep)*limstep+limstep/2] # round to nearest bin (label with bin center)
dat[f1r==(1+limstep/2),f1r:=(1-limstep/2)] # correct the highest bin so <1
dat[f1r>0.5, f1r:=(1-f1r)] # focus on MAF
dat[,absr:=floor(ABS_DIFF/alimstep)*alimstep+alimstep/2]
dat[,f1r:=floor(f1r/limstep)*limstep+limstep/2] # correct numeric error in f1r
tab <- dat[,table(absr,f1r)] # the data to compare simulations against: table of alleles by starting frequency and absolute difference (binned)

# function for use in simulations
	# simulate the sampling process of alleles from a population of known allele freqs
	# i: dummy variable so this can be used with sapply and parSapply
	# nloci: number of loci to simulate
	# f1: initial frequency
	# c1: number of chromosomes in the first sample (e.g., 2x the number of diploid individuals)
	# c2: second sample size
	# alimstep: bin size for splitting up allele frequency changes
	# gen: number of generations (if doing a Wright-Fisher drift simulation)
	# ne: effective population size, measured as # chromosomes (set to Inf for no drift)
getcounts <- function(i, nloci, f1, c1=56, c2=48, alimstep, gen, ne=Inf){ # i is a dummy argument so that it can be used with sapply and parSapply
	if(is.infinite(ne)) f2 <- f1 # for an infinite population size
	if(!is.infinite(ne)) f2 <- wf(ne, f1, gen) # Wright-Fisher sampling process

	simf1 <- replicate(nloci, mean(rbinom(1,c1,f1))/c1) # first sample
	simf2 <- replicate(nloci, mean(rbinom(1,c2,f2))/c2) # second sample

	cnts <- hist(abs(simf1-simf2), breaks=seq(0,1,by=alimstep), plot=FALSE, right=FALSE)$counts # note: not quite right to put all these in the 0.05 starting allele freq bin, since some of the 1st observations are > 0.1 (though not many) due to sampling error
	return(cnts)
}

	# wright fisher forward simulation
	# ne: effective population size in # chromosomes
	# f1: initial allele frequency
	# gen: number of generations
wf <- function(ne,f1,gen){
	currf <- f1 # current allele frequency
	for(i in 1:gen){
		currf <- mean(rbinom(1,ne,currf)/ne)
	}
	return(currf)
}


# simulations (good to run this on a cluster, take 30 min or so with 20 cores on cod)
nsims <- 1000
simsum <- array(NA,dim=c(1/limstep, 1/alimstep, nsims), dimnames=list(f1r=seq(limstep/2,1-limstep/2,by=limstep), absr=seq(alimstep/2,1-alimstep/2,by=alimstep), sim=1:nsims))
f1s <- seq(limstep/2,0.5,by=limstep) # starting allele frequencies
ne=500


if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	cl <- makeCluster(detectCores()-1) # set up cluster on my mac (or another computer), using all but one core
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	cl <- makeCluster(30) # set up 30-core cluster on a cod node
}
clusterExport(cl, c('getcounts', 'wf'))

length(f1s) # how many to do
for(j in 1:length(f1s)){
	print(j)
	nloci <- sum(tab[,j]) # total loci in this starting frequency bin
	#simsum[j,,] <- replicate(nsims, getcounts(1, nloci, f1s[j], alimstep, gen, ne)) # slow, since creating 2*nloci*1000 values and summarizing to a 20x1000 matrix
	
	simsum[j,,] <- parSapply(cl, X=1:nsims, FUN=getcounts, nloci=nloci, f1=f1s[j], c1=c1, c2=c2, alimstep=alimstep, gen=gen, ne=ne) # parallel version: slow with 3 cores, fast with 20
}

stopCluster(cl)

save(simsum, file=paste('analysis/binom_sampling_simsum_', nm, '_ne', ne, '.rdata', sep=''))
print(paste('analysis/binom_sampling_simsum_', nm, '_ne', ne, '.rdata', sep=''))

# analysis
	# for 1907-2011
load('analysis/binom_sampling_simsum_1907-2011_ne500.rdata'); ne=500

	# for 1907-2014
load('analysis/binom_sampling_simsum_1907-2014_ne500.rdata'); ne=500

	# for 2011-2014
load('analysis/binom_sampling_simsum_2011-2014_neInf.rdata'); ne=Inf
load('analysis/binom_sampling_simsum_2011-2014_ne3000.rdata'); ne=3000
load('analysis/binom_sampling_simsum_2011-2014_ne500.rdata'); ne=500

	# probability of observing as many or more loci in each starting freq & abs diff category
pvals <- as.matrix(tab)
pvals[,] <- NA

for(i in 1:ncol(tab)){
	f1nm <- colnames(tab)[i]
	for(j in 1:nrow(tab)){
		absnm <- rownames(tab)[j]
		if(tab[j,i]>0)	pvals[j,i] <- sum(simsum[f1nm,absnm,]>=tab[j,i])/dim(simsum)[3]
	}
}
pvals

pvals.adj <- matrix(data=p.adjust(pvals, method='fdr'), ncol=ncol(pvals), dimnames=list( absr=rownames(pvals), f1r=colnames(pvals)))

print(pvals.adj, digits=2)

pvalsstar <- pvals.adj
pvalsstar[pvals.adj>0.05] <- ''
print(pvalsstar, digits=2, quote=FALSE)

# plots
	# plot of observed frequency differences vs. the simuluations
xs <- as.numeric(rownames(tab)) # the x-axis for plotting the observed data
xs2 <- as.numeric(dimnames(simsum)$absr)
quartz(width=4, height=6)
# pdf(width=4, height=6, file=paste('analysis/figures/binom_samp_sims_vs_obs_ne', ne, '.pdf', sep=''))
par(mfrow=c(1/limstep,1), mai=c(0.1, 0.2, 0.1, 0.1), mgp=c(2.5,0.5,0), las=1, tcl=-0.2, omi=c(0.2,0.2,0,0))
for(i in 1:(1/limstep)){
	plot(xs, log10(tab[,i]+1), ylab='', xlab='', cex.axis=0.7, type='n') # set up the figure
	for(j in 1:nsims) lines(xs2, log10(simsum[i,,j]+1), col=rgb(0.5, 0.5, 0.5, 0.05))
	text(0.6,par('usr')[4]*0.85,labels=colnames(tab)[i], cex=0.8) # print initial allele freq in the upper right
	points(xs, log10(tab[,i]+1)) # plot to put obs on top
}
mtext('log10(# loci)',side=2,line=0.3, las=0,outer=TRUE, cex=0.8)
mtext('abs diff',side=1,line=0.5, outer=TRUE, cex=0.8)

dev.off()