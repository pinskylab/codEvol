# summarize the output of pi_by_indiv.sh

library(data.table,  lib.loc="/projects/cees/lib/R_packages/")
library(reshape2,  lib.loc="/projects/cees/lib/R_packages/")
library(RColorBrewer,  lib.loc="/projects/cees/lib/R_packages/")

base <- '/projects/cees/in_progress/historic_cod_malin/'

# read in pi data
pi <- fread(input = paste('zcat <', base, 'analysis/angsd/pi_by_indiv.tsv.gz', sep=''), na.strings='NA')
setnames(pi, '#Chromo', 'CHROM') # to match coverage data.table
setnames(pi, 'Pos', 'POS')

# read in pop data
pops <- data.frame(pop=character(0), Indiv=character(0))
popfiles <- list.files(path=paste(base, 'data/', sep=''), pattern='pop.*.txt')
for(i in 1:length(popfiles)){
	thispop <- gsub("pop_|.txt", "", popfiles[i])
	newpop <- fread(input = paste(base, 'data/', popfiles[i], sep=''), header=FALSE)
	newpop[,pop:=thispop]
	setnames(newpop, "V1", "Indiv")
	pops <- rbind(pops, newpop)
}
popnms <- pops[,sort(unique(pop))]


#############################
## Histograms of average pi 
#############################
# summarize by individual
pi.sum <- pi[,.(Watterson=mean(Watterson), Pairwise=mean(Pairwise)), by=Indiv]

# append popdata to pi summary data
setkey(pi.sum, Indiv)
setkey(pops, Indiv)
pi.sum2 <- pi.sum[pops]

# histogram of pi by population
pdf(width=5, height=14, file=paste(base, 'analysis/figures/pi_by_indiv_hist.pdf', sep=''))
par(mfrow=c(length(popnms), 1))
for(i in 1:length(popnms)){
	pi.sum2[pop==popnms[i],hist(Pairwise, col='grey', xlim=c(-20,0), breaks=seq(-20,0,length.out=30), main=popnms[i])]
}
dev.off()


#####################
## Pi vs. coverage
#####################
# read in coverage data
cov <- fread(input = paste(base, 'analysis/allindivsLG03.gdepth', sep=''))
covm <- melt(cov, id.var=c('CHROM', 'POS'), variable.name='Indiv', value.name='cov') # melt

# merge pi and coverage
setkey(pi, CHROM, POS, Indiv)
setkey(covm, CHROM, POS, Indiv)
picov <- pi[covm]
	dim(pi)
	dim(covm)
	dim(picov) # all match rows: good!
	
# merge picov and pop
setkey(picov, Indiv)
setkey(pops, Indiv)
picov <- picov[pops]
	
# plot all loci in all individuals
jpeg(units='in', res=200, width=4, height=4, file=paste(base, 'analysis/figures/pi_vs_cov_by_locusindiv.jpg', sep=''))
#picov[,plot(cov, Pairwise, xlab='Coverage', ylab='Nucleotide diversity')]
picov[sample(1:nrow(picov),100000),plot(cov, Pairwise, xlab='Coverage', ylab='Nucleotide diversity', pch=16, cex=0.5)]
dev.off()
	
# summarize by individual
picov.sum <- picov[,.(cov=mean(cov), Pairwise=mean(Pairwise)), by=.(Indiv,pop)]

# plot all individuals: pi vs. cov
cols <- brewer.pal(6, 'Dark2')

#jpeg(units='in', res=200, width=4, height=4, file=paste(base, 'analysis/figures/pi_vs_cov_by_indiv.jpg', sep=''))
pdf(width=4, height=4, file=paste(base, 'analysis/figures/pi_vs_cov_by_indiv.pdf', sep=''))
picov.sum[,plot(cov, Pairwise, xlab='Coverage', ylab='Average nucleotide diversity', pch=16, cex=0.5, main='By individual', col=cols[as.numeric(pop)])]
legend('bottomleft', legend=unique(picov.sum$pop), col=cols, pch=16, cex=0.5)
dev.off()

