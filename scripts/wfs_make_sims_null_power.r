# Make simulations for a power analysis
# to run on cod node

# set parameters
#pop <- 'Lof'; Ne <- 46000; c1 <- 46; c2 <- 44; gen=11 # see 7/5/2017 Ne calculation with wfabc_1
pop <- 'Can'; Ne<- 5900; c1 <- 32; c2 <- 40; gen <- 8 # see 8/25/2017

# load functions
source('scripts/wfs.r')
require(parallel, lib.loc="/projects/cees/lib/R_packages/")
require(data.table, lib.loc="/projects/cees/lib/R_packages/")

# parameters for the simulation
nsims <- 100 # how many sims at each level of ss and f1s
ss <- seq(0,1,length.out=21) # each selection coefficient to simulate
f1s <- c(0.001, 0.01, seq(0.05, 0.95, by=0.05), 0.99, 0.999) # each starting allele freq to simulate

# start cluster
cl <- makeCluster(30) # set up 30-core cluster on a cod node
clusterExport(cl, c('wfs'))


# make simulations (parallel way)
thisout<-matrix(NA, nrow=7,ncol=0)
for(i in 1:length(ss)){
	print(paste('i=', i, 'of', length(ss)))
	for(j in 1:length(f1s)){
		temp <- parSapply(cl, 1:nsims, FUN=wfs, f1min=f1s[j], f1max=f1s[j], smin=ss[i], smax=ss[i], c1=c1, c2=c2, gen=gen, ne=Ne, h=0.5, simplify=TRUE)
		thisout <- cbind(thisout, temp)
	}
}
	
# stop cluster
stopCluster(cl)

# how many
ncol(thisout)

# save to Frequency-style file
outfile<-paste('analysis/Frequency_table_PowerSims_', pop, '_Ne', Ne, '_cnt', c1, '_', c2,'.txt', sep='')
outfile

cat('CHROM\tPOS\tN_CHR\t{FREQ}\t{FREQ}\t\t0\n', file=outfile, append=FALSE) # header

for(i in 1:ncol(thisout)){
	if(i %% 1000 == 0) print(i)
	cat(paste('s=', thisout['s',i], ',f1=', thisout['f1',i], '\t', i, '\t', c1, '\t', round(thisout['f1samp',i],6), '\t', c2, '\t', round(thisout['f2samp',i],6), '\t', round(abs(thisout['f2samp',i]-thisout['f1samp',i]),6), '\n', sep=''), file=outfile, append=TRUE)
}

