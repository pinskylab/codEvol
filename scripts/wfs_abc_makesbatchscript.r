# write out a shell script that will submit an sbatch job for every combination of sample sizes that we need

if(!grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table)
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
}


load('analysis/wfs_targ.rdata') # targets normalized
load('analysis/wfs_sims_meansds.rdata') # useful as a list of sample sizes

# how many loci at each sample size
setkey(targ, alcnt1, alcnt2)
nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2)]
	nrow(nloci) # 185... why not 193?
	nrow(meansds)
	setdiff(meansds[,paste(alcnt1, alcnt2)], nloci[,paste(alcnt1, alcnt2)]) # 8 sample sizes not in targ (maybe because LG01, 02, etc were trimmed out?)

# sample sizes with >100 loci
nloci[nloci>5000,]
	
# write out the shell script
outfile <- 'scripts/wfs_abc_submit_all_sbatch.sh'

cat('#!/bin/bash\n', file=outfile, append=FALSE) # header

for(i in 1:nrow(nloci)){
	cat(paste('sbatch --job-name=abc', nloci[i,alcnt1], ',', nloci[i,alcnt2], ' --cpus-per-task=', min(16, nloci[i,nloci]), ' scripts/wfs_abc_sbatch.sh ', nloci[i,alcnt1], ' ', nloci[i, alcnt2], ' ', min(16, nloci[i,nloci]), '\n', sep=''), file=outfile, append=TRUE)
}