# test wfs simulations

require(foreach)
require(doParallel)
source('scripts/wfs.r')

# make cluster
cl <- makeCluster(3)
registerDoParallel(cl)


# probability of fixation of a neutral allele should equal the initial frequency
	ne <- 100
	gen=1000 # run for a long time: likely to fixaton (expected time is 4N)
	smin <- smax <- 0
	f1min <- f1max <- ne/2/ne # starting allele freq
	nsims <- 1000


	sims <- foreach(i=1:nsims, .combine=rbind) %dopar% {
		if(i %% 100 == 0) cat(i)
		wfs(i, f1min=f1min, f1max=f1max, smin=smin, smax=smax, c1=50, c2=50, gen=gen, ne=ne)
	}
	sims <- as.data.frame(sims)

	# examine drift variance and sampling variance
	par(mfrow=c(3,1))
	bks <- seq(0,1,length.out=20)
	hist(sims$f1samp, main='f1samp', breaks=bks, col='grey'); abline(v=f1min, col='red', lwd=2)
	hist(sims$f2, main='f2', breaks=bks, col='grey'); abline(v=f1min, col='red', lwd=2)
	hist(sims$f2samp-sims$f2, main='delta f2samp', breaks=seq(-0.5,0.5,length.out=20), col='grey'); abline(v=0, col='red', lwd=2)

	sum(sims$f2 %in% c(0,1))/nsims # fraction that reached fixation
	f1min # expected fraction that will reach f2=1
	sum(sims$f2 == 1)/sum(sims$f2 %in% c(0,1)) # should be about f1min. divide by only those that reached fixation
	mean(sims$f2) # mean final allele frequency (should also be about f1)


# probability of fixation of a non-neutral allele should equal s if initial frequency is 1/ne
# (usually written as 2s, but that's when fitness of waa is 1+s2. here it is 1+s)
# for large N and small s
	ne <- 1000
	gen=10000 # run for a long time: likely to fixation in 4N generations with drift
	smin <- smax <- 0.2
	f1min <- f1max <- 1/ne # starting allele freq
	nsims <- 1000


	sims <- foreach(i=1:nsims, .combine=rbind) %dopar% {
		if(i %% 100 == 0) cat(i)
		wfs(i, f1min=f1min, f1max=f1max, smin=smin, smax=smax, c1=ne, c2=ne, gen=gen, ne=ne)
	}
	sims <- as.data.frame(sims)


	sum(sims$f2 %in% c(0,1))/nsims # fraction that reached fixation: make gen larger until this =1
	smin # expected fraction that will reach f2=1 (s in my sims is the coefficient for the homozygote, not for the het)
	smin/(1-exp(-2*ne*smin)) # Eq. 2.22 in Graur & Li
	(1-exp(-smin))/(1-exp(-2*ne*smin)) # Eq. 2.21 in Graur & Li
	sum(sims$f2 == 1)/sum(sims$f2 %in% c(0,1)) # should be about f1min. divide by only those that reached fixation
	mean(sims$f2) # mean final allele frequency (should also be about f1)


## examine how much drift vs. sampling variance for large ne
	ne <- 10000 # like NEA cod
	gen=11 # like NEA cod
	smin <- smax <- 0
	f1min <- f1max <- ne/2/ne # starting allele freq
	nsims <- 1000


	sims <- foreach(i=1:nsims, .combine=rbind) %dopar% {
		if(i %% 100 == 0) cat(i)
		wfs(i, f1min=f1min, f1max=f1max, smin=smin, smax=smax, c1=50, c2=50, gen=gen, ne=ne)
	}
	sims <- as.data.frame(sims)

	# examine drift variance and sampling variance
	par(mfrow=c(3,1))
	bks <- seq(0,1,length.out=20)
	hist(sims$f1samp, main='f1samp', breaks=bks, col='grey'); abline(v=f1min, col='red', lwd=2)
	hist(sims$f2, main='f2', breaks=bks, col='grey'); abline(v=f1min, col='red', lwd=2)
	hist(sims$f2samp-sims$f2, main='delta f2samp', breaks=seq(-0.5,0.5,length.out=20), col='grey'); abline(v=0, col='red', lwd=2)



# clean up cluster
stopCluster(cl)

