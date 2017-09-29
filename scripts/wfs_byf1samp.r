## Function to do Wright-Fisher simulation with selection and a finite sample size
# using priors on Ne and s
# and specify initial sample allele frequency, f1samp (compare to wfs.r, which sets a prior on initial true allele freq f1)
# Sample at two time points

# f1samp: the initial sample allele frequency. Initial true allele freq will be sampled with this constraint.
# smin: the minimum for the selection coefficient uniform prior
# smax:
# c1: the size of the first sample (in # chromosomes)
# c2: the size of the second sample (in # chromosomes)
# gen: the number of generations
# ne: a vector of possible Ne values. Simulation picks one to use.
# h: dominance effect

wfs_byf1samp <- function(f1samp=0.5, smin=-1, smax=1, c1=58, c2=48, gen=20, ne=100, h=0.5){ 
	tol <- .Machine$double.eps^0.5 # tolerane for equality. default tolerance in all.equal

	# make sure f1samp is possible given initial sample size (c1)
	# test whether difference from a possible sample size is less than tol
	if(!any(abs(f1samp - ((0:c1)/c1)) < tol)) stop(paste('f1samp=', f1samp, 'not possible given c1=', c1))

	# choose parameters for this simulation
		s <- runif(1, min=smin, max=smax) # selection coefficient
	if(length(ne)>1) thisne <- round(sample(ne, 1)) # ne has to be integer for binomial sampling
	if(length(ne)==1) thisne <- round(ne)

	# choose initial population frequency
	# brute force method: keep picking possible f1s until one produces the correct f1samp
	thisf1 <- -100
	while(abs(thisf1 - f1samp) >= tol){
		f1 <- runif(1, min=0, max=1) # pick a starting allele frequency
		thisf1 <- rbinom(1,c1,f1)/c1 # pick a first sample allele frequency that corresponds
	}
	
#	print(paste(f1, s, thisne))
	
	p <- f1 # current allele frequency
	waa <- 1+s # relative fitness of genotype AA
	wab <- 1+s*h
	wbb <- 1
	for(i in 1:gen){
		x <- (waa*p^2 + wab*p*(1-p))/(waa*p^2 + wab*2*p*(1-p) + wbb*(1-p)^2) # probability of sampling allele A, given selection
		p <- rbinom(1,thisne,x)/thisne
#		print(paste(x,p))
	}
	f2 <- p

	f2samp <- rbinom(1,c2,f2)/c2 # second sample allele frequency

	# calculate summary statistics (Foll et al. 2015 Molecular Ecology Resources)
#	stats <- sumstats(f1samp, f2samp, c1, c2, gen)
	
	# return values
#	out <- c(ne=thisne, f1=f1, s=s, gen=gen, f2=f2, f1samp=f1samp, f2samp=f2samp, fsdprime=stats[1], fsiprime=stats[2])
	out <- c(ne=thisne, f1=f1, s=s, gen=gen, f2=f2, f1samp=f1samp, f2samp=f2samp)
	return(out)
}


sumstats <- function(f1samp, f2samp, c1, c2, gen){

	x <- min(f1samp, 1-f1samp) # minor allele frequency
	y <- min(f2samp, 1-f2samp)
	z <- (x+y)/2
	if(z>0){ # as long as one MAF > 0
		fs <- (x-y)^2/(z*(1-z))
		nstar <- 2/(1/c1 + 1/c2) # harmonic mean sample size
	
		if(f1samp>f2samp){ # focal allele is decreasing
			fsdprime <- (fs*(1-1/(2*nstar))-2/nstar)/(gen*(1+fs/4)*(1-1/c2))
			fsiprime <- 0	
		}
		if(f1samp<f2samp){ # increasing
			fsdprime <- 0
			fsiprime <- (fs*(1-1/(2*nstar))-2/nstar)/(gen*(1+fs/4)*(1-1/c2))
		}
		if(f1samp==f2samp){ # staying the same
			fsdprime <- (fs*(1-1/(2*nstar))-2/nstar)/(gen*(1+fs/4)*(1-1/c2))
			fsiprime <- (fs*(1-1/(2*nstar))-2/nstar)/(gen*(1+fs/4)*(1-1/c2))
		}
	}
	if(z==0){ # otherwise, z==0 and Fs is undefined
		fsdprime <- NA
		fsiprime <- NA
	}

	return(c(fsdprime=fsdprime, fsiprime=fsiprime))
}