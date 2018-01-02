## Function to do Wright-Fisher simulation with selection and a finite sample size
# using priors on Ne and s
# and specify initial sample allele frequency, f1samp (compare to wfs.r, which sets a prior on initial true allele freq f1)
# Sample at THREE time points

# f1samp: the initial sample allele frequency. Initial true allele freq will be sampled with this constraint.
# smin: the minimum for the selection coefficient uniform prior
# smax:
# c1: the size of the first sample (in # chromosomes)
# c2: the size of the second sample (in # chromosomes)
# c3: the size of the third sample (in # chromosomes)
# gen1: the number of generations for sample c2
# gen2: the number of generations for sample c3
# ne: a vector of possible Ne values. Simulation picks one to use.
# h: dominance effect

wfs_byf1samp_3samps <- function(f1samp=0.5, smin=-1, smax=1, c1=58, c2=48, c3=48, gen1=10, gen2=20, ne=100, h=0.5){ 
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
	for(i in 1:gen2){
		x <- (waa*p^2 + wab*p*(1-p))/(waa*p^2 + wab*2*p*(1-p) + wbb*(1-p)^2) # probability of sampling allele A, given selection
		p <- rbinom(1,thisne,x)/thisne
#		print(paste(x,p))
		if(i==gen1){
			f2 <- p # take second sample
		}
	}
	f3 <- p # third sample

	f2samp <- rbinom(1,c2,f2)/c2 # second sample allele frequency
	f3samp <- rbinom(1,c3,f3)/c3 # third sample allele frequency


	# return values
#	out <- c(ne=thisne, f1=f1, s=s, gen=gen, f2=f2, f1samp=f1samp, f2samp=f2samp, fsdprime=stats[1], fsiprime=stats[2])
	out <- c(ne=thisne, f1=f1, s=s, gen1=gen1, gen2=gen2, f2=f2, f3=f3, f1samp=f1samp, f2samp=f2samp, f3samp=f3samp)
	return(out)
}