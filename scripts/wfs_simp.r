## Function to do Wright-Fisher simulation with selection and a finite sample size
# using priors on Ne, s, and initial allele frequency
# Sample at two time points
# compared to wfs.r, does not calculate extra summary statistcs (fsd and fsi)

# i is a dummy argument so that it can be used with sapply and parSapply
# f1: the initial allele frequency
# s: the selection coefficient
# c1: the size of the first sample (in # chromosomes)
# c2: the size of the second sample (in # chromosomes)
# gen: the number of generations
# ne: Ne value
# h: dominance effect

wfs_simp <- function(i, f1=0.5, s=0, c1=50, c2=59, gen=20, ne=100, h=0.5){ 
	# choose parameters for this simulation
	
	p <- f1 # current allele frequency
	waa <- 1+s # relative fitness of genotype AA
	wab <- 1+s*h
	wbb <- 1
	for(i in 1:gen){
		x <- (waa*p^2 + wab*p*(1-p))/(waa*p^2 + wab*2*p*(1-p) + wbb*(1-p)^2) # probability of sampling allele A, given selection
		p <- rbinom(1,ne,x)/ne
	}
	f2 <- p

	f1samp <- rbinom(1,c1,f1)/c1 # first sample allele frequency
	f2samp <- rbinom(1,c2,f2)/c2 # second sample allele frequency
	
	# return values
	out <- c(f2=f2, f1samp=f1samp, f2samp=f2samp)
	return(out)
}
