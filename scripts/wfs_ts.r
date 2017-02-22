## Function to do Wright-Fisher simulation with selection and a finite sample size
# using set Ne, s, and initial allele frequency
# Sample each generation

# f1: the initial allele frequency uniform prior
# s: the selection coefficient uniform prior
# cnt: the size of the samples (in # chromosomes)
# gen: the number of generations
# ne: Ne value
# h: dominance effect

wfs_ts <- function(f1=0, s=0, cnt=50, gen=20, ne=100, h=0.5){ 
	p <- rep(NA, gen)
	fsamp <- p

	p[1] <- f1 # current allele frequency
	fsamp[1] <- rbinom(1,cnt,p[1])/cnt # first sample allele frequency
	waa <- 1+s # relative fitness of genotype AA
	wab <- 1+s*h
	wbb <- 1
	for(i in 2:gen){
		x <- (waa*p[i-1]^2 + wab*p[i-1]*(1-p[i-1]))/(waa*p[i-1]^2 + wab*2*p[i-1]*(1-p[i-1]) + wbb*(1-p[i-1])^2) # probability of sampling allele A, given selection
		p[i] <- rbinom(1,ne,x)/ne
		fsamp[i] <- rbinom(1,cnt,p[i])/cnt # first sample allele frequency
#		print(paste(x,p))
	}

	
	# return values
	out <- list(gen=1:gen, p=p, fsamp=fsamp)
	return(out)
}
