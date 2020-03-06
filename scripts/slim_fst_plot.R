# Make SLiM soft sweep simulations
# then plot Fst from sliding-window

require(data.table)


###########
# functions
###########

runslim <- function(L=30000000, ne=50, f=0.1, n1=22, n2=22, s=1.0, r=3.11e-8, g=11, o='tmp/test', i='tmp/test.vcf'){
  # notice complicated quoting needed to pass a string to slim
  # use -long to enable verbose output for troubleshooting
  system2("slim", c(paste0("-d ", c(paste0("L=", L), paste0("ne=", ne), paste0("f=", f), paste0("n1=", n1), 
                                    paste0("n2=", n2), paste0("s=", s), paste0("r=", r), 
                                    paste0("g=", g), paste0("\"o='", o, "'\""), paste0("\"i='", i, "'\""))), 
                    '-long', " scripts/slim_softsweep.slim"))
}

###########
# run slim
###########

bifiles <- list.files('analysis/', pattern = 'slim_burnin_n') # get the burn-in files
nes <- as.numeric(gsub('slim_burnin_n|_[12345]\\.vcf', '', bifiles)) # find the ne values from the file names
iter <- as.numeric(gsub('slim_burnin_n[[:digit:]]*_|\\.vcf', '', bifiles)) # find the iteraction #s from the file names

# Lofoten 07-11 parameters
paramtable <- expand.grid(iter = sort(unique(iter)), ne = sort(unique(nes)), f = c(0.1, 0.01), n1 = 22, n2 = 23, s = c(0.01, 0.03, 0.1, 0.3, 1), 
                          r = 3.11e-8, g = 11, L = 3e7)
dim(paramtable)

# run slim for Lofoten
for(i in 1:nrow(paramtable)){
  infile = with(paramtable[i,], paste0('analysis/slim_burnin_n', ne, '_', iter, '.vcf'))
  outfile = with(paramtable[i,], paste0('analysis/slim_sim_n', ne, '_s', s, '_f', f, '_i', iter))
  with(paramtable[i,], runslim(L = L, ne = ne, f = f, n1 = n1, n2 = n2, s = s, r = r, g = g, o = outfile, i = infile))
  
}



#############
# calculate fst
#############






# read in simulation output and merge populations together
p1 <- fread('tmp/slim_sw_1766427176444_1.vcf') # Ne 5000
p2 <- fread('tmp/slim_sw_1766427176444_11.vcf')

setnames(p1, 1, 'CHROM') # fix column name
setnames(p2, 1, 'CHROM')
setnames(p1, grep('i[[:digit:]]*', names(p1)), paste0('t1_', grep('i[[:digit:]]*', names(p1), value = TRUE))) # append t1 to individual names for time 1
setnames(p2, grep('i[[:digit:]]*', names(p2)), paste0('t2_', grep('i[[:digit:]]*', names(p2), value = TRUE)))

keep1 <- c(2, grep('i[[:digit:]]*', names(p1)))
keep2 <- c(2, grep('i[[:digit:]]*', names(p2)))
p <- merge(p1[, ..keep1], p2[ , ..keep2], by = 'POS')
dim(p)

# reformat genotypes for hierfstat
nms <- colnames(p)
for(i in grep('t[[:digit:]]', nms)){
  thisnm <- nms[i]
  p[, (thisnm) := as.numeric(gsub('\\|', '', get(thisnm))) + 11]
}

# sliding window FST
slide <- data.table(starts = seq(0, chrsize - wndw, by = step), # base 0. start of sliding windows
  ends = seq(wndw - 1, chrsize, by = step), fst = NA_real_)
slide[, mids := rowMeans(cbind(starts, ends))]
inds <- grep('i[[:digit:]]*', names(p)) # column indices for individuals

length(slide$starts)
for(i in 1:length(slide$starts)){
  if(i %% 100 == 0) print(i)
  
  # extract loci within this window
  locs <- t(p[POS >= slide$starts[i] & POS <= slide$ends[i], ..inds])
  
  # add population id
  pop <- rep(1, nrow(locs))
  pop[grep('t2', rownames(locs))] <- 2
  
  # make hierfstat data.frame and calculate Fst
  locs <- as.data.frame(cbind(pop, locs)) # add pop column on the left
  slide$fst[i] <- wc(locs)$FST # Weir-Cockerham Fst calculation
}


# plot
slide[, plot(mids/1e6, fst, cex = 0.6, xlab = 'Position (Mb)')]
