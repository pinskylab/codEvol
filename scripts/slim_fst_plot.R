# Plot Fst from SLiM soft sweep simulation

require(PopGenome)
require(hierfstat)
require(data.table)

# settings
wndw <- 50000 # slding window size
step <- 10000 # step size
chrsize <- 30000000 # chromosome size

# compress and tabix so PopGenome can read
# NOTE: this requires that the path to bgzip and tabix match your system
# should set this up to first make a list of all .vcf files
# system('../TabixBinary/bgzip tmp/slim_sw_1607725495122_1.vcf')
# system('../TabixBinary/bgzip tmp/slim_sw_1607725495122_11.vcf')
# 
# system('../TabixBinary/tabix -p vcf tmp/slim_sw_1607725495122_1.vcf.gz')
# system('../TabixBinary/tabix -p vcf tmp/slim_sw_1607725495122_11.vcf.gz')
# 
# # read in
# p1 <- readVCF('tmp/slim_sw_1607725495122_1.vcf.gz', approx=FALSE, tid = '1', frompos = 1, topos = 99999, include.unknown=TRUE, numcols = 100000)
# p2 <- readVCF('tmp/slim_sw_1607725495122_11.vcf.gz', approx=FALSE, tid = '1', frompos = 1, topos = 99999, include.unknown=TRUE, numcols = 100000)

# read in data and merge populations together
p1 <- fread('tmp/slim_sw_1766425993417_1.vcf')
p2 <- fread('tmp/slim_sw_1766425993417_11.vcf')

setnames(p1, 1, 'CHROM') # fix column name
setnames(p2, 1, 'CHROM')
setnames(p1, grep('i[[:digit:]]*', names(p1)), paste0('t1_', grep('i[[:digit:]]*', names(p1), value = TRUE))) # append t1 to individual names for time 1
setnames(p2, grep('i[[:digit:]]*', names(p2)), paste0('t2_', grep('i[[:digit:]]*', names(p2), value = TRUE)))

keep1 <- c(2, grep('i[[:digit:]]*', names(p1)))
keep2 <- c(2, grep('i[[:digit:]]*', names(p2)))
p <- merge(p1[, ..keep1], p2[ , ..keep2], by = 'POS')

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
