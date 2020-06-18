# run OutFLANK on slim sims

library(data.table)
library(OutFLANK)
library(ggplot2)

###########
# functions
###########

# count the number of alternative alleles at each locus
countalleles <- function(txt){
  if(length(dim(txt))>1){
    len <- nrow(txt)
  } else {
    len <- length(txt)
  }
  out <- rep(NA_real_, len)
  out[which(txt == '0|0')] <- 0
  out[which(txt == '0|1')] <- 1
  out[which(txt == '1|0')] <- 1
  out[which(txt == '1|1')] <- 2
  return(out)  
}

################
# prep s=0 files
################
files <- list.files('analysis/slim_sim', pattern = 'slim_sim_n[013]*_s0_[[:alnum:][:punct:]]*_i1_1.vcf$') # get the first iteration VCF file for the first generation for s=0
length(files)

for(i in 1:length(files)){
  print(i)
  thesefiles <- list.files('analysis/slim_sim', pattern = gsub('_i1_1', '_i[[:alnum:]]*_1', files[i]), full.names = TRUE)

  for(k in 1:length(thesefiles)){
    dat1 <- fread(thesefiles[k])
    file2 <- gsub('_1.vcf$', '_11.vcf', thesefiles[k]) # name of 2nd generation file
    dat2 <- fread(file2)
    
    setnames(dat1, '#CHROM', 'CHROM')
    setnames(dat2, '#CHROM', 'CHROM')
    
    # remove info columns
    dat1[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]
    dat2[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]
    
    # remove any duplicated positions (can happen in SLiM sims)
    dat1 <- dat1[!duplicated(cbind(CHROM, POS)), ]
    dat2 <- dat2[!duplicated(cbind(CHROM, POS)), ]
    
    # calculate allele frequencies
    for(j in 3:ncol(dat1)){
      counts <- countalleles(dat1[, ..j])    
      dat1[, (names(dat1)[j]) := counts]
      setnames(dat1, j, paste0(names(dat1)[j], '_1'))
    }
    
    for(j in 3:ncol(dat2)){
      counts <- countalleles(dat2[, ..j])    
      dat2[, (names(dat2)[j]) := counts]
      setnames(dat2, j, paste0(names(dat2)[j], '_11'))
    }
    
    
    # merge together on loci that match
    dat <- merge(dat1, dat2, by = c('CHROM', 'POS')) # append the individuals from dat2

    snpmat <- t(as.matrix(dat[,3:ncol(dat)]))
    colnames(snpmat) <- dat[, .(paste0(CHROM, '_', POS)), by = .(CHROM, POS)][,as.character(V1)]
    colnames(snpmat) <- gsub('^1', k+1, colnames(snpmat))
    
    if(k ==1) outmat <- snpmat
    if(k > 1) outmat <- cbind(outmat, snpmat)
  }
  
  outfile <- gsub('analysis/slim_sim/|_i[[:alnum:]]*_1*', '', thesefiles[k])
  outfile <- paste0('tmp/', gsub('.vcf', '_comb.rds', outfile))
  saveRDS(outmat, file = outfile)
  
}

######################################
# read in files
# and run outflank test on each one
######################################
files <- list.files('analysis/slim_sim', pattern = 'slim_sim_n[[:alnum:][:punct:]]*_1.vcf$',
                    full.names = TRUE) # get the VCF file for the first generation
length(files)


for(i in 1:length(files)){ # slow for larger Nes
  if(i %% 100 == 0) cat(paste0(i, ' '))
  dat1 <- fread(files[i])
  file2 <- gsub('_1.vcf$', '_11.vcf', files[i]) # name of 2nd generation file
  dat2 <- fread(file2)
  
  setnames(dat1, '#CHROM', 'CHROM')
  setnames(dat2, '#CHROM', 'CHROM')
  
  # remove info columns
  dat1[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]
  dat2[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]
  
  # remove any duplicated positions (can happen in SLiM sims)
  dat1 <- dat1[!duplicated(cbind(CHROM, POS)), ]
  dat2 <- dat2[!duplicated(cbind(CHROM, POS)), ]
  
  # calculate allele frequencies
  for(j in 3:ncol(dat1)){
      counts <- countalleles(dat1[, ..j])    
      dat1[, (names(dat1)[j]) := counts]
      setnames(dat1, j, paste0(names(dat1)[j], '_1'))
  }

  for(j in 3:ncol(dat2)){
    counts <- countalleles(dat2[, ..j])    
    dat2[, (names(dat2)[j]) := counts]
    setnames(dat2, j, paste0(names(dat2)[j], '_11'))
  }

  
  # merge together on loci that match
  dat <- merge(dat1, dat2, by = c('CHROM', 'POS')) # append the individuals from dat2
  rm(dat1, dat2)

  snpmat <- t(as.matrix(dat[,3:ncol(dat)]))
  locusNames <- dat[, .(paste0(CHROM, '_', POS)), by = .(CHROM, POS)][,as.character(V1)]
  colnames(snpmat) <- locusNames
  popNames <- ifelse(grepl('_11', rownames(snpmat)), '11', '1')
  
  if(i == 1 & exists('outflpow')) rm(outflpow) # remove output if it exists and this is the first iteration
  
  # run calculations on single chromosomes if enough loci
  if(nrow(dat) > 1000) {
    
    fstmat <- MakeDiploidFSTMat(snpmat, locusNames, popNames) # takes a few minutes
    
    # run outflank
    outfl <- OutFLANK(FstDataFrame=fstmat, LeftTrimFraction=0.05, RightTrimFraction=0.05, 
                      Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
    
    # summarize
    tempsum <- data.table(noutl = outfl$numberHighFstOutliers,
                          nloc = ncol(snpmat))
    tempsum[, ne := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n|_s.*', '', files[i]))] # get Ne from file name
    tempsum[, s := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n[[:digit:]]*_s|_f.*', '', files[i]))] # get s
    tempsum[, f := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n[[:digit:]]*_s[[:digit:].]*_f|_i.*', '', files[i]))] # get min f
    tempsum[, i := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\_1.vcf', '', files[i]))] # get iteration
    tempsum[, comb := 'onechrom']
    
    if(!exists('outflpow')){
      outflpow <- tempsum # if doesn't exist, this iteration will initialize the output
    } else {
      outflpow <- rbind(outflpow, tempsum) # otherwise, append
    }
  }

  # prep full genome (add 20 s=0 chromosomes)
  s0name <- gsub('analysis/slim_sim', 'tmp', files[i])
  s0name <- gsub('_i[[:alnum:]]*_1*.vcf', '_comb.rds', s0name)
  s0name <- gsub('_s[[:digit:][:punct:]]*_', '_s0_', s0name)
  s0sims <- readRDS(s0name)
  
  snpmatgen <- cbind(snpmat, s0sims)
  locusNamesgen <- colnames(snpmatgen)
  
  # run calculations on whole genome if enough loci
  if(ncol(snpmatgen) > 1000) {
    
    fstmatgen <- MakeDiploidFSTMat(snpmatgen, locusNamesgen, popNames)
    
    # run outflank
    outflgen <- OutFLANK(FstDataFrame=fstmatgen, LeftTrimFraction=0.05, RightTrimFraction=0.05, 
                      Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
    
    # summarize
    tempsum <- data.table(noutl = outflgen$numberHighFstOutliers,
                          nloc = ncol(snpmatgen))
    tempsum[, ne := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n|_s.*', '', files[i]))] # get Ne from file name
    tempsum[, s := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n[[:digit:]]*_s|_f.*', '', files[i]))] # get s
    tempsum[, f := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n[[:digit:]]*_s[[:digit:].]*_f|_i.*', '', files[i]))] # get min f
    tempsum[, i := as.numeric(gsub('analysis\\/slim_sim\\/slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\_1.vcf', '', files[i]))] # get iteration
    tempsum[, comb := 'genome']
    
    if(!exists('outflpow')){
      outflpow <- tempsum # if doesn't exist, this will initialize the output
    } else {
      outflpow <- rbind(outflpow, tempsum) # otherwise, append
    }
  }
}

# write out results
write.csv(outflpow, file = 'analysis/slim_outflank.csv', row.names = FALSE)


#######################################################
# walk through careful check of one sim
# see https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html
# assumes s0 files are still in tmp/ (from previous code)
#####################################################
library(bigsnpr)

# parameters
ne = 10000
s = 0
f = 0.05
i = 1
comb = 'genome'

# read in data 
infile1 <- paste0('analysis/slim_sim/slim_sim_n', ne, '_s', s, '_f', f, '_i', i, '_1.vcf')
dat1 <- fread(infile1)
infile2 <- gsub('_1.vcf$', '_11.vcf', infile1) # name of 2nd generation file
dat2 <- fread(infile2)

setnames(dat1, '#CHROM', 'CHROM')
setnames(dat2, '#CHROM', 'CHROM')

# remove info columns
dat1[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]
dat2[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]

# remove any duplicated positions (can happen in SLiM sims)
dat1 <- dat1[!duplicated(cbind(CHROM, POS)), ]
dat2 <- dat2[!duplicated(cbind(CHROM, POS)), ]

# calculate allele frequencies
for(j in 3:ncol(dat1)){
  counts <- countalleles(dat1[, ..j])    
  dat1[, (names(dat1)[j]) := counts]
  setnames(dat1, j, paste0(names(dat1)[j], '_1'))
}

for(j in 3:ncol(dat2)){
  counts <- countalleles(dat2[, ..j])    
  dat2[, (names(dat2)[j]) := counts]
  setnames(dat2, j, paste0(names(dat2)[j], '_11'))
}


# merge together on loci that match
dat <- merge(dat1, dat2, by = c('CHROM', 'POS')) # append the individuals from dat2
rm(dat1, dat2)

snpmat <- t(as.matrix(dat[,3:ncol(dat)]))
colnames(snpmat) <- dat[, paste0(CHROM, "_", POS)]
popNames <- ifelse(grepl('_11', rownames(snpmat)), '11', '1')

if(comb == 'genome'){
  # prep full genome (add 20 s=0 chromosomes)
  s0name <- gsub('analysis/slim_sim', 'tmp', infile)
  s0name <- gsub('_i[[:alnum:]]*_1*.vcf', '_comb.rds', s0name)
  s0name <- gsub('_s[[:digit:][:punct:]]*_', '_s0_', s0name)
  s0sims <- readRDS(s0name)
  snpmat <- cbind(snpmat, s0sims)
}
dim(snpmat)
locusNames <- colnames(snpmat)
chroms <- as.integer(sapply(strsplit(locusNames, split = '_'), '[', 1))
pos <- as.integer(sapply(strsplit(locusNames, split = '_'), '[', 2))

# calc FST
fstmat <- MakeDiploidFSTMat(snpmat, locusNames, popNames)
  
# data check: He vs. Fst
# "Here, you can see how some of the low H loci have high FST. 
# These are all neutral loci in the simulation, 
# and it is important to exclude them from the OutFLANK algorithm."
plot(fstmat$He, fstmat$FST)

# Data checks: FST vs. FSTNoCorr
# "Look for loci that deviate from the linear relationship in this plot, and remove those loci."
plot(fstmat$FST, fstmat$FSTNoCorr)
abline(0,1)

# Data prep: decide which SNPs to use for calibrating the null distribution of Fst
# uses bigsnpr
G<-add_code256(big_copy(snpmat,type="raw"),code=bigsnpr:::CODE_012)
newpc<-snp_autoSVD(G = G, infos.chr = chroms, infos.pos = pos, thr.r2 = 0.5)# prune to r2 < 0.5
which_pruned <- attr(newpc, which="subset") # Indexes of remaining SNPS after pruning
length(which_pruned)

# run outflank
outfl <- OutFLANK(FstDataFrame=fstmat[which_pruned, ], LeftTrimFraction=0.05, RightTrimFraction=0.05, 
                  Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
  
#Check the fit and make sure it looks good, especially in the right tail:
OutFLANKResultsPlotter(outfl, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)

## Zoom in on right tail
OutFLANKResultsPlotter(outfl, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = TRUE, RightZoomFraction = 0.15, titletext = NULL)

# check the P-value histogram:
# we expect this histogram to be flat and maybe have a bump near 0 for selected sites
hist(outfl$results$pvaluesRightTail)

# Using estimated neutral mean FST and df to calculate P-values for all loci
P1 <- pOutlierFinderChiSqNoCorr(fstmat, Fstbar = outfl$FSTNoCorrbar, 
                                dfInferred = outfl$dfInferred, qthreshold = 0.05, Hmin = 0.1)
my_out <- P1$OutlierFlag==TRUE
sum(my_out, na.rm = TRUE) # number of outlier loci (all loci)
outfl$numberHighFstOutliers # number of outlier loci (LD-trimmed loci)

# plot the outliers
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")

hist(P1$pvaluesRightTail) # p-values

####################################################
# analyze and plot the results from the Slim sims
####################################################
outflpow <- fread('analysis/slim_outflank.csv')

# summarize power
outflpowsum <- outflpow[, .(foutl = sum(noutl>0)/.N, noutl = mean(noutl), 
                            seoutl = sd(noutl)/sqrt(.N), nsims = .N), 
                        by = .(ne, s, f, comb)]
outflpowsum

# false positive rates
outflpowsum[s==0 & comb=='genome', ]

# plot fraction of sims with >=1 outlier detected
ggplot(outflpowsum, aes(s, foutl, group = as.factor(f), color = as.factor(f))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(comb ~ ne) +
  coord_cartesian(ylim = c(0,1)) +
  labs(y='Fraction of simulations with >0 outliers')
ggsave('figures/slim_outflank_foutl.png', width = 7, height = 4, dpi = 150)


# plot number of outliers
ggplot(outflpowsum, aes(s, noutl, group = as.factor(f), color = as.factor(f))) +
  geom_point() +
  geom_errorbar(aes(ymin=noutl-seoutl, ymax=noutl+seoutl)) + 
  geom_smooth(method = 'lm') +
  facet_grid(comb ~ ne) +
  labs(y='Number of outliers')
ggsave('figures/slim_outflank_noutl.png', width = 7, height = 4, dpi = 150)




#####################
# clean up tmp files
#####################


