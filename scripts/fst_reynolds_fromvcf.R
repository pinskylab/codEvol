# function to calculate Fst components from VCF files
# method of Reynolds et al. 1983 Genetics to match ANGSD

# Run as
#	  Rscript fst_reynolds_fromvcf.R infile1 infile2 outfile
# infile1: name of input file (vcf or vcf.gz)
# infile2: second input file (vcf or vcf.gz)
# outfile: output table
# does nothing to check that arguments are appropriate

#########
# set up
#########

require(data.table)
options(scipen = 999) # (effectively) turn off scientific notation

opts <- commandArgs(trailingOnly = TRUE)
#print(length(opts))
infile1 <- opts[1]
infile2 <- opts[2]
outfile <- opts[3]

# for troubleshooting
# print(infile1)
# print(infile2)
# print(outfile)


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

###############
# read in data
###############

dat1 <- fread(infile1)
dat2 <- fread(infile2)
  
setnames(dat1, '#CHROM', 'CHROM')
setnames(dat2, '#CHROM', 'CHROM')

# remove info columns
dat1[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]
dat2[, ':='(ID = NULL, REF = NULL, ALT = NULL, QUAL = NULL, FILTER = NULL, INFO = NULL, FORMAT = NULL)]

# calculate allele frequencies
for(i in 3:ncol(dat1)){
  if(i == 3){
    counts <- countalleles(dat1[, ..i])    
  }
  if(i > 3){
    counts <- counts + countalleles(dat1[, ..i])    
  }
}
dat1[, n1 := ncol(dat1) - 2]
dat1[, p1 := counts/n1/2]


for(i in 3:ncol(dat2)){
  if(i == 3){
    counts2 <- countalleles(dat2[, ..i])    
  }
  if(i > 3){
    counts2 <- counts2 + countalleles(dat2[, ..i])    
  }
}
dat2[, n2 := ncol(dat2) - 2]
dat2[, p2 := counts2/n2/2]


# merge together on loci that match
dat <- merge(dat1[, .(CHROM, POS, n1, p1)], dat2[, .(CHROM, POS, n2, p2)], by = c('CHROM', 'POS')) # append the individuals from dat2

rm(dat1, dat2)

###############
# calculations
###############
dat[, alpha1 := 1 - p1^2 - (1-p1)^2]
dat[, alpha2 := 1 - p2^2 - (1-p2)^2]
dat[, a := 1/2*((p1-p2)^2 + ((1-p1) - (1-p2))^2) +  # same as A from angsd (numerator)
      (n1 + n2)*(n1*alpha1 + n2*alpha2)/(4*n1*n2*(n1 + n2 - 1))]
dat[, ab := 1/2*((p1-p2)^2 + ((1-p1) - (1-p2))^2) +  # same as B from angsd (denominator)
      (4*n1*n2 - n1 - n2)*(n1*alpha1 + n2*alpha2)/(4*n1*n2*(n1 + n2 - 1))]

dat[, fst := a/ab] # theta in Reynolds

############
# write out
############

write.csv(dat[, .(CHROM, POS, A = a, B = ab)], file = gzfile(outfile), row.names = FALSE)

