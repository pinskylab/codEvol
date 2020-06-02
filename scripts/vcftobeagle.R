# simple function to convert VCF (from SLiM) to beagle format
# uses only genotypes information. assumes that only genotypes are in the vcf file.
# Run as
#	  Rscript vcftobeagle.R e infile1 (infile2) outfile
# e: genotype error rate
# infile1: name of input file. outfile will be the same, minus .vcf.gz, with .beagle.gz
# infile2 (optional): second input file to append
# outfile
# does nothing to check that arguments are appropriate

#########
# set up
#########

require(data.table)
options(scipen = 999) # (effectively) turn off scientific notation

opts <- commandArgs(trailingOnly = TRUE)
#print(length(opts))
e <- as.numeric(opts[1])
infile1 <- opts[2]
if(length(opts) == 3){
	outfile <- opts[3]
} 
if(length(opts) == 4){
	infile2 <- opts[3]
	outfile <- opts[4]
} 

# for troubleshooting
# print(e)
# print(infile1)
# if(exists('infile2')){
# 	print(infile2)
# }
# print(outfile)

###########
# functions
###########

# convert vcf genotype to beagle probabilities
# use same format as pcangsd for error rates, see https://github.com/Rosemeis/pcangsd/blob/master/reader.pyx
gttogp <- function(txt, e, nm){
  if(length(dim(txt))>1){
    len <- nrow(txt)
  } else {
    len <- length(txt)
  }
  out <- data.table(aa = rep(NA_real_, len), ab = rep(NA_real_, len), bb = rep(NA_real_, len))
  out[which(txt == '0|0'), c('aa', 'ab', 'bb') := list((1-e)*(1-e), 2*(1-e)*e, e*e)]
  out[which(txt == '0|1'), c('aa', 'ab', 'bb') := list((1-e)*e, (1-e)*(1-e) + e*e, (1-e)*e)]
  out[which(txt == '1|0'), c('aa', 'ab', 'bb') := list((1-e)*e, (1-e)*(1-e) + e*e, (1-e)*e)]
  out[which(txt == '1|1'), c('aa', 'ab', 'bb') := list(e*e, 2*(1-e)*e, (1-e)*(1-e))]

  setnames(out, 1:3, rep(nm, 3))
  return(out)  
}

###############
# read in data
###############

dat <- fread(infile1)
if(exists('infile2')){
	dat2 <- fread(infile2)
	
	# append _1 to individual names in first table, _11 in second
	setnames(dat, 10:ncol(dat), paste0(names(dat)[10:ncol(dat)], '_1'))
	setnames(dat2, 10:ncol(dat2), paste0(names(dat2)[10:ncol(dat2)], '_11'))

	# merge together on loci that match
	colnms <- colnames(dat2)[c(1,2, 10:ncol(dat2))] # just the position and genotype columns
	dat <- merge(dat, dat2[, ..colnms], by = c('#CHROM', 'POS')) # append the individuals from dat2
}



###############
# calculations
###############

out <- data.table(marker = dat[, paste(`#CHROM`, POS, sep = '_')], allele1 = 1, allele2 = 2)
  
# for each individual in the vcf file (remove 9 cols of locus data on left of vcf table)
colnms <- colnames(dat)
for(i in 10:ncol(dat)){
  out <- cbind(out, gttogp(dat[, get(colnms[i])], e, colnms[i]))
}


############
# write out
############

write.table(out, file = gzfile(outfile), quote = FALSE, sep = '\t', row.names = FALSE)

