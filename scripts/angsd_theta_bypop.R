# average theta values across sites

#############
# Parameters
#############
# Total number of loci evaluated (not just SNPs)
# nloci <- max(c(nrow(datCan40), nrow(datCan14), nrow(datLof07), 
#                nrow(datLof11), nrow(datLof14))) # not sure how many loci got called, so use the max in any pop

# nlocigatk <- max(c(nrow(datCan40gatk), nrow(datCan14gatk), nrow(datLof07gatk), 
#                nrow(datLof11gatk), nrow(datLof14gatk))) # not sure how many loci got called, so use the max in any pop

nloci <- 1 # to allow after-the fact calculations
nlocigatk <- 1



####################
# load functions
require(data.table)
require(boot) # for bootstrap CIs


# calculate pi from specified LGs for block bootstrapping across LGs
# frin https://stackoverflow.com/questions/11919808/block-bootstrap-from-subject-list
piblock <- function(lgs, indices, alldata, nloci){
  mydata <- do.call("rbind", lapply(indices, function(n) subset(alldata, Chromo==lgs[n])))
  pi <- mydata[, sum(exp(Pairwise), na.rm = TRUE)/nloci]
  return(pi)
}


######################
# Load data
######################

# load all loci theta calcs
datCan40 <- fread('analysis/thetas.Can_40.pestPG.gz')
datCan14 <- fread('analysis/thetas.Can_14.pestPG.gz')
datLof07 <- fread('analysis/thetas.Lof_07.pestPG.gz')
datLof11 <- fread('analysis/thetas.Lof_11.pestPG.gz')
datLof14 <- fread('analysis/thetas.Lof_14.pestPG.gz')

# gatk loci
datCan40gatk <- fread('analysis/thetas.Can_40.gatk.pestPG.gz')
datCan14gatk <- fread('analysis/thetas.Can_14.gatk.pestPG.gz')
datLof07gatk <- fread('analysis/thetas.Lof_07.gatk.pestPG.gz')
datLof11gatk <- fread('analysis/thetas.Lof_11.gatk.pestPG.gz')
datLof14gatk <- fread('analysis/thetas.Lof_14.gatk.pestPG.gz')

# fix name
setnames(datCan40, '#Chromo', 'Chromo')
setnames(datCan14, '#Chromo', 'Chromo')
setnames(datLof07, '#Chromo', 'Chromo')
setnames(datLof11, '#Chromo', 'Chromo')
setnames(datLof14, '#Chromo', 'Chromo')

setnames(datCan40gatk, '#Chromo', 'Chromo')
setnames(datCan14gatk, '#Chromo', 'Chromo')
setnames(datLof07gatk, '#Chromo', 'Chromo')
setnames(datLof11gatk, '#Chromo', 'Chromo')
setnames(datLof14gatk, '#Chromo', 'Chromo')

# remove unplaced
datCan40 <- datCan40[grep('Unplaced', Chromo, invert = TRUE), ]
datCan14 <- datCan14[grep('Unplaced', Chromo, invert = TRUE), ]
datLof07 <- datLof07[grep('Unplaced', Chromo, invert = TRUE), ]
datLof11 <- datLof11[grep('Unplaced', Chromo, invert = TRUE), ]
datLof14 <- datLof14[grep('Unplaced', Chromo, invert = TRUE), ]

datCan40gatk <- datCan40gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datCan14gatk <- datCan14gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datLof07gatk <- datLof07gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datLof11gatk <- datLof11gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datLof14gatk <- datLof14gatk[grep('Unplaced', Chromo, invert = TRUE), ]

################################
# Run theta calculations
################################


datCan40[, sum(exp(Pairwise), na.rm = TRUE)/nloci]
datCan14[, sum(exp(Pairwise), na.rm = TRUE)/nloci]
datLof07[, sum(exp(Pairwise), na.rm = TRUE)/nloci]
datLof11[, sum(exp(Pairwise), na.rm = TRUE)/nloci]
datLof14[, sum(exp(Pairwise), na.rm = TRUE)/nloci]

datCan40gatk[, sum(exp(Pairwise), na.rm = TRUE)/nlocigatk]
datCan14gatk[, sum(exp(Pairwise), na.rm = TRUE)/nlocigatk]
datLof07gatk[, sum(exp(Pairwise), na.rm = TRUE)/nlocigatk]
datLof11gatk[, sum(exp(Pairwise), na.rm = TRUE)/nlocigatk]
datLof14gatk[, sum(exp(Pairwise), na.rm = TRUE)/nlocigatk]



# block bootstrapping across LGs
lgs <- datCan40[, sort(unique(Chromo))]

print('Can40 all loci')
bootCan40lg <- boot(lgs, piblock, 1000,  alldata = datCan40, nloci = nloci)
print(bootCan40lg)
print(boot.ci(bootCan40lg, type = c('norm', 'basic', 'perc')))

print('Can14 all loci')
bootCan14lg <- boot(lgs, piblock, 1000,  alldata = datCan14, nloci = nloci)
print(bootCan14lg)
print(boot.ci(bootCan14lg, type = c('norm', 'basic', 'perc')))

print('Lof07 all loci')
bootLof07lg <- boot(lgs, piblock, 1000,  alldata = datLof07, nloci = nloci)
print(bootLof07lg)
print(boot.ci(bootLof07lg, type = c('norm', 'basic', 'perc')))

print('Lof11 all loci')
bootLof11lg <- boot(lgs, piblock, 1000,  alldata = datLof11, nloci = nloci)
print(bootLof11lg)
print(boot.ci(bootLof11lg, type = c('norm', 'basic', 'perc')))

print('Lof14 all loci')
bootLof14lg <- boot(lgs, piblock, 1000,  alldata = datLof14, nloci = nloci)
print(bootLof14lg)
print(boot.ci(bootLof14lg, type = c('norm', 'basic', 'perc')))

print('Can40 gatk loci')
bootCan40lggatk <- boot(lgs, piblock, 1000,  alldata = datCan40gatk, nloci = nlocigatk)
print(bootCan40lggatk)
print(boot.ci(bootCan40lggatk, type = c('norm', 'basic', 'perc')))

print('Can14 gatk loci')
bootCan14lggatk <- boot(lgs, piblock, 1000,  alldata = datCan14gatk, nloci = nlocigatk)
print(bootCan14lggatk)
print(boot.ci(bootCan14lggatk, type = c('norm', 'basic', 'perc')))

print('Lof07 gatk loci')
bootLof07lggatk <- boot(lgs, piblock, 1000,  alldata = datLof07gatk, nloci = nlocigatk)
print(bootLof07lggatk)
print(boot.ci(bootLof07lggatk, type = c('norm', 'basic', 'perc')))

print('Lof11 gatk loci')
bootLof11lggatk <- boot(lgs, piblock, 1000,  alldata = datLof11gatk, nloci = nlocigatk)
print(bootLof11lggatk)
print(boot.ci(bootLof11lggatk, type = c('norm', 'basic', 'perc')))

print('Lof14 gatk loci')
bootLof14lggatk <- boot(lgs, piblock, 1000,  alldata = datLof14gatk, nloci = nlocigatk)
print(bootLof14lggatk)
print(boot.ci(bootLof14lggatk, type = c('norm', 'basic', 'perc')))
