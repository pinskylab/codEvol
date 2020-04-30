## examine results from site reshuffled FST from ANGSD genotypes
## run after angsd_fst_siteshuffle_null.sh

#################
# parameters
#################
minloci <- 2 # should match angsd_fst_siteshuffle_null.r
winsz <- 50000 # window size
winstp <- 10000 # window step

###########################
# load functions
###########################
require(data.table)
#require(plyr)
require(ggplot2)
require(RColorBrewer)

calcp <- function(fst, null) return((sum(null > fst)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen

#####################
# read in and prep data
#####################

# max FST per genome from reshuffling (all sites)
nullcan <- fread('analysis/Can_40.Can_14.fst.siteshuffle.csv.gz')
nulllof0711 <- fread('analysis/Lof_07.Lof_11.fst.siteshuffle.csv.gz')
nulllof0714 <- fread('analysis/Lof_07.Lof_14.fst.siteshuffle.csv.gz')
nulllof1114 <- fread('analysis/Lof_11.Lof_14.fst.siteshuffle.csv.gz')

# max FST per genome from reshuffling (GATK sites, combined in linkage blocks, >1 SNP per window)
nullcangatk <- fread('analysis/Can_40.Can_14.gatk.fst.siteshuffle.csv.gz')
nulllof0711gatk <- fread('analysis/Lof_07.Lof_11.gatk.fst.siteshuffle.csv.gz')
nulllof0714gatk <- fread('analysis/Lof_07.Lof_14.gatk.fst.siteshuffle.csv.gz')
nulllof1114gatk <- fread('analysis/Lof_11.Lof_14.gatk.fst.siteshuffle.csv.gz')


# sliding window FST A and B components from ANGSD, after collapsing to unlinked loci
# header is missing the fst column, so have to skip and make our own
# need to make the all loci AB files
# can <- fread('analysis/Can_40.Can_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) # all sites
# lof0711 <- fread('analysis/Lof_07.Lof_11.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
# lof0714 <- fread('analysis/Lof_07.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
# lof1114 <- fread('analysis/Lof_11.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 

cangatk <- fread('analysis/Can_40.Can_14.gatk.ldtrim.fst.AB.csv.gz', drop = 1) # gatk sites
lof0711gatk <- fread('analysis/Lof_07.Lof_11.gatk.ldtrim.fst.AB.csv.gz', drop = 1) 
lof0714gatk <- fread('analysis/Lof_07.Lof_14.gatk.ldtrim.fst.AB.csv.gz', drop = 1) 
lof1114gatk <- fread('analysis/Lof_11.Lof_14.gatk.ldtrim.fst.AB.csv.gz', drop = 1) 


# nucleotide position for the whole genome (start position for each chr)
chrmax <- fread('data/lg_length.csv')
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, chr)

# merge nucleotide position into the frequency files
setkey(can, chr)
can <- can[chrmax[, .(chr, start)], ]
can[, posgen := midPos + start]
can[,start := NULL]

setkey(lof0711, chr)
lof0711 <- lof0711[chrmax[, .(chr, start)], ]
lof0711[, posgen := midPos + start]
lof0711[,start := NULL]

setkey(lof0714, chr)
lof0714 <- lof0714[chrmax[, .(chr, start)], ]
lof0714[, posgen := midPos + start]
lof0714[,start := NULL]

setkey(lof1114, chr)
lof1114 <- lof1114[chrmax[, .(chr, start)], ]
lof1114[, posgen := midPos + start]
lof1114[,start := NULL]

setkey(cangatk, CHROM)
cangatk <- cangatk[chrmax[, .(CHROM = chr, start)], ]
cangatk[, posgen := POS + start]
cangatk[,start := NULL]

setkey(lof0711gatk, CHROM)
lof0711gatk <- lof0711gatk[chrmax[, .(CHROM = chr, start)], ]
lof0711gatk[, posgen := POS + start]
lof0711gatk[,start := NULL]

setkey(lof0714gatk, CHROM)
lof0714gatk <- lof0714gatk[chrmax[, .(CHROM = chr, start)], ]
lof0714gatk[, posgen := POS + start]
lof0714gatk[,start := NULL]

setkey(lof1114gatk, CHROM)
lof1114gatk <- lof1114gatk[chrmax[, .(CHROM = chr, start)], ]
lof1114gatk[, posgen := POS + start]
lof1114gatk[,start := NULL]

######################
# Calc windowed FST
######################

# create new columns as indices for windows
for(j in 1:(winsz/winstp)){
  cangatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
  lof0711gatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
  lof0714gatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
  lof1114gatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
}

# calc fst and # snps per window
for(j in 1:(winsz/winstp)){
  if(j ==1){
    canfstsgatk <- cangatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
    lof0711fstsgatk <- lof0711gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
    lof0714fstsgatk <- lof0714gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
    lof1114fstsgatk <- lof1114gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
  } 
  if(j > 1){
    canfstsgatk <- rbind(canfstsgatk, cangatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
    lof0711fstsgatk <- rbind(lof0711fstsgatk, lof0711gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
    lof0714fstsgatk <- rbind(lof0714fstsgatk, lof0714gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
    lof1114fstsgatk <- rbind(lof1114fstsgatk, lof1114gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
  } 
}


#######################
## null model stats
#######################

nullcan[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0711[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0714[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof1114[, .(max = max(x), u95 = quantile(x, probs = 0.95))]

nullcangatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0711gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0714gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof1114gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]

###########################
# calc p-values per window
###########################

can[, p := calcp(fst, nullcan$x), by = .(chr, midPos)]
lof0711[, p := calcp(fst, nulllof0711$x), by = .(chr, midPos)]
lof0714[, p := calcp(fst, nulllof0714$x), by = .(chr, midPos)]
lof1114[, p := calcp(fst, nulllof1114$x), by = .(chr, midPos)]

canfstsgatk[, p := calcp(fst, nullcangatk$x), by = .(CHROM, midPos)]
lof0711fstsgatk[, p := calcp(fst, nulllof0711gatk$x), by = .(CHROM, midPos)]
lof0714fstsgatk[, p := calcp(fst, nulllof0714gatk$x), by = .(CHROM, midPos)]
lof1114fstsgatk[, p := calcp(fst, nulllof1114gatk$x), by = .(CHROM, midPos)]


######################
# combine the datasets
######################
# all loci
lof0711[, pop := 'lof0711']
lof0714[, pop := 'lof0714']
lof1114[, pop := 'lof1114']
can[, pop := 'can']

dat <- rbind(lof0711, lof0714, lof1114, can)
nrow(dat)
dat

dat[, pop := factor(pop, levels = c('can', 'lof0711', 'lof0714', 'lof1114'))]


# gatk loci
canfstsgatk[, pop := 'can']
lof0711fstsgatk[, pop := 'lof0711']
lof0714fstsgatk[, pop := 'lof0714']
lof1114fstsgatk[, pop := 'lof1114']

datgatk <- rbind(canfstsgatk, lof0711fstsgatk, lof0714fstsgatk, lof1114fstsgatk)
nrow(datgatk)
datgatk

datgatk[, pop := factor(pop, levels = c('can', 'lof0711', 'lof0714', 'lof1114'))]

# remove NAs and negative windows
datgatk <- datgatk[!is.na(fst) & midPos > 0, ]

# sort by window
setkey(datgatk, pop, CHROM, midPos)

##############
# Write out
##############

write.csv(dat, file = gzfile('output/fst_siteshuffle.angsd.csv.gz'), row.names = FALSE)
write.csv(datgatk, file = gzfile('output/fst_siteshuffle.angsd.gatk.csv.gz'), row.names = FALSE)

##############
# plots
##############

# plot p-value vs. position (all loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p1 <- ggplot(dat, aes(posgen, -log10(p), color = chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p1

ggsave(plot = p1, device = 'png', filename = 'figures/fst.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)



# plot p-value vs. position (gatk loci)
# only where nloci >= minloci
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(datgatk[nloci >= minloci, ], aes(midPos, -log10(p), color = CHROM)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p2

ggsave(plot = p2, device = 'png', filename = 'figures/fst.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)


# plot p-value vs. nloci (gatk loci)
ggplot(datgatk, aes(nloci, -log10(p), color = pop)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') + 
  scale_x_log10()

#################
# print outliers
#################


datgatk[p < 0.05, ]
