## examine results from site reshuffled FST from ANGSD genotypes
## run after angsd_fst_siteshuffle_null.sh

###########################
# load functions
###########################
require(data.table)
require(plyr)
require(ggplot2)
require(RColorBrewer)


#####################
# read in and prep data
#####################

# max FST per genome from reshuffling (all sites)
nullcan <- fread('analysis/Can_40.Can_14.fst.siteshuffle.csv.gz')
nulllof0711 <- fread('analysis/Lof_07.Lof_11.fst.siteshuffle.csv.gz')
nulllof0714 <- fread('analysis/Lof_07.Lof_14.fst.siteshuffle.csv.gz')
nulllof1114 <- fread('analysis/Lof_11.Lof_14.fst.siteshuffle.csv.gz')

# max FST per genome from reshuffling (GATK sites)
nullcangatk <- fread('analysis/Can_40.Can_14.gatk.fst.siteshuffle.csv.gz')
nulllof0711gatk <- fread('analysis/Lof_07.Lof_11.gatk.fst.siteshuffle.csv.gz')
nulllof0714gatk <- fread('analysis/Lof_07.Lof_14.gatk.fst.siteshuffle.csv.gz')
nulllof1114gatk <- fread('analysis/Lof_11.Lof_14.gatk.fst.siteshuffle.csv.gz')


# sliding window FSTs from ANGSD
# header is missing the fst column, so have to skip and make our own
can <- fread('analysis/Can_40.Can_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof0711 <- fread('analysis/Lof_07.Lof_11.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof0714 <- fread('analysis/Lof_07.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof1114 <- fread('analysis/Lof_11.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 

# sliding windows from ANGSD (GATK sites)
cangatk <- fread('analysis/Can_40.Can_14.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof0711gatk <- fread('analysis/Lof_07.Lof_11.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof0714gatk <- fread('analysis/Lof_07.Lof_14.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof1114gatk <- fread('analysis/Lof_11.Lof_14.gatk.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 


# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- can[,.(len=max(midPos)), by=chr]
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

setkey(cangatk, chr)
cangatk <- cangatk[chrmax[, .(chr, start)], ]
cangatk[, posgen := midPos + start]
cangatk[,start := NULL]

setkey(lof0711gatk, chr)
lof0711gatk <- lof0711gatk[chrmax[, .(chr, start)], ]
lof0711gatk[, posgen := midPos + start]
lof0711gatk[,start := NULL]

setkey(lof0714gatk, chr)
lof0714gatk <- lof0714gatk[chrmax[, .(chr, start)], ]
lof0714gatk[, posgen := midPos + start]
lof0714gatk[,start := NULL]

setkey(lof1114gatk, chr)
lof1114gatk <- lof1114gatk[chrmax[, .(chr, start)], ]
lof1114gatk[, posgen := midPos + start]
lof1114gatk[,start := NULL]


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
# calc p-values per site
###########################
calcp <- function(fst, null){
  return((sum(null > fst)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen
}

can[, p := calcp(fst, nullcan$x), by = .(chr, midPos)]
lof0711[, p := calcp(fst, nulllof0711$x), by = .(chr, midPos)]
lof0714[, p := calcp(fst, nulllof0714$x), by = .(chr, midPos)]
lof1114[, p := calcp(fst, nulllof1114$x), by = .(chr, midPos)]

cangatk[, p := calcp(fst, nullcangatk$x), by = .(chr, midPos)]
lof0711gatk[, p := calcp(fst, nulllof0711gatk$x), by = .(chr, midPos)]
lof0714gatk[, p := calcp(fst, nulllof0714gatk$x), by = .(chr, midPos)]
lof1114gatk[, p := calcp(fst, nulllof1114gatk$x), by = .(chr, midPos)]


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
lof0711gatk[, pop := 'lof0711']
lof0714gatk[, pop := 'lof0714']
lof1114gatk[, pop := 'lof1114']
cangatk[, pop := 'can']

datgatk <- rbind(lof0711gatk, lof0714gatk, lof1114gatk, cangatk)
nrow(datgatk)
datgatk

datgatk[, pop := factor(pop, levels = c('can', 'lof0711', 'lof0714', 'lof1114'))]

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
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(datgatk, aes(posgen, -log10(p), color = chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p2

ggsave(plot = p2, device = 'png', filename = 'figures/fst.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)



