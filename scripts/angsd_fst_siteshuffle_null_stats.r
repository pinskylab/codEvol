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
# read and prep in data
#####################

# max FST per genome from reshuffling sites
nullcan <- fread('analysis/Can_40.Can_14.fst.siteshuffle.csv.gz')
nulllof0711 <- fread('analysis/Lof_07.Lof_11.fst.siteshuffle.csv.gz')

# sliding window FSTs from ANGSD
# header is missing the fst column, so have to skip and make our own
can <- fread('analysis/Can_40.Can_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof11 <- fread('analysis/Lof_07.Lof_11.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof14 <- fread('analysis/Lof_07.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
lof1114 <- fread('analysis/Lof_11.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- can[,.(len=max(midPos)), by=chr]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, chr)

# merge nucleotide position into the frequency files
setkey(lof11, chr)
lof11 <- lof11[chrmax[, .(chr, start)], ]
lof11[, posgen := midPos + start]
lof11[,start := NULL]

setkey(lof14, chr)
lof14 <- lof14[chrmax[, .(chr, start)], ]
lof14[, posgen := midPos + start]
lof14[,start := NULL]

setkey(lof1114, chr)
lof1114 <- lof1114[chrmax[, .(chr, start)], ]
lof1114[, posgen := midPos + start]
lof1114[,start := NULL]

setkey(can, chr)
can <- can[chrmax[, .(chr, start)], ]
can[, posgen := midPos + start]
can[,start := NULL]





#######################
## null model stats
#######################

can[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
lof0711[, .(max = max(x), u95 = quantile(x, probs = 0.95))]


###########################
# calc p-values per site
###########################
calcp <- function(fst, null){
  return((sum(null > fst)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen
}

can[, p := calcp(fst, nullcan$x), by = .(chr, midPos)]
lof11[, p := calcp(fst, nulllof0711$x), by = .(chr, midPos)]
lof14[, p:= NA_real_]
lof1114[, p:= NA_real_]

######################
# combine the datasets
######################
lof11[, pop := 'lof11']
lof14[, pop := 'lof14']
lof1114[, pop := 'lof1114']
can[, pop := 'can']

dat <- rbind(lof11, lof14, lof1114, can)
nrow(dat)
dat

dat[, pop := factor(pop, levels = c('can', 'lof11', 'lof14', 'lof1114'))]



##############
# plots
##############

# plot p-value vs. position
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p1 <- ggplot(dat, aes(posgen, -log10(p), color = chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p1

ggsave(plot = p1, device = 'png', filename = 'figures/fst.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)




