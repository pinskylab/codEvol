## examine results from site reshuffled LD from ANGSD genotypes
## run after ngsLD_siteshuffle_null.sh

###########################
# load functions
###########################
require(data.table)
require(ggplot2)
require(RColorBrewer)


#####################
# read in and prep data
#####################

# max FST per genome from reshuffling (GATK sites)
nullcangatk <- fread('analysis/ld.siteshuffle.Can_40.Can_14.gatk.csv.gz')
nulllof0711gatk <- fread('analysis/ld.siteshuffle.Lof_07.Lof_11.gatk.csv.gz')
nulllof0714gatk <- fread('analysis/ld.siteshuffle.Lof_07.Lof_14.gatk.csv.gz')
nulllof1114gatk <- fread('analysis/ld.siteshuffle.Lof_11.Lof_14.gatk.csv.gz')


# sliding window LD changes from ngsLD
# header is missing the fst column, so have to skip and make our own
bins <- fread('analysis/ld_change_region_5e4_ngsLD.gatk.csv.gz', drop = 1); width='5e4'




#######################
## null model stats
#######################

nullcangatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0711gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0714gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof1114gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]

###########################
# calc p-values per site
###########################
calcp <- function(ld, null){
  return((sum(null > ld)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen
}

bins[, p_can := calcp(ld_diff_Can, nullcangatk$x), by = .(chr, pos1mid1)]
bins[, p_lof0711 := calcp(ld_diff_Lof0711, nulllof0711gatk$x), by = .(chr, pos1mid1)]
bins[, p_lof0714 := calcp(ld_diff_Lof0714, nulllof0714gatk$x), by = .(chr, pos1mid1)]
bins[, p_lof1114 := calcp(ld_diff_Lof1114, nulllof1114gatk$x), by = .(chr, pos1mid1)]

##########################
# change to long format
##########################

binsl <- melt(bins[, .(chr, mid = pos1mid1, POSgen, p_can, p_lof0711, p_lof0714, p_lof1114)], 
              id.vars = c('chr', 'mid', 'POSgen'),
              variable.name = 'pop', value.name = 'p')

##############
# plots
##############


# plot p-value vs. position (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(binsl, aes(POSgen, -log10(p), color = chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p2

ggsave(plot = p2, device = 'png', filename = 'figures/ld_change.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)


# plot p-value vs. Nsites (gatk loci)
ggplot(datgatk, aes(Nsites, -log10(p), color = pop)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')


#################
# print outliers
#################


datgatk[p < 0.05, ]
