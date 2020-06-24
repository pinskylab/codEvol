## examine results from site reshuffled theta from ANGSD genotypes
## run after angsd_theta_siteshuffle_null.sh/.r

###########################
# load functions
###########################
require(data.table)
require(plyr)
require(ggplot2)
require(RColorBrewer)


##################
# settings
##################

width = '5e4' # window size, used for naming some of the output plots
winsz = 5e4 # for scaling the windowed pi and thetaW values





#####################
# read in and prep data
#####################

# max theta per genome from reshuffling (all sites) from angsd_theta_siteshuffle_null.r
# nullcan <- fread('analysis/theta.siteshuffle.Can_40.Can_14.csv.gz')
# nulllof0711 <- fread('analysis/theta.siteshuffle.Lof_07.Lof_11.csv.gz')
# nulllof0714 <- fread('analysis/theta.siteshuffle.Lof_07.Lof_14.csv.gz')
# nulllof1114 <- fread('analysis/theta.siteshuffle.Lof_11.Lof_14.csv.gz')

# max theta per genome from reshuffling (GATK sites) from angsd_theta_siteshuffle_null.r
nullcangatk <- fread('analysis/theta.siteshuffle.Can_40.Can_14.gatk.csv.gz')
nulllof0711gatk <- fread('analysis/theta.siteshuffle.Lof_07.Lof_11.gatk.csv.gz')
nulllof0714gatk <- fread('analysis/theta.siteshuffle.Lof_07.Lof_14.gatk.csv.gz')
nulllof1114gatk <- fread('analysis/theta.siteshuffle.Lof_11.Lof_14.gatk.csv.gz')


# sliding window theta change (all sites)
# dat <- fread('analysis/theta_change_region_5e4.csv.gz', drop = 1) 

# sliding windows theta change (GATK sites) from angsd_theta_siteshuffle_null.r
datcangatk <- fread('analysis/theta_change_region_50000.Can_40.Can_14.gatk.csv.gz', drop = 1)
datlof0711gatk <- fread('analysis/theta_change_region_50000.Lof_07.Lof_11.gatk.csv.gz', drop = 1)
datlof0714gatk <- fread('analysis/theta_change_region_50000.Lof_07.Lof_14.gatk.csv.gz', drop = 1)
datlof1114gatk <- fread('analysis/theta_change_region_50000.Lof_11.Lof_14.gatk.csv.gz', drop = 1)

# combine population files
datcangatk[, pop := 'can']
datlof0711gatk[, pop := 'lof0711']
datlof0714gatk[, pop := 'lof0714']
datlof1114gatk[, pop := 'lof1114']
datgatk <- rbind(datcangatk, datlof0711gatk, datlof0714gatk, datlof1114gatk)

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- fread('data/lg_length.csv')
chrmax$start=c(0,cumsum(chrmax$length)[1:(nrow(chrmax)-1)])

# merge nucleotide position into the frequency files
# setkey(dat, Chr)
# dat <- dat[chrmax[, .(Chr = chr, start)], ]
# dat[, POSgen := WinCenter + start]
# dat[,start := NULL]

datgatk <- merge(datgatk, chrmax[, .(Chromo = chr, start)], by = "Chromo")
datgatk[, POSgen := WinCenter + start]
datgatk[,start := NULL]


#######################
## null model stats
#######################

# nullcan[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95), # using 5% and 95% to be lenient
#             tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
#             tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
# nulllof0711[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
#             tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
#             tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
# nulllof0714[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
#                 tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
#                 tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
# nulllof1114[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
#                 tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
#                 tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]

nullcangatk[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
                tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
                tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
nulllof0711gatk[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
                tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
                tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
nulllof0714gatk[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
                tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
                tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
nulllof1114gatk[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
                tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
                tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]


###########################
# calc p-values per site
###########################
calcpG <- function(thetachange, null){ # for increases in theta
  return((sum(null > thetachange)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen
}
calcpL <- function(thetachange, null){ # for decreases in theta
  return((sum(null < thetachange)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen
}

# dat[pop == 'can' & tWd > 0, tWd.p := calcpG(tWd, nullcan$maxtWd), by = .(Chr, WinCenter)] # thetaW all loci
# dat[pop == 'can' & tWd <= 0, tWd.p := calcpL(tWd, nullcan$mintWd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0711' & tWd > 0, tWd.p := calcpG(tWd, nulllof0711$maxtWd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0711' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0711$mintWd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0714' & tWd > 0, tWd.p := calcpG(tWd, nulllof0714$maxtWd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0714' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0714$mintWd), by = .(Chr, WinCenter)]
# dat[pop == 'lof1114' & tWd > 0, tWd.p := calcpG(tWd, nulllof1114$maxtWd), by = .(Chr, WinCenter)]
# dat[pop == 'lof1114' & tWd <= 0, tWd.p := calcpL(tWd, nulllof1114$mintWd), by = .(Chr, WinCenter)]

# dat[pop == 'can' & tPd > 0, tPd.p := calcpG(tPd, nullcan$maxtPd), by = .(Chr, WinCenter)] # theta pi
# dat[pop == 'can' & tPd <= 0, tPd.p := calcpL(tPd, nullcan$mintPd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0711' & tPd > 0, tPd.p := calcpG(tPd, nulllof0711$maxtPd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0711' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0711$mintPd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0714' & tPd > 0, tPd.p := calcpG(tPd, nulllof0714$maxtPd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0714' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0714$mintPd), by = .(Chr, WinCenter)]
# dat[pop == 'lof1114' & tPd > 0, tPd.p := calcpG(tPd, nulllof1114$maxtPd), by = .(Chr, WinCenter)]
# dat[pop == 'lof1114' & tPd <= 0, tPd.p := calcpL(tPd, nulllof1114$mintPd), by = .(Chr, WinCenter)]

# dat[pop == 'can' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nullcan$maxtDd), by = .(Chr, WinCenter)] # tajima's D
# dat[pop == 'can' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nullcan$mintDd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0711' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof0711$maxtDd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0711' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof0711$mintDd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0714' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof0714$maxtDd), by = .(Chr, WinCenter)]
# dat[pop == 'lof0714' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof0714$mintDd), by = .(Chr, WinCenter)]
# dat[pop == 'lof1114' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof1114$maxtDd), by = .(Chr, WinCenter)]
# dat[pop == 'lof1114' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof1114$mintDd), by = .(Chr, WinCenter)]


datgatk[pop == 'can' & tWd > 0, tWd.p := calcpG(tWd, nullcangatk$maxtWd), by = .(Chromo, WinCenter)] # thetaW gatk loci
datgatk[pop == 'can' & tWd <= 0, tWd.p := calcpL(tWd, nullcangatk$mintWd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0711' & tWd > 0, tWd.p := calcpG(tWd, nulllof0711gatk$maxtWd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0711' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0711gatk$mintWd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0714' & tWd > 0, tWd.p := calcpG(tWd, nulllof0714gatk$maxtWd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0714' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0714gatk$mintWd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof1114' & tWd > 0, tWd.p := calcpG(tWd, nulllof1114gatk$maxtWd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof1114' & tWd <= 0, tWd.p := calcpL(tWd, nulllof1114gatk$mintWd), by = .(Chromo, WinCenter)]

datgatk[pop == 'can' & tPd > 0, tPd.p := calcpG(tPd, nullcangatk$maxtPd), by = .(Chromo, WinCenter)] # theta pi
datgatk[pop == 'can' & tPd <= 0, tPd.p := calcpL(tPd, nullcangatk$mintPd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0711' & tPd > 0, tPd.p := calcpG(tPd, nulllof0711gatk$maxtPd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0711' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0711gatk$mintPd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0714' & tPd > 0, tPd.p := calcpG(tPd, nulllof0714gatk$maxtPd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0714' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0714gatk$mintPd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof1114' & tPd > 0, tPd.p := calcpG(tPd, nulllof1114gatk$maxtPd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof1114' & tPd <= 0, tPd.p := calcpL(tPd, nulllof1114gatk$mintPd), by = .(Chromo, WinCenter)]

datgatk[pop == 'can' & tDd > 0, tDd.p := calcpG(tDd, nullcangatk$maxtDd), by = .(Chromo, WinCenter)] # tajima's D
datgatk[pop == 'can' & tDd <= 0, tDd.p := calcpL(tDd, nullcangatk$mintDd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0711' & tDd > 0, tDd.p := calcpG(tDd, nulllof0711gatk$maxtDd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0711' & tDd <= 0, tDd.p := calcpL(tDd, nulllof0711gatk$mintDd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0714' & tDd > 0, tDd.p := calcpG(tDd, nulllof0714gatk$maxtDd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof0714' & tDd <= 0, tDd.p := calcpL(tDd, nulllof0714gatk$mintDd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof1114' & tDd > 0, tDd.p := calcpG(tDd, nulllof1114gatk$maxtDd), by = .(Chromo, WinCenter)]
datgatk[pop == 'lof1114' & tDd <= 0, tDd.p := calcpL(tDd, nulllof1114gatk$mintDd), by = .(Chromo, WinCenter)]


##############
# Write out
##############

write.csv(datgatk, file = gzfile('analysis/theta_siteshuffle.angsd.gatk.csv.gz'), row.names = FALSE)

#########################
# plots of region change
#########################

# plot pi change (all loci)
# cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
# p1 <- ggplot(dat, aes(POSgen, tPd/5e4, color = Chr)) + 
#   geom_point(size = 0.5, alpha = 0.3) +
#   facet_wrap(~pop, ncol = 1, scales = 'free') +
#   scale_color_manual(values = cols) +
#   ylab('Change in pi per site')
# p1
# ggsave(plot = p1, device = 'png', filename = paste0('figures/pi_change_vs_pos_NEA_CAN_runmean', width, '.angsd.png'), 
#        width = 7.5, height = 6, units = 'in', dpi = 300)
# 
# # plot theta_Waterson change (all loci)
# cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
# p2 <- ggplot(dat, aes(POSgen, tWd/5e4, color = Chr)) + 
#   geom_point(size = 0.5, alpha = 0.3) +
#   facet_wrap(~pop, ncol = 1, scales = 'free') +
#   scale_color_manual(values = cols) +
#   ylab('Change in Wattersons theta per site')
# p2
# ggsave(plot = p2, device = 'png', filename = paste0('figures/thetaW_change_vs_pos_NEA_CAN_runmean', width, '.angsd.png'), 
#        width = 7.5, height = 6, units = 'in', dpi = 300)
# 
# # plot Tajima's D change (all loci)
# cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
# p3 <- ggplot(dat, aes(POSgen, Tajimadiff/5e4, color = Chr)) + 
#   geom_point(size = 0.5, alpha = 0.3) +
#   facet_wrap(~pop, ncol = 1, scales = 'free') +
#   scale_color_manual(values = cols) +
#   ylab('Change in Tajimas D per site')
# p3
# ggsave(plot = p3, device = 'png', filename = paste0('figures/TajimasD_change_vs_pos_NEA_CAN_runmean', width, '.angsd.png'), 
#        width = 7.5, height = 6, units = 'in', dpi = 300)


# plot pi change (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p1 <- ggplot(datgatk, aes(POSgen, tPd/winsz, color = Chromo)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols) +
  ylab('Change in pi per site')
p1
ggsave(plot = p1, device = 'png', filename = paste0('figures/pi_change_vs_pos_NEA_CAN_runmean', width, '.angsd.gatk.png'), 
       width = 7.5, height = 6, units = 'in', dpi = 300)

# plot theta_Waterson change (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(datgatk, aes(POSgen, tWd/winsz, color = Chromo)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols) +
  ylab('Change in Wattersons theta per site')
p2
ggsave(plot = p2, device = 'png', filename = paste0('figures/thetaW_change_vs_pos_NEA_CAN_runmean', width, '.angsd.gatk.png'), 
       width = 7.5, height = 6, units = 'in', dpi = 300)

# plot Tajima's D change (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p3 <- ggplot(datgatk, aes(POSgen, tDd, color = Chromo)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols) +
  ylab('Change in Tajimas D per window')
p3
ggsave(plot = p3, device = 'png', filename = paste0('figures/TajimasD_change_vs_pos_NEA_CAN_runmean', width, '.angsd.gatk.png'), 
       width = 7.5, height = 6, units = 'in', dpi = 300)



#######################
# plots for p-values
#######################

# plot pi p-value vs. position (all loci)
# show decreases in pi below the x-axis
# cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
# p1 <- ggplot(dat, aes(POSgen, -log10(tPd.p)*sign(tPd), color = Chromo)) + 
#   geom_point(size = 0.2, alpha = 0.3) +
#   facet_wrap(~pop, ncol = 1) +
#   scale_color_manual(values = cols) +
#   geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
# p1
# 
# ggsave(plot = p1, device = 'png', filename = 'figures/pi_change.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)



# plot pi p-value vs. position (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(datgatk, aes(POSgen, -log10(tPd.p)*sign(tPd), color = Chromo)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p2

ggsave(plot = p2, device = 'png', filename = 'figures/pi_change.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)



# plot thetaW p-value vs. position (all loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p3 <- ggplot(dat, aes(POSgen, -log10(tWd.p)*sign(tWd), color = Chromo)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p3

ggsave(plot = p3, device = 'png', filename = 'figures/thetaW_change.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)

# plot thetaW p-value vs. position (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p4 <- ggplot(datgatk, aes(POSgen, -log10(tWd.p)*sign(tWd), color = Chromo)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p4

ggsave(plot = p4, device = 'png', filename = 'figures/thetaW_change.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)



# plot Tajama's D p-value vs. position (all loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p5 <- ggplot(dat, aes(POSgen, -log10(tDd.p)*sign(tDd), color = Chromo)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p5

ggsave(plot = p5, device = 'png', filename = 'figures/tajimasD_change.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)


# plot Tajama's D p-value vs. position (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p6 <- ggplot(datgatk, aes(POSgen, -log10(tDd.p)*sign(tDd), color = Chromo)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p6

ggsave(plot = p6, device = 'png', filename = 'figures/tajimasD_change.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)


#################
# print outliers
#################


datgatk[tPd.p < 0.05 & pop == 'can',]
datgatk[tPd.p < 0.05 & pop == 'lof0711',]
datgatk[tPd.p < 0.05 & pop == 'lof0714',]

datgatk[tWd.p < 0.05 & pop == 'can',]
datgatk[tWd.p < 0.05 & pop == 'lof0711',]
datgatk[tWd.p < 0.05 & pop == 'lof0714',]

datgatk[tDd.p < 0.05 & pop == 'can',]
datgatk[tDd.p < 0.05 & pop == 'lof0711',]
datgatk[tDd.p < 0.05 & pop == 'lof0714',]
