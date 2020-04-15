## examine results from site reshuffled theta from ANGSD genotypes
## run after angsd_theta_siteshuffle_null.sh/.r

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

# max theta per genome from reshuffling (all sites)
nullcan <- fread('analysis/theta.siteshuffle.Can_40.Can_14.csv.gz')
nulllof0711 <- fread('analysis/theta.siteshuffle.Lof_07.Lof_11.csv.gz')
nulllof0714 <- fread('analysis/theta.siteshuffle.Lof_07.Lof_14.csv.gz')
nulllof1114 <- fread('analysis/theta.siteshuffle.Lof_11.Lof_14.csv.gz')

# max theta per genome from reshuffling (GATK sites)
nullcangatk <- fread('analysis/theta.siteshuffle.Can_40.Can_14.gatk.csv.gz')
nulllof0711gatk <- fread('analysis/theta.siteshuffle.Lof_07.Lof_11.gatk.csv.gz')
nulllof0714gatk <- fread('analysis/theta.siteshuffle.Lof_07.Lof_14.gatk.csv.gz')
nulllof1114gatk <- fread('analysis/theta.siteshuffle.Lof_11.Lof_14.gatk.csv.gz')


# sliding window theta change
dat <- fread('analysis/theta_change_region_5e4.csv.gz', drop = 1) 

# sliding windows theta change (GATK sites)
datgatk <- fread('analysis/theta_change_region_5e4.gatk.csv.gz', drop = 1) 


# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- fread('data/lg_length.csv')
chrmax$start=c(0,cumsum(chrmax$length)[1:(nrow(chrmax)-1)])
setkey(chrmax, chr)

# merge nucleotide position into the frequency files
setkey(dat, Chr)
dat <- dat[chrmax[, .(Chr = chr, start)], ]
dat[, POSgen := WinCenter + start]
dat[,start := NULL]

setkey(datgatk, Chr)
datgatk <- datgatk[chrmax[, .(Chr = chr, start)], ]
datgatk[, POSgen := WinCenter + start]
datgatk[,start := NULL]


#######################
## null model stats
#######################

nullcan[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95), # using 5% and 95% to be lenient
            tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
            tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
nulllof0711[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
            tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
            tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
nulllof0714[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
                tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
                tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]
nulllof1114[, .(tWd_l95 = quantile(mintWd, 0.05), tWd_u95 = quantile(maxtWd, probs = 0.95),
                tPd_l95 = quantile(mintPd, 0.05), tPd_u95 = quantile(maxtPd, probs = 0.95),
                tDd_l95 = quantile(mintDd, 0.05), tDd_u95 = quantile(maxtDd, probs = 0.95))]

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

dat[pop == 'can' & tWd > 0, tWd.p := calcpG(tWd, nullcan$maxtWd), by = .(Chr, WinCenter)] # thetaW all loci
dat[pop == 'can' & tWd <= 0, tWd.p := calcpL(tWd, nullcan$mintWd), by = .(Chr, WinCenter)]
dat[pop == 'lof0711' & tWd > 0, tWd.p := calcpG(tWd, nulllof0711$maxtWd), by = .(Chr, WinCenter)]
dat[pop == 'lof0711' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0711$mintWd), by = .(Chr, WinCenter)]
dat[pop == 'lof0714' & tWd > 0, tWd.p := calcpG(tWd, nulllof0714$maxtWd), by = .(Chr, WinCenter)]
dat[pop == 'lof0714' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0714$mintWd), by = .(Chr, WinCenter)]
dat[pop == 'lof1114' & tWd > 0, tWd.p := calcpG(tWd, nulllof1114$maxtWd), by = .(Chr, WinCenter)]
dat[pop == 'lof1114' & tWd <= 0, tWd.p := calcpL(tWd, nulllof1114$mintWd), by = .(Chr, WinCenter)]

dat[pop == 'can' & tPd > 0, tPd.p := calcpG(tPd, nullcan$maxtPd), by = .(Chr, WinCenter)] # theta pi
dat[pop == 'can' & tPd <= 0, tPd.p := calcpL(tPd, nullcan$mintPd), by = .(Chr, WinCenter)]
dat[pop == 'lof0711' & tPd > 0, tPd.p := calcpG(tPd, nulllof0711$maxtPd), by = .(Chr, WinCenter)]
dat[pop == 'lof0711' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0711$mintPd), by = .(Chr, WinCenter)]
dat[pop == 'lof0714' & tPd > 0, tPd.p := calcpG(tPd, nulllof0714$maxtPd), by = .(Chr, WinCenter)]
dat[pop == 'lof0714' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0714$mintPd), by = .(Chr, WinCenter)]
dat[pop == 'lof1114' & tPd > 0, tPd.p := calcpG(tPd, nulllof1114$maxtPd), by = .(Chr, WinCenter)]
dat[pop == 'lof1114' & tPd <= 0, tPd.p := calcpL(tPd, nulllof1114$mintPd), by = .(Chr, WinCenter)]

dat[pop == 'can' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nullcan$maxtDd), by = .(Chr, WinCenter)] # tajima's D
dat[pop == 'can' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nullcan$mintDd), by = .(Chr, WinCenter)]
dat[pop == 'lof0711' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof0711$maxtDd), by = .(Chr, WinCenter)]
dat[pop == 'lof0711' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof0711$mintDd), by = .(Chr, WinCenter)]
dat[pop == 'lof0714' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof0714$maxtDd), by = .(Chr, WinCenter)]
dat[pop == 'lof0714' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof0714$mintDd), by = .(Chr, WinCenter)]
dat[pop == 'lof1114' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof1114$maxtDd), by = .(Chr, WinCenter)]
dat[pop == 'lof1114' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof1114$mintDd), by = .(Chr, WinCenter)]


datgatk[pop == 'can' & tWd > 0, tWd.p := calcpG(tWd, nullcangatk$maxtWd), by = .(Chr, WinCenter)] # thetaW gatk loci
datgatk[pop == 'can' & tWd <= 0, tWd.p := calcpL(tWd, nullcangatk$mintWd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0711' & tWd > 0, tWd.p := calcpG(tWd, nulllof0711gatk$maxtWd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0711' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0711gatk$mintWd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0714' & tWd > 0, tWd.p := calcpG(tWd, nulllof0714gatk$maxtWd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0714' & tWd <= 0, tWd.p := calcpL(tWd, nulllof0714gatk$mintWd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof1114' & tWd > 0, tWd.p := calcpG(tWd, nulllof1114gatk$maxtWd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof1114' & tWd <= 0, tWd.p := calcpL(tWd, nulllof1114gatk$mintWd), by = .(Chr, WinCenter)]

datgatk[pop == 'can' & tPd > 0, tPd.p := calcpG(tPd, nullcangatk$maxtPd), by = .(Chr, WinCenter)] # theta pi
datgatk[pop == 'can' & tPd <= 0, tPd.p := calcpL(tPd, nullcangatk$mintPd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0711' & tPd > 0, tPd.p := calcpG(tPd, nulllof0711gatk$maxtPd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0711' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0711gatk$mintPd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0714' & tPd > 0, tPd.p := calcpG(tPd, nulllof0714gatk$maxtPd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0714' & tPd <= 0, tPd.p := calcpL(tPd, nulllof0714gatk$mintPd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof1114' & tPd > 0, tPd.p := calcpG(tPd, nulllof1114gatk$maxtPd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof1114' & tPd <= 0, tPd.p := calcpL(tPd, nulllof1114gatk$mintPd), by = .(Chr, WinCenter)]

datgatk[pop == 'can' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nullcangatk$maxtDd), by = .(Chr, WinCenter)] # tajima's D
datgatk[pop == 'can' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nullcangatk$mintDd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0711' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof0711gatk$maxtDd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0711' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof0711gatk$mintDd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0714' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof0714gatk$maxtDd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof0714' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof0714gatk$mintDd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof1114' & Tajimadiff > 0, tDd.p := calcpG(Tajimadiff, nulllof1114gatk$maxtDd), by = .(Chr, WinCenter)]
datgatk[pop == 'lof1114' & Tajimadiff <= 0, tDd.p := calcpL(Tajimadiff, nulllof1114gatk$mintDd), by = .(Chr, WinCenter)]



##############
# plots
##############

# plot pi p-value vs. position (all loci)
# show decreases in pi below the x-axis
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p1 <- ggplot(dat, aes(POSgen, -log10(tPd.p)*sign(tPd), color = Chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p1

ggsave(plot = p1, device = 'png', filename = 'figures/pi_change.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)



# plot pi p-value vs. position (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(datgatk, aes(POSgen, -log10(tPd.p)*sign(tPd), color = Chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p2

ggsave(plot = p2, device = 'png', filename = 'figures/pi_change.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)



# plot thetaW p-value vs. position (all loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p3 <- ggplot(dat, aes(POSgen, -log10(tWd.p)*sign(tWd), color = Chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p3

ggsave(plot = p3, device = 'png', filename = 'figures/thetaW_change.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)

# plot thetaW p-value vs. position (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p4 <- ggplot(datgatk, aes(POSgen, -log10(tWd.p)*sign(tWd), color = Chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p4

ggsave(plot = p4, device = 'png', filename = 'figures/thetaW_change.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)



# plot Tajama's D p-value vs. position (all loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p5 <- ggplot(dat, aes(POSgen, -log10(tDd.p)*sign(Tajimadiff), color = Chr)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p5

ggsave(plot = p5, device = 'png', filename = 'figures/tajimasD_change.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)


# plot Tajama's D p-value vs. position (gatk loci)
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p6 <- ggplot(datgatk, aes(POSgen, -log10(tDd.p)*sign(Tajimadiff), color = Chr)) + 
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


datgatk[tWd.p < 0.05 & pop == 'can',]
datgatk[tWd.p < 0.05 & pop == 'lof0711',]
datgatk[tWd.p < 0.05 & pop == 'lof0714',]
