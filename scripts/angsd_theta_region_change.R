# calculate change in pi and Tajima's D through time
# for ANGSD-called loci
# run after angsd_theta.sh


######################
# calculate theta change
######################
require(data.table)

# read in data (all loci)
windsz <- '5e4'
datCan40 <- fread('analysis/thetas.windowed.Can_40.pestPG')
datCan14 <- fread('analysis/thetas.windowed.Can_14.pestPG')
dat07 <- fread('analysis/thetas.windowed.Lof_07.pestPG')
dat11 <- fread('analysis/thetas.windowed.Lof_11.pestPG')
dat14 <- fread('analysis/thetas.windowed.Lof_14.pestPG')

# merge
can <- merge(datCan40[, .(Chr, WinCenter, tW1 = tW, tP1 = tP, Tajima1 = Tajima)], 
             datCan14[, .(Chr, WinCenter, tW2 = tW, tP2 = tP, Tajima2 = Tajima)])
lof0711 <- merge(dat07[, .(Chr, WinCenter, tW1 = tW, tP1 = tP, Tajima1 = Tajima)], 
                 dat11[, .(Chr, WinCenter, tW2 = tW, tP2 = tP, Tajima2 = Tajima)])
lof0714 <- merge(dat07[, .(Chr, WinCenter, tW1 = tW, tP1 = tP, Tajima1 = Tajima)], 
                 dat14[, .(Chr, WinCenter, tW2 = tW, tP2 = tP, Tajima2 = Tajima)])
lof1114 <- merge(dat11[, .(Chr, WinCenter, tW1 = tW, tP1 = tP, Tajima1 = Tajima)], 
                 dat14[, .(Chr, WinCenter, tW2 = tW, tP2 = tP, Tajima2 = Tajima)])

# calculate change
can[, ':='(tWd = tW2 - tW1, tPd = tP2 - tP1, Tajimadiff = Tajima2 - Tajima1)]
lof0711[, ':='(tWd = tW2 - tW1, tPd = tP2 - tP1, Tajimadiff = Tajima2 - Tajima1)]
lof0714[, ':='(tWd = tW2 - tW1, tPd = tP2 - tP1, Tajimadiff = Tajima2 - Tajima1)]
lof1114[, ':='(tWd = tW2 - tW1, tPd = tP2 - tP1, Tajimadiff = Tajima2 - Tajima1)]

# merge populations together
can[, pop := 'can']
lof0711[, pop := 'lof0711']
lof0714[, pop := 'lof0714']
lof1114[, pop := 'lof1114']
bins <- rbind(can, lof0711, lof0714, lof1114)


# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- datCan40[, .(len = max(WinCenter) + 25000), by = Chr]
chrmax$start = c(0, cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, Chr)

setkey(bins, Chr)
bins <- merge(bins, chrmax[,.(Chr, start)])
bins[, POSgen := WinCenter + start]
bins[, start := NULL]



# save region change data
filenm <- paste('analysis/theta_change_region_', windsz, '.csv.gz', sep='')
filenm
write.csv(bins, file = gzfile(filenm))





###########################################
# plot theta change from regions
###########################################
require(data.table)
require(RColorBrewer)
require(ggplot2)

# read in
bins <- fread('analysis/theta_change_region_5e4.csv.gz', drop = 1); width='5e4'


# plot pi change
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p1 <- ggplot(bins, aes(POSgen, tPd/5e4, color = Chr)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols) +
  ylab('Change in pi per site')
p1
ggsave(plot = p1, device = 'png', filename = paste0('figures/pi_change_vs_pos_NEA_CAN_runmean', width, '.angsd.png'), 
       width = 7.5, height = 6, units = 'in', dpi = 300)

# plot theta_Waterson change
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(bins, aes(POSgen, tWd/5e4, color = Chr)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols) +
  ylab('Change in Wattersons theta per site')
p2
ggsave(plot = p2, device = 'png', filename = paste0('figures/thetaW_change_vs_pos_NEA_CAN_runmean', width, '.angsd.png'), 
       width = 7.5, height = 6, units = 'in', dpi = 300)

# plot Tajima's D change
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p3 <- ggplot(bins, aes(POSgen, Tajimadiff/5e4, color = Chr)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols) +
  ylab('Change in Tajimas D per site')
p3
ggsave(plot = p3, device = 'png', filename = paste0('figures/TajimasD_change_vs_pos_NEA_CAN_runmean', width, '.angsd.png'), 
       width = 7.5, height = 6, units = 'in', dpi = 300)

