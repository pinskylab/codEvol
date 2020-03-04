# Plot output from PCAngsd selection scan

library(RcppCNPy)
library(data.table)

# read in files
datcan <- as.data.table(npyLoad('analysis/pcangsd_can.selection.npy')) # selection statistics along a tested PC and they are χ²-distributed with 1 degree of freedom
datlof0711 <- as.data.table(npyLoad('analysis/pcangsd_lof0711.selection.npy'))
datlof0714 <- as.data.table(npyLoad('analysis/pcangsd_lof0714.selection.npy'))
datlof1114 <- as.data.table(npyLoad('analysis/pcangsd_lof1114.selection.npy'))

datcangt <- as.data.table(npyLoad('analysis/pcangsd_can_gatk.selection.npy')) # only gatk loci
datlof0711gt <- as.data.table(npyLoad('analysis/pcangsd_lof0711_gatk.selection.npy'))
datlof0714gt <- as.data.table(npyLoad('analysis/pcangsd_lof0714_gatk.selection.npy'))
datlof1114gt <- as.data.table(npyLoad('analysis/pcangsd_lof1114_gatk.selection.npy'))

datcangtnd <- as.data.table(npyLoad('analysis/pcangsd_can_gatk_no_dam.selection.npy')) # only gatk loci no aDNA damaga
datlof0711gtnd <- as.data.table(npyLoad('analysis/pcangsd_lof0711_gatk_no_dam.selection.npy'))
datlof0714gtnd <- as.data.table(npyLoad('analysis/pcangsd_lof0714_gatk_no_dam.selection.npy'))
datlof1114gtnd <- as.data.table(npyLoad('analysis/pcangsd_lof1114_gatk_no_dam.selection.npy'))

datcangtei <- as.data.table(npyLoad('analysis/pcangsd_can_gatk_ex_inv.selection.npy')) # only gatk loci no aDNA damaga
datlof0711gtei <- as.data.table(npyLoad('analysis/pcangsd_lof0711_gatk_ex_inv.selection.npy'))
datlof0714gtei <- as.data.table(npyLoad('analysis/pcangsd_lof0714_gatk_ex_inv.selection.npy'))
datlof1114gtei <- as.data.table(npyLoad('analysis/pcangsd_lof1114_gatk_ex_inv.selection.npy'))

sitescan <- fread('analysis/pcangsd_can.sites', header = FALSE) # sites to go with matrix
siteslof0711 <- fread('analysis/pcangsd_lof0711.sites', header = FALSE)
siteslof0714 <- fread('analysis/pcangsd_lof0714.sites', header = FALSE)
siteslof1114 <- fread('analysis/pcangsd_lof1114.sites', header = FALSE)

sitescangt <- fread('analysis/pcangsd_can_gatk.sites', header = FALSE) # sites to go with matrix
siteslof0711gt <- fread('analysis/pcangsd_lof0711_gatk.sites', header = FALSE)
siteslof0714gt <- fread('analysis/pcangsd_lof0714_gatk.sites', header = FALSE)
siteslof1114gt <- fread('analysis/pcangsd_lof1114_gatk.sites', header = FALSE)

sitescangtnd <- fread('analysis/pcangsd_can_gatk_no_dam.sites', header = FALSE) # sites to go with matrix
siteslof0711gtnd <- fread('analysis/pcangsd_lof0711_gatk_no_dam.sites', header = FALSE)
siteslof0714gtnd <- fread('analysis/pcangsd_lof0714_gatk_no_dam.sites', header = FALSE)
siteslof1114gtnd <- fread('analysis/pcangsd_lof1114_gatk_no_dam.sites', header = FALSE)

sitescangtei <- fread('analysis/pcangsd_can_gatk_ex_inv.sites', header = FALSE) # sites to go with matrix
siteslof0711gtei <- fread('analysis/pcangsd_lof0711_gatk_ex_inv.sites', header = FALSE)
siteslof0714gtei <- fread('analysis/pcangsd_lof0714_gatk_ex_inv.sites', header = FALSE)
siteslof1114gtei <- fread('analysis/pcangsd_lof1114_gatk_ex_inv.sites', header = FALSE)


gatk <- fread('data_31_01_20/GATK_filtered_SNP_set.tab') # snps that pass GATK filters
gatk2 <- fread('data_31_01_20/GATK_filtered_SNP_set_no_Dam.tab') # another snp set that pass filters


# label and combine files
datcan[, ':='(pop = 'can', type = 'all')]
datlof0711[, ':='(pop = 'lof0711', type = 'all')]
datlof0714[, ':='(pop = 'lof0714', type = 'all')]
datlof1114[, ':='(pop = 'lof1114', type = 'all')]

datcangt[, ':='(pop = 'can', type = 'gatk')]
datlof0711gt[, ':='(pop = 'lof0711', type = 'gatk')]
datlof0714gt[, ':='(pop = 'lof0714', type = 'gatk')]
datlof1114gt[, ':='(pop = 'lof1114', type = 'gatk')]

datcangtnd[, ':='(pop = 'can', type = 'gatk_nodam')]
datlof0711gtnd[, ':='(pop = 'lof0711', type = 'gatk_nodam')]
datlof0714gtnd[, ':='(pop = 'lof0714', type = 'gatk_nodam')]
datlof1114gtnd[, ':='(pop = 'lof1114', type = 'gatk_nodam')]

datcangtei[, ':='(pop = 'can', type = 'gatk_exinv')]
datlof0711gtei[, ':='(pop = 'lof0711', type = 'gatk_exinv')]
datlof0714gtei[, ':='(pop = 'lof0714', type = 'gatk_exinv')]
datlof1114gtei[, ':='(pop = 'lof1114', type = 'gatk_exinv')]

sitescan[, ':='(pop = 'can', type = 'all')]
siteslof0711[, ':='(pop = 'lof0711', type = 'all')]
siteslof0714[, ':='(pop = 'lof0714', type = 'all')]
siteslof1114[, ':='(pop = 'lof1114', type = 'all')]

sitescangt[, ':='(pop = 'can', type = 'gatk')]
siteslof0711gt[, ':='(pop = 'lof0711', type = 'gatk')]
siteslof0714gt[, ':='(pop = 'lof0714', type = 'gatk')]
siteslof1114gt[, ':='(pop = 'lof1114', type = 'gatk')]

sitescangtnd[, ':='(pop = 'can', type = 'gatk_nodam')]
siteslof0711gtnd[, ':='(pop = 'lof0711', type = 'gatk_nodam')]
siteslof0714gtnd[, ':='(pop = 'lof0714', type = 'gatk_nodam')]
siteslof1114gtnd[, ':='(pop = 'lof1114', type = 'gatk_nodam')]

sitescangtei[, ':='(pop = 'can', type = 'gatk_exinv')]
siteslof0711gtei[, ':='(pop = 'lof0711', type = 'gatk_exinv')]
siteslof0714gtei[, ':='(pop = 'lof0714', type = 'gatk_exinv')]
siteslof1114gtei[, ':='(pop = 'lof1114', type = 'gatk_exinv')]

dat <- rbind(datcan, datlof0711, datlof0714, datlof1114, 
             datcangt, datlof0711gt, datlof0714gt, datlof1114gt, 
             datcangtnd, datlof0711gtnd, datlof0714gtnd, datlof1114gtnd, 
             datcangtei, datlof0711gtei, datlof0714gtei, datlof1114gtei)

sites <- rbind(sitescan, siteslof0711, siteslof0714, siteslof1114, 
             sitescangt, siteslof0711gt, siteslof0714gt, siteslof1114gt, 
             sitescangtnd, siteslof0711gtnd, siteslof0714gtnd, siteslof1114gtnd, 
             sitescangtei, siteslof0711gtei, siteslof0714gtei, siteslof1114gtei)

rm(datcan, datlof0711, datlof0714, datlof1114, 
   datcangt, datlof0711gt, datlof0714gt, datlof1114gt, 
   datcangtnd, datlof0711gtnd, datlof0714gtnd, datlof1114gtnd, 
   datcangtei, datlof0711gtei, datlof0714gtei, datlof1114gtei)

rm(sitescan, siteslof0711, siteslof0714, siteslof1114, 
      sitescangt, siteslof0711gt, siteslof0714gt, siteslof1114gt, 
      sitescangtnd, siteslof0711gtnd, siteslof0714gtnd, siteslof1114gtnd, 
      sitescangtei, siteslof0711gtei, siteslof0714gtei, siteslof1114gtei)

# separate CHROM and POS
sites[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]
sites[, POS := as.integer(POS)]

# merge selection test with site names
nrow(dat)
nrow(sites)
dat <- cbind(sites[, .(CHROM, POS)], dat)

rm(sites)

# merge in snp sets
gatk[, gatk1 := 1] # 1 marks that a locus is in the snp set
gatk2[, gatk2 := 1]

nrow(dat)
dat <- merge(dat, gatk, by = c('CHROM', 'POS'), all.x = TRUE)
dat <- merge(dat, gatk2, by = c('CHROM', 'POS'), all.x = TRUE)
nrow(dat)

dat[is.na(gatk1), gatk1 := 0] # turn NA to 0
dat[is.na(gatk2), gatk2 := 0]

nrow(gatk) # 720987
dat[, sum(gatk1), by = .(pop, type)] # fewer than in snp set

nrow(gatk2) # 265784
dat[, sum(gatk2), by = .(pop, type)] # slightly fewer

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- dat[, .(len = max(POS)), by = CHROM]
chrmax$start=c(0, cumsum(chrmax$len)[1:(nrow(chrmax) - 1)])
setkey(chrmax, CHROM)
setkey(dat, CHROM)

dat <- merge(dat, chrmax[,.(CHROM, start)], by='CHROM')

dat[,POSgen := POS + start]
dat[,start := NULL]


# calculate p-value
dat[, p := pchisq(q = V1, df = 1, lower.tail = FALSE)]

# FDR correct by pop and type
dat[, pfdr := p.adjust(p, method = 'fdr'), by = .(pop, type)]


# add a vector for color by LG
cols <- c('#a6cee3aa', '#1f78b4aa') # light blue, blue, partially transparent: for alternating LGs
dat[, lgcol := cols[1]]
dat[CHROM %in% chrmax$CHROM[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]

# plot and save out
# all loci
ylims <- c(0, max(-log10(c(datcan$pfdr, datlof0711$pfdr, datlof0714$pfdr, datlof1114$pfdr))))

png(filename = paste0('figures/pcangsd_selscan.png'), width = 20, height = 12, units = "in", res = 150)
par(mfrow = c(4,1))

dat[pop == 'can' & type == 'all', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)", 
              main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0711' & type == 'all', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                  main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0714' & type == 'all', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                  main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof1114' & type == 'all', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                  main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()


# first set of GATK-filtered loci
ylims <- c(0, max(-log10(dat[type == 'gatk', pfdr])))

png(filename = paste0('figures/pcangsd_selscan_gatk.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))

dat[pop == 'can' & type == 'gatk', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)", 
                        main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0711' & type == 'gatk', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0714' & type == 'gatk', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof1114' & type == 'gatk', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()


# GATK-filtered loci no aDNA damage
ylims <- c(0, max(-log10(dat[type == 'gatk_nodam', pfdr])))

png(filename = paste0('figures/pcangsd_selscan_gatkNoDam.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))

dat[pop == 'can' & type == 'gatk_nodam', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)", 
                                        main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0711' & type == 'gatk_nodam', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                                            main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0714' & type == 'gatk_nodam', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                                            main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof1114' & type == 'gatk_nodam', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                                            main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()

# GATK-filtered loci no inversion LGs
ylims <- c(0, max(-log10(dat[type == 'gatk_exinv', pfdr])))

png(filename = paste0('figures/pcangsd_selscan_gatkExInv.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))

dat[pop == 'can' & type == 'gatk_exinv', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)", 
                                              main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0711' & type == 'gatk_exinv', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                                                  main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof0714' & type == 'gatk_exinv', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                                                  main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dat[pop == 'lof1114' & type == 'gatk_exinv', plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                                                  main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()



#####################
##Look at the PCAs
pcacan <- fread('analysis/pcangsd_can.cov') # pca covariance matrix. cols are individuals, rows are PCs
pcalof0711 <- fread('analysis/pcangsd_lof0711.cov')
pcalof0714 <- fread('analysis/pcangsd_lof0714.cov')
pcalof1114 <- fread('analysis/pcangsd_lof1114.cov')

pcacangt <- fread('analysis/pcangsd_can_gatk.cov') # only GATK snps
pcalof0711gt <- fread('analysis/pcangsd_lof0711_gatk.cov')
pcalof0714gt <- fread('analysis/pcangsd_lof0714_gatk.cov')
pcalof1114gt <- fread('analysis/pcangsd_lof1114_gatk.cov')

pcacangtnd <- fread('analysis/pcangsd_can_gatk_no_dam.cov') # only GATK no aDNA damage snps
pcalof0711gtnd <- fread('analysis/pcangsd_lof0711_gatk_no_dam.cov')
pcalof0714gtnd <- fread('analysis/pcangsd_lof0714_gatk_no_dam.cov')
pcalof1114gtnd <- fread('analysis/pcangsd_lof1114_gatk_no_dam.cov')

pcacangtei <- fread('analysis/pcangsd_can_gatk_ex_inv.cov') # only GATK w/out inversions
pcalof0711gtei <- fread('analysis/pcangsd_lof0711_gatk_ex_inv.cov')
pcalof0714gtei <- fread('analysis/pcangsd_lof0714_gatk_ex_inv.cov')
pcalof1114gtei <- fread('analysis/pcangsd_lof1114_gatk_ex_inv.cov')

# eigenvector decomposition of the covariance matrix
ecan <- eigen(pcacan)
elof0711 <- eigen(pcalof0711)
elof0714 <- eigen(pcalof0714)
elof1114 <- eigen(pcalof1114)

ecangt <- eigen(pcacangt)
elof0711gt <- eigen(pcalof0711gt)
elof0714gt <- eigen(pcalof0714gt)
elof1114gt <- eigen(pcalof1114gt)

ecangtnd <- eigen(pcacangtnd)
elof0711gtnd <- eigen(pcalof0711gtnd)
elof0714gtnd <- eigen(pcalof0714gtnd)
elof1114gtnd <- eigen(pcalof1114gtnd)

ecangtei <- eigen(pcacangtei)
elof0711gtei <- eigen(pcalof0711gtei)
elof0714gtei <- eigen(pcalof0714gtei)
elof1114gtei <- eigen(pcalof1114gtei)

# plot PCAs
# red is early, black is late
pdf('figures/pcangsd_pca.pdf', width = 6, height = 6)
par(mfrow=c(2,2), mai = c(0.6, 0.7, 0.3, 0.1), las = 1, mgp = c(2, 0.7, 0))
# All loci
plot(ecan$vectors[, 1:2], col = rep(c('red', 'black'), c(21,24)), xlab = 'PC1', ylab = 'PC2',
     main = 'Canada all loci')

plot(elof0711$vectors[, 1:2], col = rep(c('red', 'black'), c(22,23)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-11 all loci')

plot(elof0714$vectors[, 1:2], col = rep(c('red', 'black'), c(22,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-14 all loci')

plot(elof1114$vectors[, 1:2], col = rep(c('red', 'black'), c(23,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 11-14 all loci')

legend('top', pch = 1, col = c('red', 'black'), legend = c('early', 'late'))



# GATK loci
plot(ecangt$vectors[, 1:2], col = rep(c('red', 'black'), c(21,24)), xlab = 'PC1', ylab = 'PC2',
     main = 'Canada GATK')

plot(elof0711gt$vectors[, 1:2], col = rep(c('red', 'black'), c(22,23)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-11 GATK')

plot(elof0714gt$vectors[, 1:2], col = rep(c('red', 'black'), c(22,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-14 GATK')

plot(elof1114gt$vectors[, 1:2], col = rep(c('red', 'black'), c(23,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 11-14 GATK')


legend('top', pch = 1, col = c('red', 'black'), legend = c('early', 'late'))

# GATK no damage
plot(ecangtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(21,24)), xlab = 'PC1', ylab = 'PC2',
     main = 'Canada GATK no dam')

plot(elof0711gtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(22,23)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-11 GATK no dam')

plot(elof0714gtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(22,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-14 GATK no dam')

plot(elof1114gtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(23,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 11-14 GATK no dam')


# GATK no inv
plot(ecangtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(21,24)), xlab = 'PC1', ylab = 'PC2',
     main = 'Canada GATK ex inv')

plot(elof0711gtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(22,23)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-11 GATK ex inv')

plot(elof0714gtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(22,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-14 GATK ex inv')

plot(elof1114gtnd$vectors[, 1:2], col = rep(c('red', 'black'), c(23,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 11-14 GATK ex inv')


dev.off()