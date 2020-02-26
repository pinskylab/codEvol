# Plot output from PCAngsd selection scan

library(RcppCNPy)
library(data.table)

# read in files
datcan <- as.data.table(npyLoad('analysis/pcangsd_can.selection.npy')) # selection statistics along a tested PC and they are χ²-distributed with 1 degree of freedom
datlof0711 <- as.data.table(npyLoad('analysis/pcangsd_lof0711.selection.npy'))
datlof0714 <- as.data.table(npyLoad('analysis/pcangsd_lof0714.selection.npy'))
datlof1114 <- as.data.table(npyLoad('analysis/pcangsd_lof1114.selection.npy'))

sitescan <- fread('analysis/pcangsd_can.sites', header = FALSE) # sites to go with matrix
siteslof0711 <- fread('analysis/pcangsd_lof0711.sites', header = FALSE)
siteslof0714 <- fread('analysis/pcangsd_lof0714.sites', header = FALSE)
siteslof1114 <- fread('analysis/pcangsd_lof1114.sites', header = FALSE)

gatk <- fread('data_31_01_20/GATK_filtered_SNP_set.tab') # snps that pass GATK filters
gatk2 <- fread('data_31_01_20/GATK_filtered_SNP_set_no_Dam.tab') # another snp set that pass filters


# separate CHROM and POS
sitescan[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]
siteslof0711[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]
siteslof0714[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]
siteslof1114[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]

sitescan[, POS := as.integer(POS)]
siteslof0711[, POS := as.integer(POS)]
siteslof0714[, POS := as.integer(POS)]
siteslof1114[, POS := as.integer(POS)]

# calculate p-value
datcan[, p := pchisq(q = V1, df = 1, lower.tail = FALSE)]
datlof0711[, p := pchisq(q = V1, df = 1, lower.tail = FALSE)]
datlof0714[, p := pchisq(q = V1, df = 1, lower.tail = FALSE)]
datlof1114[, p := pchisq(q = V1, df = 1, lower.tail = FALSE)]

# FDR correct
datcan[, pfdr := p.adjust(p, method = 'fdr')]
datlof0711[, pfdr := p.adjust(p, method = 'fdr')]
datlof0714[, pfdr := p.adjust(p, method = 'fdr')]
datlof1114[, pfdr := p.adjust(p, method = 'fdr')]

# merge selection test with site names
datcan <- cbind(sitescan[, .(CHROM, POS)], datcan)
datlof0711 <- cbind(siteslof0711[, .(CHROM, POS)], datlof0711)
datlof0714 <- cbind(siteslof0714[, .(CHROM, POS)], datlof0714)
datlof1114 <- cbind(siteslof1114[, .(CHROM, POS)], datlof1114)

# merge in snp sets
gatk[, gatk1 := 1] # 1 marks that a locus is in the snp set
gatk2[, gatk2 := 1]

nrow(datcan)
datcan <- merge(datcan, gatk, by = c('CHROM', 'POS'), all.x = TRUE)
datcan <- merge(datcan, gatk2, by = c('CHROM', 'POS'), all.x = TRUE)
nrow(datcan)

nrow(datlof0711)
datlof0711 <- merge(datlof0711, gatk, by = c('CHROM', 'POS'), all.x = TRUE)
datlof0711 <- merge(datlof0711, gatk2, by = c('CHROM', 'POS'), all.x = TRUE)
nrow(datlof0711)

nrow(datlof0714)
datlof0714 <- merge(datlof0714, gatk, by = c('CHROM', 'POS'), all.x = TRUE)
datlof0714 <- merge(datlof0714, gatk2, by = c('CHROM', 'POS'), all.x = TRUE)
nrow(datlof0714)

nrow(datlof1114)
datlof1114 <- merge(datlof1114, gatk, by = c('CHROM', 'POS'), all.x = TRUE)
datlof1114 <- merge(datlof1114, gatk2, by = c('CHROM', 'POS'), all.x = TRUE)
nrow(datlof1114)

datcan[is.na(gatk1), gatk1 := 0] # turn NA to 0
datcan[is.na(gatk2), gatk2 := 0]
datlof0711[is.na(gatk1), gatk1 := 0]
datlof0711[is.na(gatk2), gatk2 := 0]
datlof0714[is.na(gatk1), gatk1 := 0]
datlof0714[is.na(gatk2), gatk2 := 0]
datlof1114[is.na(gatk1), gatk1 := 0]
datlof1114[is.na(gatk2), gatk2 := 0]

nrow(gatk) # 720987
datcan[, sum(gatk1)] # 542865 (fewer than in snp set)
datlof0711[, sum(gatk1)] # 561950 (slightly fewer than in snp set)
datlof0714[, sum(gatk1)] # 556264 (slightly fewer than in snp set)
datlof1114[, sum(gatk1)] # 565334 (slightly fewer than in snp set)

nrow(gatk2) # 265784
datcan[, sum(gatk2)] # 200745 (slightly fewer)
datlof0711[, sum(gatk2)] # 210453 (slightly fewer than in snp set)
datlof0714[, sum(gatk2)] # 208494 (slightly fewer than in snp set)
datlof1114[, sum(gatk2)] # 211484 (slightly fewer than in snp set)


# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- datcan[, .(len = max(POS)), by = CHROM]
chrmax$start=c(0, cumsum(chrmax$len)[1:(nrow(chrmax) - 1)])
setkey(chrmax, CHROM)

setkey(datcan, CHROM)
setkey(datlof0711, CHROM)
setkey(datlof0714, CHROM)
setkey(datlof1114, CHROM)

datcan <- merge(datcan, chrmax[,.(CHROM, start)], by='CHROM')
datlof0711 <- merge(datlof0711, chrmax[,.(CHROM, start)], by='CHROM')
datlof0714 <- merge(datlof0714, chrmax[,.(CHROM, start)], by='CHROM')
datlof1114 <- merge(datlof1114, chrmax[,.(CHROM, start)], by='CHROM')

datcan[,POSgen := POS + start]
datlof0711[,POSgen := POS + start]
datlof0714[,POSgen := POS + start]
datlof1114[,POSgen := POS + start]

datcan[,start := NULL]
datlof0711[,start := NULL]
datlof0714[,start := NULL]
datlof1114[,start := NULL]


# add a vector for color by LG
cols <- c('#a6cee3aa', '#1f78b4aa') # light blue, blue, partially transparent: for alternating LGs

datcan[, lgcol := cols[1]]
datlof0711[, lgcol := cols[1]]
datlof0714[, lgcol := cols[1]]
datlof1114[, lgcol := cols[1]]

datcan[CHROM %in% chrmax$CHROM[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]
datlof0711[CHROM %in% chrmax$CHROM[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]
datlof0714[CHROM %in% chrmax$CHROM[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]
datlof1114[CHROM %in% chrmax$CHROM[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]

# plot and save out
# all loci
ylims <- c(0, max(-log10(c(datcan$pfdr, datlof0711$pfdr, datlof0714$pfdr, datlof1114$pfdr))))

png(filename = paste0('figures/pcangsd_selscan.png'), width = 20, height = 12, units = "in", res = 150)
par(mfrow = c(4,1))

datcan[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)", 
              main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof0711[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                  main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof0714[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                  main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof1114[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                  main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()


# first set of GATK-filtered loci
ylims <- c(0, max(-log10(c(datcan[gatk1 == 1, pfdr], datlof0711[gatk1 == 1, pfdr], 
                           datlof0714[gatk1 == 1, pfdr], datlof1114[gatk1 == 1, pfdr]))))

png(filename = paste0('figures/pcangsd_selscan_gatk.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))

datcan[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)", 
                        main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof0711[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof0714[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof1114[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()


# second set of GATK-filtered loci
ylims <- c(0, max(-log10(c(datcan[gatk2 == 1, pfdr], datlof0711[gatk2 == 1, pfdr], 
                           datlof0714[gatk2 == 1, pfdr], datlof1114[gatk2 == 1, pfdr]))))

png(filename = paste0('figures/pcangsd_selscan_gatkNoDam.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))

datcan[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)", 
                        main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof0711[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof0714[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

datlof1114[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(PCANGSD p, FDR-adjusted)",
                            main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()



#####################
##Look at the PCAs
pcacan <- fread('analysis/pcangsd_can.cov') # pca covariance matrix. cols are individuals, rows are PCs
pcalof0711 <- fread('analysis/pcangsd_lof0711.cov')
pcalof0714 <- fread('analysis/pcangsd_lof0714.cov')
pcalof1114 <- fread('analysis/pcangsd_lof1114.cov')

# eigenvector decomposition of the covariance matrix
ecan <- eigen(pcacan)
elof0711 <- eigen(pcalof0711)
elof0714 <- eigen(pcalof0714)
elof1114 <- eigen(pcalof1114)

# plot PCAs
# red is early, black is late
pdf('figures/pcangsd_pca.pdf', width = 6, height = 6, units = "in")
par(mfrow=c(2,2), mai = c(0.6, 0.7, 0.3, 0.1), las = 1, mgp = c(2, 0.7, 0))
plot(ecan$vectors[, 1:2], col = rep(c('red', 'black'), c(21,24)), xlab = 'PC1', ylab = 'PC2',
     main = 'Canada')

plot(elof0711$vectors[, 1:2], col = rep(c('red', 'black'), c(22,23)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-11')

plot(elof0714$vectors[, 1:2], col = rep(c('red', 'black'), c(22,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 07-14')

plot(elof1114$vectors[, 1:2], col = rep(c('red', 'black'), c(23,21)), xlab = 'PC1', ylab = 'PC2',
     main = 'Lofoten 11-14')


legend('top', pch = 1, col = c('red', 'black'), legend = c('early', 'late'))

dev.off()