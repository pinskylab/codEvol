# plot pFst output

require(ggplot2)
require(data.table)

# read in files
datcan <-fread('gunzip -c analysis/pfst_can.out.gz', header=FALSE, col.names = c('CHROM', 'POS', 'p'))
datlof0711 <-fread('gunzip -c analysis/pfst_lof0711.out.gz', header=FALSE, col.names = c('CHROM', 'POS', 'p'))
datlof0714 <-fread('gunzip -c analysis/pfst_lof0714.out.gz', header=FALSE, col.names = c('CHROM', 'POS', 'p'))
datlof1114 <-fread('gunzip -c analysis/pfst_lof1114.out.gz', header=FALSE, col.names = c('CHROM', 'POS', 'p'))
gatk <- fread('data_31_01_20/GATK_filtered_SNP_set.tab') # snps that pass GATK filters
gatk2 <- fread('data_31_01_20/GATK_filtered_SNP_set_no_Dam.tab') # another snp set that pass filters

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
datcan[, sum(gatk1)] # 720542 (slightly fewer than in snp set)
datlof0711[, sum(gatk1)] # 717978 (slightly fewer than in snp set)
datlof0714[, sum(gatk1)] # 717978 (slightly fewer than in snp set)
datlof1114[, sum(gatk1)] # 717977 (slightly fewer than in snp set)

nrow(gatk2) # 265784
datcan[, sum(gatk2)] # 265603 (slightly fewer)
datlof0711[, sum(gatk2)] # 264474 (slightly fewer than in snp set)
datlof0714[, sum(gatk2)] # 264474 (slightly fewer than in snp set)
datlof1114[, sum(gatk2)] # 264474 (slightly fewer than in snp set)


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

# fix p = 0 to min p
minp <- min(c(datcan[p > 0, min(p)], datlof0711[p > 0, min(p)], datlof0714[p > 0, min(p)], datlof1114[p > 0, min(p)]))
datcan[p == 0, p := minp]
datlof0711[p == 0, p := minp]
datlof0714[p == 0, p := minp]
datlof1114[p == 0, p := minp]


# fdr-adjust the p-values
datcan[, pfdr := p.adjust(p, method = 'fdr')]
datlof0711[, pfdr := p.adjust(p, method = 'fdr')]
datlof0714[, pfdr := p.adjust(p, method = 'fdr')]
datlof1114[, pfdr := p.adjust(p, method = 'fdr')]

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

png(filename = paste0('figures/pFst.png'), width = 20, height = 12, units = "in", res = 150)
par(mfrow = c(4,1))

datcan[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)", 
                        main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]

datlof0711[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                            main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]

datlof0714[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                            main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]

datlof1114[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                            main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]

dev.off()

# first set of GATK-filtered loci
ylims <- c(0, max(-log10(c(datcan$pfdr, datlof0711$pfdr, datlof0714$pfdr, datlof1114$pfdr))))

png(filename = paste0('figures/pFst_gatk.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))

datcan[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)", 
                        main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]

datlof0711[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                        main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]

datlof0714[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                           main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]

datlof1114[gatk1 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                           main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]

dev.off()

# second set of GATK-filtered loci
ylims <- c(0, max(-log10(c(datcan$pfdr, datlof0711$pfdr, datlof0714$pfdr, datlof1114$pfdr))))

png(filename = paste0('figures/pFst_gatkNoDam.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))

datcan[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)", 
                        main = 'Canada', col = lgcol, cex = 0.5, ylim = ylims)]

datlof0711[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                            main = 'Lof 07-11', col = lgcol, cex = 0.5, ylim = ylims)]

datlof0714[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                            main = 'Lof 07-14', col = lgcol, cex = 0.5, ylim = ylims)]

datlof1114[gatk2 == 1, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
                            main = 'Lof 11-14', col = lgcol, cex = 0.5, ylim = ylims)]

dev.off()
