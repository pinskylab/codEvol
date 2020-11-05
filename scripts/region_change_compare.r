# Compare freq change, LD change, Tajima's D change, pi change
# set up for 50k windows

require(data.table)
require(vioplot)

################################
# read in data
################################

# fst: sliding window fst from ANGSD (GATK nodam2 sites)
dat1 <- fread('analysis/Can_40.Can_14.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1) # output by angsd_fst.sh. skip headers.
dat2 <- fread('analysis/Lof_07.Lof_11.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1)
dat3 <- fread('analysis/Lof_07.Lof_14.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1)
dat4 <- fread('analysis/Lof_11.Lof_14.gatk.slide', col.names = c('region', 'CHROM', 'midPos', 'nloci', 'fst'), skip = 1)

dat1[, pop := 'can']
dat2[, pop := 'lof0711']
dat3[, pop := 'lof0714']
dat4[, pop := 'lof1114']

fst <- rbind(dat1, dat2, dat3, dat4)
fst[, region := NULL]
setnames(fst, c('midPos', 'CHROM'), c('WinCenter', 'chr'))
fst[, WinCenter := WinCenter - 1] # adjust to match the others

# fst: sliding window fst from ANGSD (GATK nodam2 unlinked sites)
fst.unl <- fread('output/fst_siteshuffle.angsd.gatk.csv.gz') # output by angsd_fst_siteshuffle_null_stats.r
setnames(fst.unl, c('CHROM', 'midPos', 'fst', 'nloci', 'p'), c('chr', 'WinCenter', 'fst.unl', 'nloci.unl', 'fst.p'))

# change in ld
ldwide <- fread('analysis/ld_change_region_5e4_ngsLD.gatk.csv.gz') # from ngsLD_region_change.r
ld <- rbind(ldwide[, .(pop = 'can', chr, WinCenter = pos1mid1, ld_diff = ld_diff_Can)], # reformat to wide
            ldwide[, .(pop = 'lof0711', chr, WinCenter = pos1mid1, ld_diff = ld_diff_Lof0711)],
            ldwide[, .(pop = 'lof0714', chr, WinCenter = pos1mid1, ld_diff = ld_diff_Lof0714)],
            ldwide[, .(pop = 'lof1114', chr, WinCenter = pos1mid1, ld_diff = ld_diff_Lof1114)])

# change in pi and D: sliding window pi from ANGSD (GATK nodam2 unlinked sites)
dat1 <- fread('analysis/theta_change_region_50000.Can_40.Can_14.gatk.unlinked.csv.gz', drop = 1) # from angsd_theta_siteshuffle_null.r
dat2 <- fread('analysis/theta_change_region_50000.Lof_07.Lof_11.gatk.unlinked.csv.gz', drop = 1)
dat3 <- fread('analysis/theta_change_region_50000.Lof_07.Lof_14.gatk.unlinked.csv.gz', drop = 1)
dat4 <- fread('analysis/theta_change_region_50000.Lof_11.Lof_14.gatk.unlinked.csv.gz', drop = 1)

dat1[, pop := 'can']
dat2[, pop := 'lof0711']
dat3[, pop := 'lof0714']
dat4[, pop := 'lof1114']

piD <- rbind(dat1, dat2, dat3, dat4)
piD[, ':='(tW1 = NULL, tP1 = NULL, tD1 = NULL, tW2 = NULL, tP2 = NULL, tD2 = NULL, tWd = NULL)]
setnames(piD, c('Chromo'), c('chr'))


# pi and Tajima's D p-values
piD.p <- fread('analysis/theta_siteshuffle.angsd.gatk.csv.gz') # output by angsd_theta_siteshuffle_null_stats.r
piD.p <- piD.p[, .(chr = Chromo, WinCenter, tPd.p, tDd.p, pop)]


# merge
bins <- merge(fst, fst.unl, by=c('chr', 'WinCenter', 'pop'), all=TRUE)
bins <- merge(bins, ld, by=c('chr', 'WinCenter', 'pop'), all=TRUE)
bins <- merge(bins, piD, by=c('chr', 'WinCenter', 'pop'), all=TRUE)
bins <- merge(bins, piD.p, by=c('chr', 'WinCenter', 'pop'), all=TRUE)

nrow(bins)
nrow(fst)
nrow(fst.unl)
nrow(ld)
nrow(piD)
nrow(piD.p)

###################################################################
## prep across populations
## trim to non-overlapping windows, reshape, and calc percentiles
###################################################################
# fst and fst unlinked
binswidefst <- dcast(bins[(WinCenter %% 25000) == 0, .(chr, WinCenter, pop, nloci, nloci.unl, fst, fst.unl, fst.p)], chr + WinCenter ~ pop, value.var = c('fst', 'fst.unl', 'fst.p', 'nloci', 'nloci.unl'))
nrow(binswidefst)
binswidefst[, fstperc.can := ecdf(fst_can)(fst_can)]
binswidefst[, fstperc.lof0711 := ecdf(fst_lof0711)(fst_lof0711)]
binswidefst[, fstperc.lof0714 := ecdf(fst_lof0714)(fst_lof0714)]
binswidefst[, fstperc.lof1114 := ecdf(fst_lof1114)(fst_lof1114)]

binswidefst[, fstunlperc.can := ecdf(fst.unl_can)(fst.unl_can)]
binswidefst[, fstunlperc.lof0711 := ecdf(fst.unl_lof0711)(fst.unl_lof0711)]
binswidefst[, fstunlperc.lof0714 := ecdf(fst.unl_lof0714)(fst.unl_lof0714)]
binswidefst[, fstunlperc.lof1114 := ecdf(fst.unl_lof1114)(fst.unl_lof1114)]

# pi and theta change
binswidepiD <- dcast(bins[(WinCenter %% 25000) == 0, .(chr, WinCenter, pop, nloci.unl, tPd, tDd, tPd.p, tDd.p)], chr + WinCenter ~ pop, value.var = c('tPd', 'tDd', 'tPd.p', 'tDd.p'))
nrow(binswidepiD)
binswidepiD[, tPdperc.can := ecdf(tPd_can)(tPd_can)]
binswidepiD[, tPdperc.lof0711 := ecdf(tPd_lof0711)(tPd_lof0711)]
binswidepiD[, tPdperc.lof0714 := ecdf(tPd_lof0714)(tPd_lof0714)]
binswidepiD[, tPdperc.lof1114 := ecdf(tPd_lof1114)(tPd_lof1114)]

binswidepiD[, tDdperc.can := ecdf(tDd_can)(tDd_can)]
binswidepiD[, tDdperc.lof0711 := ecdf(tDd_lof0711)(tDd_lof0711)]
binswidepiD[, tDdperc.lof0714 := ecdf(tDd_lof0714)(tDd_lof0714)]
binswidepiD[, tDdperc.lof1114 := ecdf(tDd_lof1114)(tDd_lof1114)]

# ld change
# pi and theta
binswideld <- dcast(bins[(WinCenter %% 25000) == 0, .(chr, WinCenter, pop, nloci.unl, ld_diff)], chr + WinCenter ~ pop, value.var = c('ld_diff'))
setnames(binswideld, c('can', 'lof0711', 'lof0714', 'lof1114'), c('ld_diff_can', 'ld_diff_lof0711', 'ld_diff_lof0714', 'ld_diff_lof1114'))
nrow(binswideld)
binswideld[, ld_diffperc.can := ecdf(ld_diff_can)(ld_diff_can)]
binswideld[, ld_diffperc.lof0711 := ecdf(ld_diff_lof0711)(ld_diff_lof0711)]
binswideld[, ld_diffperc.lof0714 := ecdf(ld_diff_lof0714)(ld_diff_lof0714)]
binswideld[, ld_diffperc.lof1114 := ecdf(ld_diff_lof1114)(ld_diff_lof1114)]

# combine all together
binswide <- merge(binswidefst, binswidepiD, by = c('chr', 'WinCenter'), all = TRUE)
binswide <- merge(binswide, binswideld, by = c('chr', 'WinCenter'), all = TRUE)

nrow(binswide)
nrow(binswidefst)
nrow(binswidepiD)
nrow(binswideld)

######################
## Biplots of change
## Within same kind of metric, across populations
######################

# make plots of shared outliers
cols <- c('grey', 'purple', 'red', 'black')

#quartz(height=6, width=6)
outfile=paste('figures/region_change_across_pops_5e4.png', sep='')
outfile
png(height=6, width=6, units='in', res=300, file=outfile)
par(mfrow=c(4, 3), las=1, mai=c(0.5, 0.6, 0.1, 0.1), tcl=-0.2, mgp=c(1.5,0.4,0), cex.axis=0.7)

# fst
binswide[, plot(fst.unl_can, fst.unl_lof0711, cex=0.2, col=cols[1], pch=16)]
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99, points(fst.unl_can, fst.unl_lof0711, pch=16, cex=0.4, col = cols[2])]
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99 & fstunlperc.lof0714 > 0.99, points(fst.unl_can, fst.unl_lof0711, pch=16, cex=0.4, col = cols[3])]
binswide[fstunlperc.lof1114 > 0.99, points(fst.unl_can, fst.unl_lof0711, pch=16, cex=0.4, col = cols[4])]

legend('bottomright', col=cols, pch=16, legend=c('Locus', '1% outlier in both comparisons', '1% outlier in 3 comparisons', '1% outlier in lof1114'), cex = 0.3)

binswide[, plot(fst.unl_can, fst.unl_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0714 > 0.99, points(fst.unl_can, fst.unl_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99 & fstunlperc.lof0714 > 0.99, points(fst.unl_can, fst.unl_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[fstunlperc.lof1114 > 0.99, points(fst.unl_can, fst.unl_lof0711, pch=16, cex=0.4, col = cols[4])]

binswide[, plot(fst.unl_lof0711, fst.unl_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[fstunlperc.lof0711 > 0.99 & fstunlperc.lof0714 > 0.99, points(fst.unl_lof0711, fst.unl_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99 & fstunlperc.lof0714 > 0.99, points(fst.unl_lof0711, fst.unl_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[fstunlperc.lof1114 > 0.99, points(fst.unl_lof0711, fst.unl_lof0714, pch=16, cex=0.4, col = cols[4])]

# pi change
binswide[, plot(tPd_can, tPd_lof0711, cex=0.2, col=cols[1], pch=16)]
binswide[(tPdperc.can > 0.995 & tPdperc.lof0711 > 0.995) | (tPdperc.can < 0.005 & tPdperc.lof0711 < 0.005), 
            points(tPd_can, tPd_lof0711, pch=16, cex=0.4, col = cols[2])]
binswide[(tPdperc.can > 0.995 & tPdperc.lof0711 > 0.995 & tPdperc.lof0714 > 0.995) | (tPdperc.can < 0.005 & tPdperc.lof0711 < 0.005 & tPdperc.lof0714 < 0.005), 
            points(tPd_can, tPd_lof0711, pch=16, cex=0.4, col = cols[3])]
binswide[tPdperc.lof1114 > 0.995 | tPdperc.lof1114 < 0.005, points(tPd_can, tPd_lof0711, pch=16, cex=0.4, col = cols[4])]

binswide[, plot(tPd_can, tPd_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[(tPdperc.can > 0.995 & tPdperc.lof0714 > 0.995) | (tPdperc.can < 0.005 & tPdperc.lof0714 < 0.005), 
            points(tPd_can, tPd_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[(tPdperc.can > 0.995 & tPdperc.lof0711 > 0.995 & tPdperc.lof0714 > 0.995) | (tPdperc.can < 0.005 & tPdperc.lof0711 < 0.005 & tPdperc.lof0714 < 0.005), 
            points(tPd_can, tPd_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[tPdperc.lof1114 > 0.995 | tPdperc.lof1114 < 0.005, points(tPd_can, tPd_lof0711, pch=16, cex=0.4, col = cols[4])]

binswide[, plot(tPd_lof0711, tPd_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[(tPdperc.lof0711 > 0.995 & tPdperc.lof0714 > 0.995) | (tPdperc.lof0711 < 0.005 & tPdperc.lof0714 < 0.005), points(tPd_lof0711, tPd_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[(tPdperc.can > 0.995 & tPdperc.lof0711 > 0.995 & tPdperc.lof0714 > 0.995) | (tPdperc.can < 0.005 & tPdperc.lof0711 < 0.005 & tPdperc.lof0714 < 0.005), 
            points(tPd_lof0711, tPd_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[tPdperc.lof1114 > 0.995 | tPdperc.lof1114 < 0.005, points(tPd_lof0711, tPd_lof0714, pch=16, cex=0.4, col = cols[4])]

# tajd change
binswide[, plot(tDd_can, tDd_lof0711, cex=0.2, col=cols[1], pch=16)]
binswide[(tDdperc.can > 0.995 & tDdperc.lof0711 > 0.995) | (tDdperc.can < 0.005 & tDdperc.lof0711 < 0.005), 
            points(tDd_can, tDd_lof0711, pch=16, cex=0.4, col = cols[2])]
binswide[(tDdperc.can > 0.995 & tDdperc.lof0711 > 0.995 & tDdperc.lof0714 > 0.995) | (tDdperc.can < 0.005 & tDdperc.lof0711 < 0.005 & tDdperc.lof0714 < 0.005), 
            points(tDd_can, tDd_lof0711, pch=16, cex=0.4, col = cols[3])]
binswide[tDdperc.lof1114 > 0.995 | tDdperc.lof1114 < 0.005, points(tDd_can, tDd_lof0711, pch=16, cex=0.4, col = cols[4])]

binswide[, plot(tDd_can, tDd_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[(tDdperc.can > 0.995 & tDdperc.lof0714 > 0.995) | (tDdperc.can < 0.005 & tDdperc.lof0714 < 0.005), 
            points(tDd_can, tDd_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[(tDdperc.can > 0.995 & tDdperc.lof0711 > 0.995 & tDdperc.lof0714 > 0.995) | (tDdperc.can < 0.005 & tDdperc.lof0711 < 0.005 & tDdperc.lof0714 < 0.005), 
            points(tDd_can, tDd_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[tDdperc.lof1114 > 0.995 | tDdperc.lof1114 < 0.005, points(tDd_can, tDd_lof0711, pch=16, cex=0.4, col = cols[4])]

binswide[, plot(tDd_lof0711, tDd_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[(tDdperc.lof0711 > 0.995 & tDdperc.lof0714 > 0.995) | (tDdperc.lof0711 < 0.005 & tDdperc.lof0714 < 0.005), 
            points(tDd_lof0711, tDd_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[(tDdperc.can > 0.995 & tDdperc.lof0711 > 0.995 & tDdperc.lof0714 > 0.995) | (tDdperc.can < 0.005 & tDdperc.lof0711 < 0.005 & tDdperc.lof0714 < 0.005), 
            points(tDd_lof0711, tDd_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[tDdperc.lof1114 > 0.995 | tDdperc.lof1114 < 0.005, points(tDd_lof0711, tDd_lof0714, pch=16, cex=0.4, col = cols[4])]

# LD change
binswide[, plot(ld_diff_can, ld_diff_lof0711, cex=0.2, col=cols[1], pch=16)]
binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0711 > 0.99, points(ld_diff_can, ld_diff_lof0711, pch=16, cex=0.4, col = cols[2])]
binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0711 > 0.99 & ld_diffperc.lof0714 > 0.99, points(ld_diff_can, ld_diff_lof0711, pch=16, cex=0.4, col = cols[3])]
binswide[ld_diffperc.lof1114 > 0.99, points(ld_diff_can, ld_diff_lof0711, pch=16, cex=0.4, col = cols[4])]

binswide[, plot(ld_diff_can, ld_diff_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0714 > 0.99, points(ld_diff_can, ld_diff_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0711 > 0.99 & ld_diffperc.lof0714 > 0.99, points(ld_diff_can, ld_diff_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[ld_diffperc.lof1114 > 0.99, points(ld_diff_can, ld_diff_lof0711, pch=16, cex=0.4, col = cols[4])]

binswide[, plot(ld_diff_lof0711, ld_diff_lof0714, cex=0.2, col=cols[1], pch=16)]
binswide[ld_diffperc.lof0711 > 0.99 & ld_diffperc.lof0714 > 0.99, points(ld_diff_lof0711, ld_diff_lof0714, pch=16, cex=0.4, col = cols[2])]
binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0711 > 0.99 & ld_diffperc.lof0714 > 0.99, points(ld_diff_lof0711, ld_diff_lof0714, pch=16, cex=0.4, col = cols[3])]
binswide[ld_diffperc.lof1114 > 0.99, points(ld_diff_lof0711, ld_diff_lof0714, pch=16, cex=0.4, col = cols[4])]

dev.off()

##########################################
## Examine and output genomic regions
## with shared changes across populations
##########################################

## fst
# check outlier regions
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99 & fstunlperc.lof0714 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)),
            .(chr, WinCenter, nloci.unl_can, nloci.unl_lof0711, fst.unl_can, fst.unl_lof0711, fst.unl_lof0714, fst.unl_lof1114, fst.p_can, fst.p_lof0711, fst.p_lof0714, fst.p_lof1114)]

# binomial test for shared outlier regions
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)), .(chr, WinCenter, nloci_can, nloci_lof0711)] # shared outliers can 0711
binswide[fstunlperc.can > 0.99 & fstunlperc.lof0714 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)), .(chr, WinCenter, nloci_can, nloci_lof0711)] # shared outliers can 0714

nrow(binswide)
n12 <- binswide[!is.na(fstunlperc.can) & !is.na(fstunlperc.lof0711), .N] # number of evaluated regions (can to lof0711)
n13 <- binswide[!is.na(fstunlperc.can) & !is.na(fstunlperc.lof0714), .N] # number of evaluated regions (can to lof0714)
x1 <- binswide[fstunlperc.can > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)), .N] # outliers in can
x2 <- binswide[fstunlperc.lof0711 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)), .N] # outliers in lof0711
x3 <- binswide[fstunlperc.lof0714 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)), .N] # outliers in lof0714
x12 <- binswide[fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)), .N] # shared outliers can 0711
x13 <- binswide[fstunlperc.can > 0.99 & fstunlperc.lof0714 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114)), .N] # shared outliers can 0714
x1/n12 * x2/n12 * n12 # expected number shared can 0711
x12 # number shared
binom.test(x12, n12, p = x1/n12 * x2/n12, alternative = 'two.sided')

x1/n13 * x3/n13 * n13 # expected number shared can 0714
x13 # number shared
binom.test(x13, n13, p = x1/n13 * x3/n13, alternative = 'two.sided')

## pi change
# check outlier regions
binswide[(tPdperc.can > 0.995 & tPdperc.lof0711 > 0.995 & tPdperc.lof0714 > 0.995 & (tPdperc.lof1114 < 0.995 | is.na(tPdperc.lof1114))) |
              (tPdperc.can < 0.005 & tPdperc.lof0711 < 0.005 & tPdperc.lof0714 < 0.005 & (tPdperc.lof1114 > 0.005 | is.na(tPdperc.lof1114))),
            .(chr, WinCenter, nloci.unl_can, nloci.unl_lof0711, tPd_can, tPd_lof0711, tPd_lof0714, tPd_lof1114, tPd.p_can, tPd.p_lof0711, tPd.p_lof0714, tPd.p_lof1114)]

# binomial test for shared outlier regions
nrow(binswide)
n12 <- binswide[!is.na(tPdperc.can) & !is.na(tPdperc.lof0711), .N] # number of evaluated regions (can to lof0711)
n13 <- binswide[!is.na(tPdperc.can) & !is.na(tPdperc.lof0714), .N] # number of evaluated regions (can to lof0714)
x1 <- binswide[(tPdperc.can > 0.995 & (tPdperc.lof1114 < 0.995 | is.na(tPdperc.lof1114))) |
                    tPdperc.can < 0.005 & (tPdperc.lof1114 > 0.005 | is.na(tPdperc.lof1114)), .N] # outliers in can
x2 <- binswide[(tPdperc.lof0711 > 0.995 & (tPdperc.lof1114 < 0.995 | is.na(tPdperc.lof1114))) |
                    (tPdperc.lof0711 < 0.005 & (tPdperc.lof1114 > 0.005 | is.na(tPdperc.lof1114))), .N] # outliers in lof0711
x3 <- binswide[(tPdperc.lof0714 > 0.995 & (tPdperc.lof1114 < 0.995 | is.na(tPdperc.lof1114))) |
                    (tPdperc.lof0714 < 0.005 & (tPdperc.lof1114 > 0.005 | is.na(tPdperc.lof1114))), .N] # outliers in lof0714
x12 <- binswide[(tPdperc.can > 0.995 & tPdperc.lof0711 > 0.995 & (tPdperc.lof1114 < 0.995 | is.na(tPdperc.lof1114))) |
                     (tPdperc.can < 0.005 & tPdperc.lof0711 < 0.005 & (tPdperc.lof1114 > 0.005 | is.na(tPdperc.lof1114))), .N] # shared outliers can 0711
x13 <- binswide[(tPdperc.can > 0.995 & tPdperc.lof0714 > 0.995 & (tPdperc.lof1114 < 0.995 | is.na(tPdperc.lof1114))) |
                     (tPdperc.can < 0.005 & tPdperc.lof0714 < 0.005 & (tPdperc.lof1114 > 0.005 | is.na(tPdperc.lof1114))), .N] # shared outliers can 0714
x1/n12 * x2/n12 * n12 # expected number shared can 0711
x12 # number shared
binom.test(x12, n12, p = x1/n12 * x2/n12, alternative = 'two.sided')

x1/n13 * x3/n13 * n13 # expected number shared can 0714
x13 # number shared
binom.test(x13, n13, p = x1/n13 * x3/n13, alternative = 'two.sided')


## D change
# check outlier regions
binswide[(tDdperc.can > 0.995 & tDdperc.lof0711 > 0.995 & tDdperc.lof0714 > 0.995 & (tDdperc.lof1114 < 0.995 | is.na(tDdperc.lof1114))) |
              (tDdperc.can < 0.005 & tDdperc.lof0711 < 0.005 & tDdperc.lof0714 < 0.005 & (tDdperc.lof1114 > 0.005 | is.na(tDdperc.lof1114))),
            .(chr, WinCenter, nloci.unl_can, nloci.unl_lof0711, tDd_can, tDd_lof0711, tDd_lof0714, tDd_lof1114, tDd.p_can, tDd.p_lof0711, tDd.p_lof0714, tDd.p_lof1114)]

# binomial test for shared outlier regions
nrow(binswide)
n12 <- binswide[!is.na(tDdperc.can) & !is.na(tDdperc.lof0711), .N] # number of evaluated regions (can to lof0711)
n13 <- binswide[!is.na(tDdperc.can) & !is.na(tDdperc.lof0714), .N] # number of evaluated regions (can to lof0714)
x1 <- binswide[(tDdperc.can > 0.995 & (tDdperc.lof1114 < 0.995 | is.na(tDdperc.lof1114))) |
                    tDdperc.can < 0.005 & (tDdperc.lof1114 > 0.005 | is.na(tDdperc.lof1114)), .N] # outliers in can
x2 <- binswide[(tDdperc.lof0711 > 0.995 & (tDdperc.lof1114 < 0.995 | is.na(tDdperc.lof1114))) |
                    (tDdperc.lof0711 < 0.005 & (tDdperc.lof1114 > 0.005 | is.na(tDdperc.lof1114))), .N] # outliers in lof0711
x3 <- binswide[(tDdperc.lof0714 > 0.995 & (tDdperc.lof1114 < 0.995 | is.na(tDdperc.lof1114))) |
                    (tDdperc.lof0714 < 0.005 & (tDdperc.lof1114 > 0.005 | is.na(tDdperc.lof1114))), .N] # outliers in lof0714
x12 <- binswide[(tDdperc.can > 0.995 & tDdperc.lof0711 > 0.995 & (tDdperc.lof1114 < 0.995 | is.na(tDdperc.lof1114))) |
                     (tDdperc.can < 0.005 & tDdperc.lof0711 < 0.005 & (tDdperc.lof1114 > 0.005 | is.na(tDdperc.lof1114))), .N] # shared outliers can 0711
x13 <- binswide[(tDdperc.can > 0.995 & tDdperc.lof0714 > 0.995 & (tDdperc.lof1114 < 0.995 | is.na(tDdperc.lof1114))) |
                     (tDdperc.can < 0.005 & tDdperc.lof0714 < 0.005 & (tDdperc.lof1114 > 0.005 | is.na(tDdperc.lof1114))), .N] # shared outliers can 0714
x1/n12 * x2/n12 * n12 # expected number shared can 0711
x12 # number shared
binom.test(x12, n12, p = x1/n12 * x2/n12, alternative = 'two.sided')

x1/n13 * x3/n13 * n13 # expected number shared can 0714
x13 # number shared
binom.test(x13, n13, p = x1/n13 * x3/n13, alternative = 'two.sided')


## LD Change
# check outlier regions
binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0711 > 0.99 & ld_diffperc.lof0714 > 0.99 & (ld_diffperc.lof1114 < 0.99 | is.na(ld_diffperc.lof1114)),
            .(chr, WinCenter, nloci.unl_can, nloci.unl_lof0711, ld_diff_can, ld_diff_lof0711, ld_diff_lof0714, ld_diff_lof1114)]

# binomial test for shared outlier regions
nrow(binswide)
n12 <- binswide[!is.na(ld_diffperc.can) & !is.na(ld_diffperc.lof0711), .N] # number of evaluated regions (can to lof0711)
n13 <- binswide[!is.na(ld_diffperc.can) & !is.na(ld_diffperc.lof0714), .N] # number of evaluated regions (can to lof0714)
x1 <- binswide[ld_diffperc.can > 0.99 & (ld_diffperc.lof1114 < 0.99 | is.na(ld_diffperc.lof1114)), .N] # outliers in can
x2 <- binswide[ld_diffperc.lof0711 > 0.99 & (ld_diffperc.lof1114 < 0.99 | is.na(ld_diffperc.lof1114)), .N] # outliers in lof0711
x3 <- binswide[ld_diffperc.lof0714 > 0.99 & (ld_diffperc.lof1114 < 0.99 | is.na(ld_diffperc.lof1114)), .N] # outliers in lof0714
x12 <- binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0711 > 0.99 & (ld_diffperc.lof1114 < 0.99 | is.na(ld_diffperc.lof1114)), .N] # shared outliers can 0711
x13 <- binswide[ld_diffperc.can > 0.99 & ld_diffperc.lof0714 > 0.99 & (ld_diffperc.lof1114 < 0.99 | is.na(ld_diffperc.lof1114)), .N] # shared outliers can 0714
x1/n12 * x2/n12 * n12 # expected number shared can 0711
x12 # number shared
binom.test(x12, n12, p = x1/n12 * x2/n12, alternative = 'two.sided')

x1/n13 * x3/n13 * n13 # expected number shared can 0714
x13 # number shared
binom.test(x13, n13, p = x1/n13 * x3/n13, alternative = 'two.sided')

## output FST outlier regions
out <- binswide[(fstunlperc.can > 0.99 & fstunlperc.lof0711 > 0.99 & fstunlperc.lof0714 > 0.99 & (fstunlperc.lof1114 < 0.99 | is.na(fstunlperc.lof1114))), ]

write.csv(out, file = gzfile('analysis/outlier_50kregions_shared_07-11-14_Can.csv.gz'), row.names = FALSE)

######################
## Biplots of change
## Within same population, across metrics
######################

cols <- c('grey', 'purple', 'red')

quartz(height=6, width=6)
outfile=paste('figures/region_change_across_metrics_', wdnm, '.png', sep='')
outfile
# png(height=6, width=6, units='in', res=300, file=outfile)
par(mfrow=c(4, 3), las=1, mai=c(0.5, 0.6, 0.1, 0.1), tcl=-0.2, mgp=c(2.8,0.4,0), cex.axis=0.7)
# frequency change vs. pi change
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[freq_region_perc0711<=0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, plot(freq_diff_0711, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region_perc0711>0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, points(freq_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(freq_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(freq_diff_0711, pi_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region_perc0714<=0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, plot(freq_diff_0714, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region_perc0714>0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, points(freq_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714<=0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(freq_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714>0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(freq_diff_0714, pi_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region_percCAN<=0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, plot(freq_diff_CAN, pi_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_percCAN>0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, points(freq_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN<=0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(freq_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN>0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(freq_diff_CAN, pi_diff_CAN, cex=0.4, col=cols[3], pch=16)]

# freq vs. D
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(D_diff_CAN, D_diff_0711, D_diff_0714, na.rm=TRUE)]
bins[freq_region_perc0711<=0.99 & D_region_perc0711<=0.995 & D_region_perc0711>=0.005, plot(freq_diff_0711, D_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region_perc0711>0.99 & D_region_perc0711<=0.995 & D_region_perc0711>=0.005, points(freq_diff_0711, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & (D_region_perc0711>0.995 | D_region_perc0711<0.005), points(freq_diff_0711, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & (D_region_perc0711>0.995 | D_region_perc0711<0.005), points(freq_diff_0711, D_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region_perc0714<=0.99 & D_region_perc0714<=0.995 & D_region_perc0714>=0.005, plot(freq_diff_0714, D_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region_perc0714>0.99 & D_region_perc0714<=0.995 & D_region_perc0714>=0.005, points(freq_diff_0714, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714<=0.99 & (D_region_perc0714>0.995 | D_region_perc0714<0.005), points(freq_diff_0714, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714>0.99 & (D_region_perc0714>0.995 | D_region_perc0714<0.005), points(freq_diff_0714, D_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region_percCAN<=0.99 & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, plot(freq_diff_CAN, D_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_percCAN>0.99 & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, points(freq_diff_CAN, D_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN<=0.99 & (D_region_percCAN>0.995 | D_region_percCAN<0.005), points(freq_diff_CAN, D_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN>0.99 & (D_region_percCAN>0.995 | D_region_percCAN<0.005), points(freq_diff_CAN, D_diff_CAN, cex=0.4, col=cols[3], pch=16)]

# freq vs. ld
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
bins[freq_region_perc0711<=0.99 & ld_region_perc0711<=0.99, plot(freq_diff_0711, ld_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region_perc0711>0.99 & ld_region_perc0711<=0.99, points(freq_diff_0711, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & ld_region_perc0711>0.99, points(freq_diff_0711, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & ld_region_perc0711>0.99, points(freq_diff_0711, ld_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region_perc0714<=0.99 & ld_region_perc0714<=0.99, plot(freq_diff_0714, ld_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region_perc0714>0.99 & ld_region_perc0714<=0.99, points(freq_diff_0714, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714<=0.99 & ld_region_perc0714>0.99, points(freq_diff_0714, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714>0.99 & ld_region_perc0714>0.99, points(freq_diff_0714, ld_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region_percCAN<=0.99 & ld_region_percCAN<=0.99, plot(freq_diff_CAN, ld_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_percCAN>0.99 & ld_region_percCAN<=0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN<=0.99 & ld_region_percCAN>0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN>0.99 & ld_region_percCAN>0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.4, col=cols[3], pch=16)]


# ld vs. pi
xlims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[ld_region_perc0711<=0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, plot(ld_diff_0711, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0711', line=1.5)
bins[ld_region_perc0711>0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, points(ld_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711<=0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(ld_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711>0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(ld_diff_0711, pi_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[ld_region_perc0714<=0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, plot(ld_diff_0714, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0714', line=1.5)
bins[ld_region_perc0714>0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, points(ld_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0714<=0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(ld_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0714>0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(ld_diff_0714, pi_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[ld_region_percCAN<=0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, plot(ld_diff_CAN, pi_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_CAN', line=1.5)
bins[ld_region_percCAN>0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, points(ld_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_percCAN<=0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(ld_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_percCAN>0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(ld_diff_CAN, pi_diff_CAN, cex=0.4, col=cols[3], pch=16)]

dev.off()



	
	
###############################
## Plot FST in outlier regions
###############################
# read in Fsts
fstLof1 <- fread('analysis/Lof_07.Lof_11.fst.AB.gz', col.names = c('CHROM', 'POS', 'A', 'B')) # from angsd_fst.sh
fstLof2 <- fread('analysis/Lof_07.Lof_14.fst.AB.gz', col.names = c('CHROM', 'POS', 'A', 'B'))
fstCan <- fread('analysis/Can_40.Can_14.fst.AB.gz', col.names = c('CHROM', 'POS', 'A', 'B'))

# trim fsts to nodam2 loci
gatk <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab', col.names = c('CHROM', 'POS', 'REF', 'ALT'))
fstLof1 <- merge(fstLof1, gatk[, .(CHROM, POS)])
fstLof2 <- merge(fstLof2, gatk[, .(CHROM, POS)])
fstCan <- merge(fstCan, gatk[, .(CHROM, POS)])

# calc fst
fstCan[, fst := A/B]
fstLof1[, fst := A/B]
fstLof2[, fst := A/B]


# list outlier 50kb regions to plot
dat <- data.table(CHROM = c('LG10', 'LG11', 'LG14', 'LG20'), POS = c( 13375000, 175000, 8425000, 16725000))

# make plot
ncol = 3
nrow = ceiling(nrow(dat)/ncol)
cols = c('black', 'blue', 'green') # can, 0711, 0714
rng <- 25000

par(mfrow = c(nrow, ncol), mai = c(0.3, 0.3, 0.4, 0.1), omi = c(0.3, 0.3, 0, 0))
for(i in 1:nrow(dat)){
  fstCan[CHROM == dat$CHROM[i] & abs(POS - dat$POS[i]) < rng, plot(POS/1e6, fst, xlab = 'Mb', ylab = 'Fst', main = dat$CHROM[i], col = cols[1])]
  fstLof1[CHROM == dat$CHROM[i] & abs(POS - dat$POS[i]) < rng, points(POS/1e6, fst, xlab = 'Mb', ylab = 'Fst', main = dat$CHROM[i], col = cols[2])]
  fstLof2[CHROM == dat$CHROM[i] & abs(POS - dat$POS[i]) < rng, points(POS/1e6, fst, xlab = 'Mb', ylab = 'Fst', main = dat$CHROM[i], col = cols[3])]
}
