# run OutFLANK

library(data.table)
library(OutFLANK)

# function to convert 3 beagle columns to 012
beagleto012 <- function(dat){
  setnames(dat, c('a', 'b', 'c'))
  out <- dat[, which.max(c(a, b, c)) - 1, by = seq_len(nrow(dat))]
  return(out[,V1])
}


# beagle genotypes
genos <- fread('data_2020.05.07/GATK_no_dam2.beagle.gz')
inds <- fread('data_31_01_20/List_to_malin.tab')
unlinkcan <- fread('analysis/ld.unlinked.Can.gatk.nodam.csv.gz')
unlinklof <- fread('analysis/ld.unlinked.Lof.gatk.nodam.csv.gz')

# find indices of unlinked
unlinkcan[, beaglenames := paste0(CHROM, '_', POS), by = seq_len(nrow(unlinkcan))]
unlinklof[, beaglenames := paste0(CHROM, '_', POS), by = seq_len(nrow(unlinklof))]
trim_can <- which(genos$marker %in% unlinkcan$beaglenames)
trim_lof <- which(genos$marker %in% unlinklof$beaglenames)


# convert beagle to 012 format by using highest prob genotype
ncols <- ncol(genos)
indnms <- colnames(genos)[seq(4, ncols, by = 3)]
newnms <- c('marker', 'allele1', 'allele2', paste0(colnames(genos)[4:ncols], c('.1', '.2', '.3')))
setnames(genos, newnms)
genos012 <- genos[, .(marker)]
for(i in 1:length(indnms)){
  cat(i)
  thesecols <- paste0(indnms[i], c('.1', '.2', '.3'))
  newnm <- indnms[i]
  genos012[[newnm]] <- beagleto012(genos[, ..thesecols])
}

# separate into pops
inds[, ind := paste0('Ind', seq(0, nrow(inds)-1))]
caninds <- inds[Pop %in% c('Historic_Canada', 'Canada_2013'), ind]
lof0711inds <- inds[Pop %in% c('Historic_Lofoten', 'Lofoten_2011'), ind]
lof0714inds <- inds[Pop %in% c('Historic_Lofoten', 'Lofoten_2014'), ind]
lof1114inds <- inds[Pop %in% c('Lofoten_2011', 'Lofoten_2014'), ind]

canpops <- inds[Pop %in% c('Historic_Canada', 'Canada_2013'), Pop]
lof0711pops <- inds[Pop %in% c('Historic_Lofoten', 'Lofoten_2011'), Pop]
lof0714pops <- inds[Pop %in% c('Historic_Lofoten', 'Lofoten_2014'), Pop]
lof1114pops <- inds[Pop %in% c('Lofoten_2011', 'Lofoten_2014'), Pop]

can <- genos012[, ..caninds]
lof0711 <- genos012[, ..lof0711inds]
lof0714 <- genos012[, ..lof0714inds]
lof1114 <- genos012[, ..lof1114inds]

# prep data for outflank
matcan <- t(as.matrix(can))
matlof0711 <- t(as.matrix(lof0711))
matlof0714 <- t(as.matrix(lof0714))
matlof1114 <- t(as.matrix(lof1114))

locusNames <- genos[, marker]
colnames(matcan) <- locusNames
colnames(matlof0711) <- locusNames
colnames(matlof0714) <- locusNames
colnames(matlof1114) <- locusNames

# calc fst
fstcan <- MakeDiploidFSTMat(matcan, locusNames, canpops)
fstlof0711 <- MakeDiploidFSTMat(matlof0711, locusNames, lof0711pops)
fstlof0714 <- MakeDiploidFSTMat(matlof0714, locusNames, lof0714pops)
fstlof1114 <- MakeDiploidFSTMat(matlof1114, locusNames, lof1114pops)

# data check: He vs. Fst
# want to remove loci with low H loci but high FST. 
plot(fstcan$He, fstcan$FST)
plot(fstlof0711$He, fstlof0711$FST)
plot(fstlof0714$He, fstlof0714$FST)
plot(fstlof1114$He, fstlof1114$FST)

# Data checks: FST vs. FSTNoCorr
# "Look for loci that deviate from the linear relationship in this plot, and remove those loci."
plot(fstcan$FST, fstcan$FSTNoCorr); abline(0,1)
plot(fstlof0711$FST, fstlof0711$FSTNoCorr); abline(0,1)
plot(fstlof0714$FST, fstlof0714$FSTNoCorr); abline(0,1)
plot(fstlof1114$FST, fstlof1114$FSTNoCorr); abline(0,1)


# run outflank on trimmed loci
outcan <- OutFLANK(FstDataFrame=fstcan[trim_can,], LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05) 
outlof0711 <- OutFLANK(FstDataFrame=fstlof0711[trim_lof, ], LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
outlof0714 <- OutFLANK(FstDataFrame=fstlof0714[trim_lof, ], LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
outlof1114 <- OutFLANK(FstDataFrame=fstlof1114[trim_lof, ], LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)


#Check the fit and make sure it looks good, especially in the right tail:
OutFLANKResultsPlotter(outcan, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)
OutFLANKResultsPlotter(outlof0711, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)
OutFLANKResultsPlotter(outlof0714, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)
OutFLANKResultsPlotter(outlof1114, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)

## Zoom in on right tail
OutFLANKResultsPlotter(outcan, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = TRUE, RightZoomFraction = 0.15, titletext = NULL)
OutFLANKResultsPlotter(outlof0711, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = TRUE, RightZoomFraction = 0.15, titletext = NULL)
OutFLANKResultsPlotter(outlof0714, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = TRUE, RightZoomFraction = 0.15, titletext = NULL)
OutFLANKResultsPlotter(outlof1114, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, 
                       Zoom = TRUE, RightZoomFraction = 0.15, titletext = NULL)

# check the P-value histogram:
# we expect this histogram to be flat and maybe have a bump near 0 for selected sites
hist(outcan$results$pvaluesRightTail)
hist(outlof0711$results$pvaluesRightTail)
hist(outlof0714$results$pvaluesRightTail)
hist(outlof1114$results$pvaluesRightTail)

# Using estimated neutral mean FST and df to calculate P-values for all loci
Pcan <- pOutlierFinderChiSqNoCorr(fstcan, Fstbar = outcan$FSTNoCorrbar, 
                                  dfInferred = outcan$dfInferred, qthreshold = 0.05, Hmin = 0.1)
Plof0711 <- pOutlierFinderChiSqNoCorr(fstlof0711, Fstbar = outlof0711$FSTNoCorrbar, 
                                  dfInferred = outlof0711$dfInferred, qthreshold = 0.05, Hmin = 0.1)
Plof0714 <- pOutlierFinderChiSqNoCorr(fstlof0714, Fstbar = outlof0714$FSTNoCorrbar, 
                                      dfInferred = outlof0714$dfInferred, qthreshold = 0.05, Hmin = 0.1)
Plof1114 <- pOutlierFinderChiSqNoCorr(fstlof1114, Fstbar = outlof1114$FSTNoCorrbar, 
                                      dfInferred = outlof1114$dfInferred, qthreshold = 0.05, Hmin = 0.1)
can_out <- Pcan$OutlierFlag==TRUE
lof0711_out <- Plof0711$OutlierFlag==TRUE
lof0714_out <- Plof0714$OutlierFlag==TRUE
lof1114_out <- Plof1114$OutlierFlag==TRUE

sum(can_out, na.rm = TRUE) # number of outlier loci (all loci)
outcan$numberHighFstOutliers # number of outlier loci (LD-trimmed loci)

sum(lof0711_out, na.rm = TRUE)
outlof0711$numberHighFstOutliers

sum(lof0714_out, na.rm = TRUE)
outlof0714$numberHighFstOutliers

sum(lof1114_out, na.rm = TRUE)
outlof1114$numberHighFstOutliers

# plot the outliers
plot(Pcan$He, Pcan$FST, pch=19, col=rgb(0,0,0,0.1))
points(Pcan$He[can_out], Pcan$FST[can_out], col="blue")

hist(Pcan$pvaluesRightTail) # p-values

# write out
write.csv(Pcan, file = gzfile('analysis/angsd_outflank.Can.csv.gz'))
write.csv(Plof0711, file = gzfile('analysis/angsd_outflank.Lof0711.csv.gz'))
write.csv(Plof0714, file = gzfile('analysis/angsd_outflank.Lof0714.csv.gz'))
write.csv(Plof1114, file = gzfile('analysis/angsd_outflank.Lof1114.csv.gz'))


#####################
# plot
#####################
library(data.table)

# read in results
Pcan <- fread('analysis/angsd_outflank.Can.csv.gz')
Plof0711 <- fread('analysis/angsd_outflank.Lof0711.csv.gz')
Plof0714 <- fread('analysis/angsd_outflank.Lof0714.csv.gz')
Plof1114 <- fread('analysis/angsd_outflank.Lof1114.csv.gz')

# add POS and genome pos
Pcan[, c('chr', 'pos') := tstrsplit(LocusName, "_")]
Plof0711[, c('chr', 'pos') := tstrsplit(LocusName, "_")]
Plof0714[, c('chr', 'pos') := tstrsplit(LocusName, "_")]
Plof1114[, c('chr', 'pos') := tstrsplit(LocusName, "_")]

chrmax <- fread('data/lg_length.csv')
chrmax$start=c(0, cumsum(chrmax$len)[1:(nrow(chrmax) - 1)])

Pcan <- merge(Pcan, chrmax[,.(chr, start)], by='chr')
Pcan[,POSgen := as.numeric(pos) + start]
Pcan[,start := NULL]

Plof0711 <- merge(Plof0711, chrmax[,.(chr, start)], by='chr')
Plof0711[,POSgen := as.numeric(pos) + start]
Plof0711[,start := NULL]

Plof0714 <- merge(Plof0714, chrmax[,.(chr, start)], by='chr')
Plof0714[,POSgen := as.numeric(pos) + start]
Plof0714[,start := NULL]

Plof1114 <- merge(Plof1114, chrmax[,.(chr, start)], by='chr')
Plof1114[,POSgen := as.numeric(pos) + start]
Plof1114[,start := NULL]

# add a vector for color by LG
cols <- c('#a6cee3aa', '#1f78b4aa') # light blue, blue, partially transparent: for alternating LGs
Pcan[, lgcol := cols[1]]
Pcan[chr %in% chrmax$chr[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]

Plof0711[, lgcol := cols[1]]
Plof0711[chr %in% chrmax$chr[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]

Plof0714[, lgcol := cols[1]]
Plof0714[chr %in% chrmax$chr[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]

Plof1114[, lgcol := cols[1]]
Plof1114[chr %in% chrmax$chr[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]

# plot
png(filename = paste0('figures/angsd_OutFLANK_gatkNoDam.png'), width = 20, height = 12, units = "in", res = 300)
par(mfrow = c(4,1), mai = c(0.7, 1, 0.2, 0.1))
Pcan[, plot(POSgen/1e6, -log10(qvalues), xlab = 'Position (Mb)', ylab = "-log10(FDR-adjusted p)", 
            main = 'Canada', col = lgcol, cex = 0.5)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

Plof0711[, plot(POSgen/1e6, -log10(qvalues), xlab = 'Position (Mb)', ylab = "-log10(FDR-adjusted p)",
                main = 'Lof 07-11', col = lgcol, cex = 0.5)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

Plof0714[, plot(POSgen/1e6, -log10(qvalues), xlab = 'Position (Mb)', ylab = "-log10(FDR-adjusted p)",
                main = 'Lof 07-14', col = lgcol, cex = 0.5)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

Plof1114[, plot(POSgen/1e6, -log10(qvalues), xlab = 'Position (Mb)', ylab = "-log10(FDR-adjusted p)",
                main = 'Lof 11-14', col = lgcol, cex = 0.5)]
abline(h = -log10(0.05), lty = 2, col = 'grey')

dev.off()
