# run OutFLANK

library(data.table)
library(OutFLANK)

snpmat <- fread('gzcat data_2018.09.05/All_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_T5PERC_150bpVAR.012.gz')
loc <- fread('gzcat data_2018.09.05/All_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_T5PERC_150bpVAR.012.pos.gz')
indivs <- fread('data_2018.09.05/All_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_T5PERC_150bpVAR.012.indv', header=FALSE)

# drop first column: just a row name
snpmat[,V1:=NULL]

# turn -1 to 9 for missing data
snpmat2 <- as.matrix(snpmat)
snpmat2[snpmat2 == -1] <- 9

sort(unique(as.numeric(snpmat2))) # check that only 0,1,2,9

# make loc names
loc[,locnm := paste(V1, V2, sep='_')]

# make populations
pops <- rep('', nrow(indivs))
pops[grepl('BM_2', indivs$V1)] <- 'LOF_07'
pops[grepl('BM_0|BM_1', indivs$V1)] <- 'CAN40'
pops[grepl('LOF11', indivs$V1)] <- 'LOF_11'
pops[grepl('LOF14', indivs$V1)] <- 'LOF_14'
pops[grepl('TWI', indivs$V1)] <- 'CANMod'
sort(unique(pops))

# find loci that have >=10 individuals genotyped in each population
ss07 <- apply(snpmat2[pops=='LOF_07',], MARGIN=2, FUN=function(x) sum(x!=9)) # lof07 sample size
ss11 <- apply(snpmat2[pops=='LOF_11',], MARGIN=2, FUN=function(x) sum(x!=9))
ss14 <- apply(snpmat2[pops=='LOF_14',], MARGIN=2, FUN=function(x) sum(x!=9))
loclof <- ss07>=10 & ss11>=10 & ss14>=10; sum(loclof) # combined across populations: 514411

ss40 <- apply(snpmat2[pops=='CAN40',], MARGIN=2, FUN=function(x) sum(x!=9)) # Can
ssMod <- apply(snpmat2[pops=='CANMod',], MARGIN=2, FUN=function(x) sum(x!=9))
loccan <- ss40>=10 & ssMod>=10; sum(loccan) # combined: 447146


# make fst matrix
indlof <- pops %in% c('LOF_07', 'LOF_11', 'LOF_14'); sum(indlof)
indcan <- pops %in% c('CAN40', 'CANMod'); sum(indcan)

fstmatLof <- MakeDiploidFSTMat(snpmat2[indlof,loclof], loc$locnm[loclof], pops[indlof]) # takes a few minutes
fstmatCan <- MakeDiploidFSTMat(snpmat2[indcan,loccan], loc$locnm[loccan], pops[indcan])

# run outflank
outlof <- OutFLANK(FstDataFrame=fstmatLof, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=3, qthreshold=0.2)
outcan <- OutFLANK(FstDataFrame=fstmatCan, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.1) 


# plot neutral FST

OutFLANKResultsPlotter(outlof, withOutliers=TRUE, NoCorr=TRUE, Hmin=0.1, binwidth=0.005, Zoom=FALSE, RightZoomFraction=0.05, titletext=NULL)

OutFLANKResultsPlotter(outlof, withOutliers=TRUE, NoCorr=TRUE, Hmin=0.1, binwidth=0.005, Zoom=TRUE, RightZoomFraction=0.005, titletext=NULL) # just right-hand tail


OutFLANKResultsPlotter(outcan, withOutliers=TRUE, NoCorr=TRUE, Hmin=0.1, binwidth=0.005, Zoom=FALSE, RightZoomFraction=0.05, titletext=NULL)

OutFLANKResultsPlotter(outcan, withOutliers=TRUE, NoCorr=TRUE, Hmin=0.1, binwidth=0.005, Zoom=TRUE, RightZoomFraction=0.005, titletext=NULL) # just right-hand tail


# examine
outlof$dfInferred
outlof$numberHighFstOutliers # 3
outlofinds <- outlof$results$OutlierFlag & !is.na(outlof$results$OutlierFlag) # locus indices for outliers
outlof$results[outlofinds,] # the inferred outliers
	ss07[loclof][outlofinds]
	ss11[loclof][outlofinds]
	ss14[loclof][outlofinds]
	
outcan$dfInferred
outcan$numberHighFstOutliers # 1998
outcaninds <- outcan$results$OutlierFlag & !is.na(outcan$results$OutlierFlag) # locus indices for outliers
outlof$results[outcaninds,] # the inferred outliers
	ss40[loccan][outcaninds]
	ssMod[loccan][outcaninds]


# write out
saveRDS(outlof, file='analysis/outflankLOF.rds')
saveRDS(outcan, file='analysis/outflankCAN.rds')