# Compare Canada, 1907-2011, 1907-2014, and 2011-2014 outlier loci from null model analysis (25k)
# to be run on my mac

require(data.table)

##################
## Read in data
##################

# read in Lof
locnms11 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof11.txt', header=TRUE) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	setnames(locnms11, 3:7, c('cnt07', 'Freq_07', 'cnt11', 'Freq_11', 'ABS_DIFF0711'))
dat11 <- readRDS('analysis/wfs_nullmodel_pos&pvals_07-11.rds') # has the p-values from 1907-2011
	setnames(dat11, c('n', 'p'), c('nLof0711', 'pLof0711'))
locnms14 <- fread('data_2019_03_18/Frequency_table_Lof07_Lof14.txt', header=TRUE)
	setnames(locnms14, 3:7, c('cnt07', 'Freq_07', 'cnt14', 'Freq_14', 'ABS_DIFF0714'))
dat14 <- readRDS('analysis/wfs_nullmodel_pos&pvals_07-14.rds')
	setnames(dat14, c('n', 'p'), c('nLof0714', 'pLof0714'))

	# merge
setkey(dat11, CHROM, POS)
setkey(dat14, CHROM, POS)
datLof <- merge(dat11, dat14, all=TRUE)
	nrow(dat11)
	nrow(dat14)
	nrow(datLof)

setkey(datLof, CHROM, POS)
setkey(locnms11, CHROM, POS)
setkey(locnms14, CHROM, POS)
datLof <- merge(datLof, locnms11[, .(CHROM, POS, cnt07, Freq_07, cnt11, Freq_11)])
datLof <- merge(datLof, locnms14[, .(CHROM, POS, cnt14, Freq_14)])
	nrow(locnms11)
	nrow(locnms14)
	nrow(datLof)

# Read in Canada
locnmsCan <- fread('data_2019_03_18/Frequency_table_CAN_40_TGA.txt', header=TRUE)
	setnames(locnmsCan, 3:7, c('cntCan40', 'Freq_Can40', 'cntCanMod', 'Freq_CanMod', 'ABS_DIFFCan'))
	locnmsCan <- locnmsCan[cntCan40>0 & cntCanMod>0,] # trim to genotypes >0%
datCan <- readRDS('analysis/wfs_nullmodel_pos&pvals_Can.rds')
	setnames(datCan, c('n', 'p'), c('nCan', 'pCan'))

	# compare
	nrow(datCan) # now has 2019 dataset
	nrow(locnmsCan) # matches

	# merge
setkey(locnmsCan, CHROM, POS)
setkey(datCan, CHROM, POS)
datCan <- merge(datCan, locnmsCan)
	nrow(datCan)

##################
# merge Lof and Can
##################

setkey(datLof, CHROM, POS)
setkey(datCan, CHROM, POS)
dat <- merge(datCan, datLof, all=TRUE)[, .(CHROM, POS, cnt07, cnt11, cnt14, cntCan40, cntCanMod, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0711, pLof0714, pCan)]

	# make sure the merge worked
	nrow(dat)
	nrow(dat11)
	nrow(dat14)
	dat[,.(sum(!is.na(Freq_07)), sum(!is.na(Freq_11)), sum(!is.na(Freq_14)))] # should match nrow(dat11) and nrow(dat14)
	nrow(datCan)
	dat[,.(sum(!is.na(Freq_Can40)), sum(!is.na(Freq_CanMod)))] # should match nrow(datCan)
	dat[,sum(!is.na(Freq_07) & !is.na(Freq_Can40))] # number of loci genotyped in both
	

###############################################
# combine p-values for individual loci
# https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/
###############################################

dat[,p.comb0711Can := pchisq(-2 * sum(log(c(pLof0711, pCan))),df=2*2,lower=FALSE), by=1:nrow(dat)] # Lof 07-11 and Can
dat[,p.comb0714Can := pchisq(-2 * sum(log(c(pLof0714, pCan))),df=2*2,lower=FALSE), by=1:nrow(dat)] # Lof 07-14 and Can

##################################
# FDR-correct p-values
##################################

# FDR-correct each population separately
# after masking out unplaced
dat[,q2.Lof0711:=as.double(NA)]
dat[,q2.Lof0714:=as.double(NA)]
dat[,q2.Can:=as.double(NA)]
dat[!(CHROM %in% c('Unplaced')), q2.Lof0711 := p.adjust(pLof0711, method='fdr')] 
dat[!(CHROM %in% c('Unplaced')), q2.Lof0714 := p.adjust(pLof0714, method='fdr')] 
dat[!(CHROM %in% c('Unplaced')), q2.Can := p.adjust(pCan, method='fdr')]

# FDR-correct each population separately
# after masking out inversions and unplaced
dat[,q3.Lof0711:=as.double(NA)]
dat[,q3.Lof0714:=as.double(NA)]
dat[,q3.Can:=as.double(NA)]
dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),q3.Lof0711 := p.adjust(pLof0711, method='fdr')] 
dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),q3.Lof0714 := p.adjust(pLof0714, method='fdr')] 
dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),q3.Can := p.adjust(pCan, method='fdr')]

# For combined p-value
# after masking out unplaced
dat[,q2.comb0711Can:=as.double(NA)]
dat[,q2.comb0714Can:=as.double(NA)]
dat[!(CHROM %in% c('Unplaced')), q2.comb0711Can := p.adjust(p.comb0711Can, method='fdr')]
dat[!(CHROM %in% c('Unplaced')), q2.comb0714Can := p.adjust(p.comb0714Can, method='fdr')]

# For combined p-value
# after masking out inversions and unplaced
dat[,q3.comb0711Can:=as.double(NA)]
dat[,q3.comb0714Can:=as.double(NA)]
dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q3.comb0711Can := p.adjust(p.comb0711Can, method='fdr')]
dat[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q3.comb0714Can := p.adjust(p.comb0714Can, method='fdr')]


###################################################
# mark outliers based on fdr-corrected q-values
# or low p-values in both populations
###################################################
# qthresh <- 0.6
# 
# # in each pop
# # with inversions
# #dat[,outlierLof0714_q2 := 0]
# #dat[q2.Lof0714<qthresh, outlierLof0714_q2 := 1]
# #	dat[,sum(outlierLof0714_q2, na.rm=TRUE)]
# 
# dat[,outlierLof071114_q2 := 0]
# dat[q2.Lof071114<qthresh, outlierLof071114_q2 := 1]
# 	dat[,sum(outlierLof071114_q2, na.rm=TRUE)]
# 
# dat[,outlierCan_q2 := 0]
# dat[q2.Can<qthresh, outlierCan_q2 := 1]
# 	dat[,sum(outlierCan_q2, na.rm=TRUE)]
# 
# # in each pop
# # without inversions
# #dat[,outlierLof0714_q3 := 0]
# #dat[q3.Lof0714<qthresh, outlierLof0714_q3 := 1]
# #	dat[,sum(outlierLof0714_q3, na.rm=TRUE)]
# 
# dat[,outlierLof071114_q3 := 0]
# dat[q3.Lof071114<qthresh, outlierLof071114_q3 := 1]
# 	dat[,sum(outlierLof071114_q3, na.rm=TRUE)]
# 
# 	dat[, sum(q3.Lof071114<0.3, na.rm=TRUE)] # higher FDR threshold
# 
# 
# dat[,outlierCan_q3 := 0]
# dat[q3.Can<qthresh, outlierCan_q3 := 1]
# 	dat[,sum(outlierCan_q3, na.rm=TRUE)]
# 
# 	dat[,sum(outlierLof071114_q3 & outlierCan_q3, na.rm=TRUE)] # Lof and Can?
# 
# 
# # combined Lof and Can
# # with inversions
# #dat[,outlierLof0714_Can_q2 := 0]
# #dat[q2.comb0714Can<qthresh, outlierLof0714_Can_q2 := 1]
# #	dat[,sum(outlierLof0714_Can_q2, na.rm=TRUE)]
# 
# dat[,outlierLof071114_Can_q2 := 0]
# dat[q2.comb071114Can<qthresh, outlierLof071114_Can_q2 := 1]
# 	dat[,sum(outlierLof071114_Can_q2, na.rm=TRUE)]
# 
# # combined Lof and Can
# # without inversions
# #dat[,outlierLof0714_Can_q3 := 0]
# #dat[q3.comb0714Can<qthresh, outlierLof0714_Can_q3 := 1]
# #	dat[,sum(outlierLof0714_Can_q3, na.rm=TRUE)]
# 
# dat[,outlierLof071114_Can_q3 := 0]
# dat[q3.comb071114Can<qthresh, outlierLof071114_Can_q3 := 1]
# 	dat[,sum(outlierLof071114_Can_q3, na.rm=TRUE)]
# 
# 	dat[,sum(outlierLof071114_q3 & outlierLof071114_Can_q3, na.rm=TRUE)] # Lof and combined?
# 	dat[,sum(q3.Lof071114<0.3 & outlierLof071114_Can_q3, na.rm=TRUE)] # Lof<0.03 and combined?
# 	dat[,sum(outlierCan_q3 & outlierLof071114_Can_q3, na.rm=TRUE)] # Can and combined?
# 	
# 
# # loci with low p-values in both populations
# # with inversions
# #dat[,outlierLof0714andCan_p2 := 0]
# #dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('Unplaced')) & pLof0714<0.001 & pCan<0.001, outlierLof0714andCan_p2 := 1]
# #	dat[,sum(outlierLof0714andCan_p2, na.rm=TRUE)]
# 
# dat[,outlierLof071114andCan_p2 := 0]
# dat[kmer25==1 & dpFlag==1 & densFlag==1 & !(CHROM %in% c('Unplaced')) & pLof071114<0.001 & pCan<0.001, outlierLof071114andCan_p2 := 1]
# 	dat[,sum(outlierLof071114andCan_p2, na.rm=TRUE)]
# 	
# 
# # loci with low p-values in both populations
# # without inversions
# #dat[,outlierLof0714andCan_p3 := 0]
# #dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pLof0714<0.001 & pCan<0.001, outlierLof0714andCan_p3 := 1]
# #	dat[,sum(outlierLof0714andCan_p3, na.rm=TRUE)]
# 
# dat[,outlierLof071114andCan_p3 := 0]
# dat[kmer25==1 & dpFlag==1 & densFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pLof071114<0.001 & pCan<0.001, outlierLof071114andCan_p3 := 1]
# 	dat[,sum(outlierLof071114andCan_p3, na.rm=TRUE)]
# 
# 	dat[,sum(outlierLof071114_Can_q3 & outlierLof071114andCan_p3, na.rm=TRUE)]
# 

#############
# write out
#############
# write out
outfile <- paste('analysis/wfs_nullmodel_outliers_07-11-14_Can.csv.gz', sep='')
outfile
write.table(dat, file=gzfile(outfile), sep=',', row.names=FALSE, quote=FALSE)


#write out for Udi: just the ranked p-values, dropping Unplaced, inversions
# have to write out each population separately so that locus trimming is done appropriately	
# NEED TO REDO WITH 07-11 AND 07-14
outfile1 <- paste('analysis/wfs_nullmodel_outliers_Udi_Lof.tsv.gz', sep='')
outfile1
out1 <- dat[!(CHROM %in% c("Unplaced", "LG01", 'LG02', 'LG07', 'LG12')) & !is.na(q3.Lof0711), .(CHROM, POS, q3.Lof071114)]
nrow(out1)
write.table(out1, file=gzfile(outfile1), sep='\t', row.names=FALSE, quote=FALSE)

outfile2 <- paste('analysis/wfs_nullmodel_outliers_Udi_Can.tsv.gz', sep='')
outfile2
out2 <- dat[!(CHROM %in% c("Unplaced", "LG01", 'LG02', 'LG07', 'LG12')) & !is.na(q3.Can), .(CHROM, POS, q3.Can)]
nrow(out2)
write.table(out2, file=gzfile(outfile2), sep='\t', row.names=FALSE, quote=FALSE)

outfile3 <- paste('analysis/wfs_nullmodel_outliers_Udi_comb.tsv.gz', sep='')
outfile3
out3 <- dat[!(CHROM %in% c("Unplaced", "LG01", 'LG02', 'LG07', 'LG12')) & !is.na(q3.comb071114Can), .(CHROM, POS, q3.comb071114Can)]
nrow(out3)
write.table(out3, file=gzfile(outfile3), sep='\t', row.names=FALSE, quote=FALSE)

outfile4 <- paste('analysis/wfs_nullmodel_outliers_Udi_union.tsv.gz', sep='')
outfile4
out4 <- dat[!(CHROM %in% c("Unplaced", "LG01", 'LG02', 'LG07', 'LG12')) & (!is.na(q3.Lof0711) | !is.na(q3.Lof0714) | !is.na(q3.Can) | !is.na(q3.comb071114Can)), .(CHROM, POS, q3.union = pmin(q3.Lof071114, q3.Can, q3.comb071114Can, na.rm=TRUE))]
nrow(out4)
write.table(out4, file=gzfile(outfile4), sep='\t', row.names=FALSE, quote=FALSE)

#########################
# examine
#########################
## number of loci
nrow(dat) # total
dat[,.(sum(!is.na(Freq_11)), sum(!is.na(Freq_11))/.N)] # Lof 07-11
dat[,.(sum(!is.na(Freq_14)), sum(!is.na(Freq_14))/.N)] # Lof 07-14
dat[,.(sum(!is.na(Freq_Can40)), sum(!is.na(Freq_Can40))/.N)] # Can
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_14)), sum(!is.na(Freq_11) & !is.na(Freq_14))/.N)] # both
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_Can40)), sum(!is.na(Freq_11) & !is.na(Freq_Can40))/.N)] # both
dat[,.(sum(!is.na(Freq_14) & !is.na(Freq_Can40)), sum(!is.na(Freq_14) & !is.na(Freq_Can40))/.N)] # both
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40)), sum(!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40))/.N)] # all

	# not in unplaced
dat[,.(sum(!is.na(Freq_11) & !(CHROM %in% c('Unplaced'))), round(sum(!is.na(Freq_11) & !(CHROM %in% c('Unplaced')))/.N, 2))] # Lof 07-11
dat[,.(sum(!is.na(Freq_14) & !(CHROM %in% c('Unplaced'))), round(sum(!is.na(Freq_14) & !(CHROM %in% c('Unplaced')))/.N, 2))] # Lof 07-14
dat[,.(sum(!is.na(Freq_Can40) & !(CHROM %in% c('Unplaced'))), round(sum(!is.na(Freq_Can40) & !(CHROM %in% c('Unplaced')))/.N,2))] # Can
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_14) & !(CHROM %in% c('Unplaced'))), round(sum(!is.na(Freq_11) & !is.na(Freq_14) & !(CHROM %in% c('Unplaced')))/.N,2))] # both
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_Can40) & !(CHROM %in% c('Unplaced'))), round(sum(!is.na(Freq_11) & !is.na(Freq_Can40) & !(CHROM %in% c('Unplaced')))/.N, 2))] # both
dat[,.(sum(!is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% c('Unplaced'))), round(sum(!is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% c('Unplaced')))/.N,2))] # both
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% c('Unplaced'))), round(sum(!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% c('Unplaced')))/.N,2))] # all

	# not in unplaced or inversion LGs
locs <- c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')
dat[,.(sum(!is.na(Freq_11) & !(CHROM %in% locs)), round(sum(!is.na(Freq_11) & !(CHROM %in% locs))/.N, 2))] # Lof 07-11
dat[,.(sum(!is.na(Freq_14) & !(CHROM %in% locs)), round(sum(!is.na(Freq_14) & !(CHROM %in% locs))/.N, 2))] # Lof 07-14
dat[,.(sum(!is.na(Freq_Can40) & !(CHROM %in% locs)), round(sum(!is.na(Freq_Can40) & !(CHROM %in% locs))/.N,2))] # Can
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_14) & !(CHROM %in% locs)), round(sum(!is.na(Freq_11) & !is.na(Freq_14) & !(CHROM %in% locs))/.N,2))] # both
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_Can40) & !(CHROM %in% locs)), round(sum(!is.na(Freq_11) & !is.na(Freq_Can40) & !(CHROM %in% locs))/.N, 2))] # both
dat[,.(sum(!is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% locs)), round(sum(!is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% locs))/.N,2))] # both
dat[,.(sum(!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% locs)), round(sum(!is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !(CHROM %in% locs))/.N,2))] # all

## number of outliers: each population separately
	# without unplaced
dat[,sum(q2.Lof0711<0.2, na.rm=TRUE)]
	dat[q2.Lof0711<0.2, table(CHROM)]
	dat[q2.Lof0711<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof0711, q2.Lof0711)]

dat[,sum(q2.Lof0714<0.2, na.rm=TRUE)]
	dat[q2.Lof0714<0.2, table(CHROM)]
	dat[q2.Lof0714<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof0714, q2.Lof0714)]

dat[,sum(q2.Can<0.2, na.rm=TRUE)]
	dat[q2.Can<0.2, table(CHROM)]
	dat[q2.Can<0.2, .(CHROM, POS, Freq_Can40, Freq_CanMod, pCan, q2.Can)]

		# overlap
	dat[,sum(q2.Lof0711<0.2 & q2.Lof0714<0.2, na.rm=TRUE)]
	dat[,sum(q2.Lof0711<0.2 & q2.Can<0.2, na.rm=TRUE)]
	dat[,sum(q2.Lof0714<0.2 & q2.Can<0.2, na.rm=TRUE)]

	# without inversions or unplaced
dat[,sum(q3.Lof0711<0.2, na.rm=TRUE)]
	dat[q3.Lof0711<0.2, table(CHROM)]
	dat[q3.Lof0711<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof0711, q3.Lof0711)]

dat[,sum(q3.Lof0714<0.2, na.rm=TRUE)]
	dat[q3.Lof0714<0.2, table(CHROM)]
	dat[q3.Lof0714<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof0714, q3.Lof0714)]

dat[,sum(q3.Can<0.2, na.rm=TRUE)]
	dat[q3.Can<0.2, table(CHROM)]
	dat[q3.Can<0.2, .(CHROM, POS, Freq_Can40, Freq_CanMod, pCan, q3.Can)]

		# overlap
	dat[,sum(q3.Lof0711<0.2 & q3.Lof0714<0.2, na.rm=TRUE)]
	dat[,sum(q3.Lof0711<0.2 & q3.Can<0.2, na.rm=TRUE)]
	dat[,sum(q3.Lof0714<0.2 & q3.Can<0.2, na.rm=TRUE)]


# number of outliers: combined p-values
	# with inversions
dat[, sum(q2.comb0711Can<0.2, na.rm=TRUE)]
	dat[q2.comb0711Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0711, pCan)]

dat[, sum(q2.comb0714Can<0.2, na.rm=TRUE)]
	dat[q2.comb0714Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0714, pCan)]

		# overlap
	dat[,sum(q2.comb0711Can<0.2 & q2.comb0714Can<0.2, na.rm=TRUE)]

	# without inversions
dat[, sum(q3.comb0711Can<0.2, na.rm=TRUE)]
	dat[q3.comb0711Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0711, pCan)]

dat[, sum(q3.comb0714Can<0.2, na.rm=TRUE)]
	dat[q3.comb0714Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0714, pCan)]

		# overlap
	dat[,sum(q3.comb0711Can<0.2 & q3.Lof0711<0.2, na.rm=TRUE)]
	dat[,sum(q3.comb0711Can<0.2 & q3.Lof0714<0.2, na.rm=TRUE)]
	dat[,sum(q3.comb0711Can<0.2 & q3.Can<0.2, na.rm=TRUE)]

	dat[,sum(q3.comb0711Can<0.2 & q3.comb0714Can<0.2, na.rm=TRUE)]
	dat[,sum(q3.comb0714Can<0.2 & q3.Lof0711<0.2, na.rm=TRUE)]
	dat[,sum(q3.comb0714Can<0.2 & q3.Lof0714<0.2, na.rm=TRUE)]
	dat[,sum(q3.comb0714Can<0.2 & q3.Can<0.2, na.rm=TRUE)]


####################################################################
## do the same SNPs appear as outliers in both? p value approach
####################################################################

	# Lof 1907-2014 and Can
	# with inversions
# ntot <- dat[!(CHROM %in% c('Unplaced')) & !is.na(q2.Can) & !is.na(q2.Lof0714),.N]
# nLof <- dat[!(CHROM %in% c('Unplaced')) & !is.na(pCan),sum(q2.Lof0714<0.2, na.rm=TRUE)]
# nCan <- dat[!(CHROM %in% c('Unplaced')) & !is.na(q2.Lof0714),sum(q2.Can<0.2, na.rm=TRUE)]
# ncomb <- dat[!(CHROM %in% c('Unplaced')),sum(q2.Can<0.2 & q2.Lof0714<0.2, na.rm=TRUE)]
# nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
# nLof
# nCan
# ncomb # number in both
# nexp
# 
# dat[!(CHROM %in% c('Unplaced')) & q2.Can<0.2 & q2.Lof0714<0.2,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof0714)]
# 
# binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
# binom.test(x=5, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately


	# Lof 1907-2014 and Can
	# without inversions
# ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan) & !is.na(pLof0714),.N]
# nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan),sum(pLof0714<1e-3)]
# nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pLof0714),sum(pCan<1e-3)]
# ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(pCan<1e-3 & pLof0714<1e-3)]
# nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
# nLof
# nCan
# ncomb # number in both
# nexp
# 
# dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pCan<1e-3 & pLof0714<1e-3,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof0714)]
# 
# binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
# binom.test(x=4, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately

	# Lof 1907-2011-2014 and Can
	# with inversions
#ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & densFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pCan) & !is.na(pLof071114), .N]
#nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & densFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pCan), sum(pLof071114<1e-3)]
#nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & densFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pLof071114), sum(pCan<1e-3)]
#ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & densFlag==1 & !(CHROM %in% c('Unplaced')), sum(pCan<1e-3 & pLof071114<1e-3)]
#nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
#nLof
#nCan
#ncomb # number in both
#nexp
#
#dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & densFlag==1 & !(CHROM %in% c('Unplaced')) & pCan<1e-3 & pLof071114<1e-3,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof071114)]
#
#binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
#binom.test(x=8, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately

	# Lof 1907-2011-2014 and Can
	# without inversions
# pthresh <- 1e-2
# ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & densFlag==1 & binomp>0.05 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan) & !is.na(pLof071114),.N]
# nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & densFlag==1 & binomp>0.05 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan),sum(pLof071114<pthresh)]
# nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & densFlag==1 & binomp>0.05 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pLof071114),sum(pCan<pthresh)]
# ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & densFlag==1 & binomp>0.05 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(pCan<1e-2 & pLof071114<pthresh)]
# nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
# nLof
# nCan
# ncomb # number in both
# nexp
# 
# dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & densFlag==1 & binomp>0.05 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pCan<pthresh & pLof071114<pthresh,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof071114)] # look at loci
# dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & densFlag==1 & binomp>0.05 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pCan<pthresh & pLof071114<pthresh,table(CHROM)] # look at number of LGs
# 
# binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
# binom.test(x=7, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately



###########################
# other examination code
###########################

# examine combined p-values in various classes
dat[,range(p.comb, na.rm=TRUE)]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb, na.rm=TRUE)]

hist(dat[,-log10(p.comb)], breaks=40, col='grey')

	require(RColorBrewer)
	cols <- brewer.pal(3, 'Set1')
	bks <- seq(0,1,by=0.01)
	hist14_25k <- dat14_25k[,hist(ABS_DIFF_0714, breaks=bks, plot=FALSE)]
	hist11_25k <- dat11_25k[,hist(ABS_DIFF_0711, breaks=bks, plot=FALSE)]
	hist11_14_25k <- dat11_14_25k[,hist(ABS_DIFF_1114, breaks=bks, plot=FALSE)]

	quartz(width=9, height=3)
	# pdf(width=9, height=3, file='analysis/figures/abs_diff_hist_25k&150k.pdf')
	par(mfrow=c(1,3), las=1, cex.axis=0.8, tcl=-0.2, mgp=c(2.5,0.5, 0), mai=c(0.5, 0.5, 0.2, 0.1))

	plot(hist11_25k$mids, hist11_25k$density/sum(hist11_25k$density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='Proportion of loci', log='y', main='1907-2011', ylim=c(1e-6, 1e-1), xlim=c(0,0.6))
	lines(hist11_150k$mids, hist11_150k$density/sum(hist11_150k$density), type='o', col=cols[2], cex=0.5, main='1907-2011')

	plot(hist14_25k$mids, hist14_25k$density/sum(hist14_25k$density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='', log='y', main='1907-2014', ylim=c(1e-6, 1e-1), xlim=c(0,0.6))
	lines(hist14_150k$mids, hist14_150k$density/sum(hist14_150k$density), type='o', col=cols[2], cex=0.5)

	plot(hist11_14_25k$mids, hist11_14_25k$density/sum(hist11_14_25k$density), type='o', col=cols[1], cex=0.5, xlab='Frequency change', ylab='', log='y', main='2011-2014', ylim=c(1e-6, 1e-1), xlim=c(0,0.6))
	lines(hist11_14_150k$mids, hist11_14_150k$density/sum(hist11_14_150k$density), type='o', col=cols[2], cex=0.5)

	legend('topright', legend=c('25 kmer', '150 kmer'), lwd=1, pch=1, col=cols[1:2], cex=0.8)


# examine adjusted combined p-values
dat[,range(p.comb.adj, na.rm=TRUE)]
dat[,range(p.comb.adj2, na.rm=TRUE)]

dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb.adj, na.rm=TRUE)]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb.adj2, na.rm=TRUE)]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),range(p.comb.adj3, na.rm=TRUE)]

hist(dat[,-log10(p.comb.adj)], breaks=40, col='grey')
hist(dat[,-log10(p.comb.adj2)], breaks=40, col='grey')

	# zoom
hist(dat[,-log10(p.comb.adj)], breaks=200, col='grey', xlim=c(2,6), ylim=c(0,50))
hist(dat[,-log10(p.comb.adj2)], breaks=200, col='grey', xlim=c(2,6), ylim=c(0,50))



# plot outliers based on allele frequency change in both populations
par(mfrow=c(1,2))
dat[kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),][sample(.N,10000), plot(abs((Freq_11 +Freq_14)/2 - Freq_07), abs(Freq_CanMod - Freq_Can40), main='all that pass filters\n(sample of 5000)', xlim=c(0,0.5), ylim=c(0,0.5), col=rgb(0,0,0,0.2), cex=0.5)]
dat[p.comb.adj3<0.2 & kmer25==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), plot(abs((Freq_11 +Freq_14)/2 - Freq_07), abs(Freq_CanMod - Freq_Can40), main='outliers p.comb.adj2<0.05', xlim=c(0,0.5), ylim=c(0,0.5))]
