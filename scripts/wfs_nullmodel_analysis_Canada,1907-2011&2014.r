# Compare Canada, 1907-2011, 1907-2014, and 2011-2014 outlier loci from null model analysis (25k)
# to be run on my mac

require(data.table)

# read in data
	# Lof 1907-2011-2014 and 1907-2014
locnms <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
locnms2 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms2, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	setkey(locnms, CHROM, POS, N_CHR_07, Freq_07)
	setkey(locnms2, CHROM, POS, N_CHR_07, Freq_07)
	locnms <- locnms[locnms2,]
load('analysis/wfs_nullmodel_pvals_07-11-14.rdata') # dat. has the p-values from 1907-2011-2014
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		datLof <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_07, Freq_11, Freq_14)], by='locusnum')
		setnames(datLof, c('p', 'p.adj', 'cnt1', 'cnt2', 'cnt3'), c('pLof071114', 'p.adjLof071114', 'cnt07', 'cnt11', 'cnt14'))

		# Lof 1907-2014
	load('analysis/wfs_nullmodel_pvals_07-14.rdata') # dat. has the p-values from 1907-2014

			# make sure it matches
			nrow(dat)
			nrow(datLof)

			# merge
			dat <- as.data.table(dat)
			datLof <- merge(dat[,.(locusnum, p, p.adj)], datLof, by='locusnum')
			setnames(datLof, c('p', 'p.adj'), c('pLof0714', 'p.adjLof0714'))

	# Canada
locnms <- fread('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_Can40', 'Freq_Can40', 'N_CHR_CanMod', 'Freq_CanMod', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
load('analysis/wfs_nullmodel_pvals_Can.rdata') # dat
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging

		# make sure it matches
		nrow(dat)
		nrow(locnms)

		# merge
		dat <- as.data.table(dat)
		datCan <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_Can40, Freq_CanMod)], by='locusnum')
		setnames(datCan, c('p', 'p.adj', 'cnt1', 'cnt2'), c('pCan', 'p.adjCan', 'cntCan40', 'cntCanMod'))

##################
# merge Lof and Can
##################

setkey(datLof, CHROM, POS)
setkey(datCan, CHROM, POS)
dat <- merge(datCan, datLof, all=TRUE)[, .(CHROM, POS, cnt07, cnt11, cnt14, cntCan40, cntCanMod, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0714, pLof071114, pCan)]

	# make sure the merge worked
	nrow(dat)
	nrow(datLof)
	dat[,.(sum(!is.na(Freq_07)), sum(!is.na(Freq_11)), sum(!is.na(Freq_14)))] # should match nrow(datLof)
	nrow(datCan)
	dat[,.(sum(!is.na(Freq_Can40)), sum(!is.na(Freq_CanMod)))] # should match nrow(datCan)
	dat[,sum(!is.na(Freq_07) & !is.na(Freq_Can40))] # number of loci genotyped in both
	
################
# add kmer25
################
kmer='25' # flag for file names

kmer25Lof <- fread('data_2017.11.24/Norway_25K_mer_positions.txt') # loci that pass ther kmer25 mapping filter
kmer25Can <- fread('data_2017.11.24/Canada_25K_mer_positions.txt')

setkey(kmer25Lof, CHROM, POS)
setkey(kmer25Can, CHROM, POS)
kmer25 <- merge(kmer25Lof, kmer25Can, all=TRUE) # merge and keep all
	nrow(kmer25Lof)
	nrow(kmer25Can)
	nrow(kmer25)

kmer25[,kmer25 := 1] # add a flag for loci that pass kmer25 filter
dat <- merge(dat, kmer25, by=c('CHROM', 'POS'), all.x=TRUE)
dat[is.na(kmer25), kmer25:=0] # set NAs to 0

	# make sure the merge worked
	nrow(dat)
	dat[,sum(kmer25)] == nrow(kmer25) # should match
	dat[!is.na(Freq_07), sum(kmer25)] == nrow(kmer25Lof)
	dat[!is.na(Freq_Can40), sum(kmer25)] == nrow(kmer25Can)
	dat[!is.na(Freq_07) & !is.na(Freq_Can40),sum(kmer25)] # number of loci genotyped in both populations that pass kmer filter

#########################
# add depth statistic
#########################
dpLof <- fread('analysis/Outlier_sequencing_depth_statistic_All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz.txt')
dpCan <- fread('analysis/Outlier_sequencing_depth_statistic_All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz.txt')
	setnames(dpLof, c('CHROM', 'POS', 'dpstatLof'))
	setnames(dpCan, c('CHROM', 'POS', 'dpstatCan'))
depthstat <- merge(dpLof, dpCan, all=TRUE)
	nrow(dpLof)
	nrow(dpCan)
	nrow(depthstat)

dat <- merge(dat, depthstat, by=c('CHROM', 'POS'), all.x=TRUE)
	nrow(dat)

dat[,dpLofFlag := dpstatLof < quantile(dpstatLof, na.rm=TRUE, probs=0.95)] # TRUE if not a depth outlier (outlier are top 5%)
	dat[is.na(dpLofFlag), dpLofFlag := FALSE] # FALSE if locus not genotyped
dat[,dpCanFlag := dpstatCan < quantile(dpstatCan, na.rm=TRUE, probs=0.95)]
	dat[is.na(dpCanFlag), dpCanFlag := FALSE] 

	dat[,sum(dpLofFlag)]
	dat[,sum(dpLofFlag)]/dat[,sum(!is.na(dpstatLof))] # should be 0.95
	dat[,sum(dpCanFlag)]
	dat[,sum(dpCanFlag)]/dat[,sum(!is.na(dpstatCan))]

dat[,dpFlag := dpLofFlag & dpCanFlag]
	dat[,sum(dpFlag)] # number of loci genotyped in both that pass kmer and depth filter in both

###############################################
# combine p-values for individual loci
# https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/
###############################################

dat[,p.comb0714Can := pchisq(-2 * sum(log(c(pLof0714, pCan))),df=2*2,lower=FALSE), by=1:nrow(dat)] # Lof 07-14 and Can
dat[,p.comb071114Can := pchisq(-2 * sum(log(c(pLof071114, pCan))),df=2*2,lower=FALSE), by=1:nrow(dat)] # Lof 07-11-14 and Can

##################################
# FDR-correct p-values
##################################

# FDR-correct each population separately
# after masking out unplaced and loci failing kmer and depth
dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('Unplaced')),q2.Lof0714 := p.adjust(pLof0714, method='fdr')] 
dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('Unplaced')),q2.Lof071114 := p.adjust(pLof071114, method='fdr')] 
dat[kmer25==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')),q2.Can := p.adjust(pCan, method='fdr')]

# FDR-correct each population separately
# after masking out inversions and unplaced,  and loci failing kmer and depth
dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),q3.Lof0714 := p.adjust(pLof0714, method='fdr')] 
dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),q3.Lof071114 := p.adjust(pLof071114, method='fdr')] 
dat[kmer25==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),q3.Can := p.adjust(pCan, method='fdr')]

# For combined p-value
# no masking
dat[,q.comb0714Can := p.adjust(p.comb0714Can, method='fdr')]
dat[,q.comb071114Can := p.adjust(p.comb071114Can, method='fdr')]

# For combined p-value
# after masking out unplaced and loci failing kmer and depth
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('Unplaced')),q2.comb0714Can := p.adjust(p.comb0714Can, method='fdr')]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('Unplaced')),q2.comb071114Can := p.adjust(p.comb071114Can, method='fdr')]

# For combined p-value
# after masking out inversions and unplaced,  and loci failing kmer and depth
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')), q3.comb0714Can := p.adjust(p.comb0714Can, method='fdr')]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),q3.comb071114Can := p.adjust(p.comb071114Can, method='fdr')]


###################################################
# mark outliers based on fdr-corrected q-values
# or low p-values in both populations
###################################################

# in each pop
# with inversions
dat[,outlierLof0714_q2 := 0]
dat[q2.Lof0714<0.2, outlierLof0714_q2 := 1]
	dat[,sum(outlierLof0714_q2, na.rm=TRUE)]

dat[,outlierLof071114_q2 := 0]
dat[q2.Lof071114<0.2, outlierLof071114_q2 := 1]
	dat[,sum(outlierLof071114_q2, na.rm=TRUE)]

dat[,outlierCan_q2 := 0]
dat[q2.Can<0.2, outlierCan_q2 := 1]
	dat[,sum(outlierCan_q2, na.rm=TRUE)]

# in each pop
# without inversions
dat[,outlierLof0714_q3 := 0]
dat[q3.Lof0714<0.2, outlierLof0714_q3 := 1]
	dat[,sum(outlierLof0714_q3, na.rm=TRUE)]

dat[,outlierLof071114_q3 := 0]
dat[q3.Lof071114<0.2, outlierLof071114_q3 := 1]
	dat[,sum(outlierLof071114_q3, na.rm=TRUE)]

dat[,outlierCan_q3 := 0]
dat[q3.Can<0.2, outlierCan_q3 := 1]
	dat[,sum(outlierCan_q3, na.rm=TRUE)]


# combined Lof and Can
# with inversions
dat[,outlierLof0714_Can_q2 := 0]
dat[q2.comb0714Can<0.2, outlierLof0714_Can_q2 := 1]
	dat[,sum(outlierLof0714_Can_q2, na.rm=TRUE)]

dat[,outlierLof071114_Can_q2 := 0]
dat[q2.comb071114Can<0.2, outlierLof071114_Can_q2 := 1]
	dat[,sum(outlierLof071114_Can_q2, na.rm=TRUE)]

# combined Lof and Can
# without inversions
dat[,outlierLof0714_Can_q3 := 0]
dat[q3.comb0714Can<0.2, outlierLof0714_Can_q3 := 1]
	dat[,sum(outlierLof0714_Can_q3, na.rm=TRUE)]

dat[,outlierLof071114_Can_q3 := 0]
dat[q3.comb071114Can<0.2, outlierLof071114_Can_q3 := 1]
	dat[,sum(outlierLof071114_Can_q3, na.rm=TRUE)]



# loci with low p-values in both populations
# with inversions
dat[,outlierLof0714andCan_p2 := 0]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('Unplaced')) & pLof0714<0.001 & pCan<0.001, outlierLof0714andCan_p2 := 1]
	dat[,sum(outlierLof0714andCan_p2, na.rm=TRUE)]

dat[,outlierLof071114andCan_p2 := 0]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('Unplaced')) & pLof071114<0.001 & pCan<0.001, outlierLof071114andCan_p2 := 1]
	dat[,sum(outlierLof071114andCan_p2, na.rm=TRUE)]
	

# loci with low p-values in both populations
# without inversions
dat[,outlierLof0714andCan_p3 := 0]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pLof0714<0.001 & pCan<0.001, outlierLof0714andCan_p3 := 1]
	dat[,sum(outlierLof0714andCan_p3, na.rm=TRUE)]

dat[,outlierLof071114andCan_p3 := 0]
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pLof071114<0.001 & pCan<0.001, outlierLof071114andCan_p3 := 1]
	dat[,sum(outlierLof071114andCan_p3, na.rm=TRUE)]


#############
# write out
#############
# write out tab-separated so Bastiaan can read easily
outfile <- paste('analysis/wfs_nullmodel_outliers_07-11-14_Can_', kmer, 'k.tsv.gz', sep='')
outfile
write.table(dat, file=gzfile(outfile), sep='\t', row.names=FALSE, quote=FALSE)


#write out for Udi: just the ranked p-values, dropping Unplaced and SNPs that fail kmer25 filter and SNPs that fail depth filter in both populations
# have to write out each population separately so that locus trimming is done appropriately	
outfile1 <- paste('analysis/wfs_nullmodel_outliers_Udi_Lof.tsv.gz', sep='')
outfile1
out1 <- dat[kmer25==1 & dpLofFlag==TRUE & CHROM != "Unplaced" & !is.na(pLof0714), .(CHROM, POS, pLof0714, pLof071114)]
nrow(out1)
write.table(out1, file=gzfile(outfile1), sep='\t', row.names=FALSE, quote=FALSE)

outfile2 <- paste('analysis/wfs_nullmodel_outliers_Udi_Can.tsv.gz', sep='')
outfile2
out2 <- dat[kmer25==1 & dpCanFlag==TRUE & CHROM != "Unplaced" & !is.na(pCan), .(CHROM, POS, pCan)]
nrow(out2)
write.table(out2, file=gzfile(outfile2), sep='\t', row.names=FALSE, quote=FALSE)

outfile3 <- paste('analysis/wfs_nullmodel_outliers_Udi_comb.tsv.gz', sep='')
outfile3
out3 <- dat[kmer25==1 & dpCanFlag==TRUE & dpLofFlag==TRUE & CHROM != "Unplaced" & !is.na(p.comb071114Can), .(CHROM, POS, p.comb0714Can, p.comb071114Can)]
nrow(out3)
write.table(out3, file=gzfile(outfile3), sep='\t', row.names=FALSE, quote=FALSE)

#########################
# examine
#########################
## number of loci
nrow(dat) # total
dat[,sum(!is.na(cnt07))] # Lof
dat[,sum(!is.na(cntCan40))] # Can
dat[,sum(!is.na(cnt07) & !is.na(cntCan40))] # both

	# that pass kmer filter
nrow(dat[kmer25==1,])
dat[kmer25==1,sum(!is.na(cnt07))] # Lof
dat[kmer25==1,sum(!is.na(cntCan40))] # Can
dat[kmer25==1,sum(!is.na(cnt07) & !is.na(cntCan40))] # both

	# that pass kmer and depth filters
dat[kmer25==1 & dpLofFlag==1,sum(!is.na(cnt07))] # Lof
dat[kmer25==1 & dpCanFlag==1,sum(!is.na(cntCan40))] # Can
dat[kmer25==1 & dpFlag==1,sum(!is.na(cnt07) & !is.na(cntCan40))] # both

	# that pass kmer and depth filters, not in unplaced
dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('Unplaced')),sum(!is.na(cnt07))] # Lof
dat[kmer25==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')),sum(!is.na(cntCan40))] # Can
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('Unplaced')),sum(!is.na(cnt07) & !is.na(cntCan40))] # both

	# that pass kmer and depth filters, not in unplaced or inversion LGs
dat[kmer25==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(!is.na(cnt07))] # Lof
dat[kmer25==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(!is.na(cntCan40))] # Cam
dat[kmer25==1 & dpFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(!is.na(cnt07) & !is.na(cntCan40))] # both


## number of outliers: each population separately
	# with inversions
dat[,sum(outlierLof0714_q2)]
	dat[outlierLof0714_q2==TRUE, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof0714, q2.Lof0714)]

dat[,sum(outlierLof071114_q2)]
	dat[outlierLof071114_q2==TRUE, table(CHROM)]
	dat[outlierLof071114_q2==TRUE, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof071114, q2.Lof071114)]

dat[,sum(outlierCan_q2)]
	dat[outlierCan_q2==TRUE, table(CHROM)]
	dat[outlierCan_q2==TRUE, .(CHROM, POS, Freq_Can40, Freq_CanMod, pCan, q2.Can)]

	# without inversions
dat[,sum(outlierLof0714_q3)]
	dat[outlierLof0714_q3==TRUE, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof0714, q3.Lof0714)]

dat[,sum(outlierLof071114_q3)]
	dat[outlierLof071114_q3==TRUE, table(CHROM)]
	dat[outlierLof071114_q3==TRUE, .(CHROM, POS, Freq_07, Freq_11, Freq_14, pLof071114, q3.Lof071114)]

dat[,sum(outlierCan_q3)]
	dat[outlierCan_q3==TRUE, table(CHROM)]
	dat[outlierCan_q3==TRUE, .(CHROM, POS, Freq_Can40, Freq_CanMod, pCan, q3.Can)]


# number of outliers: combined p-values
	# with inversions
dat[, sum(q2.comb0714Can<0.2, na.rm=TRUE)]
	dat[q2.comb0714Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0714, pCan, outlierLof0714_Can_q2, outlierLof071114_Can_q2)]

dat[, sum(q2.comb071114Can<0.2, na.rm=TRUE)]
	dat[q2.comb071114Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0714, pCan, outlierLof0714_Can_q2, outlierLof071114_Can_q2)]

	# without inversions
dat[, sum(q3.comb0714Can<0.2, na.rm=TRUE)]
	dat[q3.comb0714Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0714, pCan, outlierLof071114_Can_q3)]

dat[, sum(q3.comb071114Can<0.2, na.rm=TRUE)]
	dat[q3.comb071114Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pLof0714, pCan, outlierLof0714_Can_q3)]




####################################################################
## do the same SNPs appear as outliers in both? p value approach
####################################################################

	# Lof 1907-2014 and Can
	# with inversions
ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pCan) & !is.na(pLof0714),.N]
nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pCan),sum(pLof0714<1e-3)]
nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pLof0714),sum(pCan<1e-3)]
ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('Unplaced')),sum(pCan<1e-3 & pLof0714<1e-3)]
nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
nLof
nCan
ncomb # number in both
nexp

dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('Unplaced')) & pCan<1e-3 & pLof0714<1e-3,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof0714)]

binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
binom.test(x=5, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately


	# Lof 1907-2014 and Can
	# without inversions
ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan) & !is.na(pLof0714),.N]
nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan),sum(pLof0714<1e-3)]
nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pLof0714),sum(pCan<1e-3)]
ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(pCan<1e-3 & pLof0714<1e-3)]
nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
nLof
nCan
ncomb # number in both
nexp

dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pCan<1e-3 & pLof0714<1e-3,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof0714)]

binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
binom.test(x=4, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately

	# Lof 1907-2011-2014 and Can
	# with inversions
ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pCan) & !is.na(pLof071114),.N]
nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pCan),sum(pLof071114<1e-3)]
nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('Unplaced')) & !is.na(pLof071114),sum(pCan<1e-3)]
ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('Unplaced')),sum(pCan<1e-3 & pLof071114<1e-3)]
nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
nLof
nCan
ncomb # number in both
nexp

dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('Unplaced')) & pCan<1e-3 & pLof071114<1e-3,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof071114)]

binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
binom.test(x=6, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately

	# Lof 1907-2011-2014 and Can
	# without inversions
ntot <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan) & !is.na(pLof071114),.N]
nLof <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pCan),sum(pLof071114<1e-3)]
nCan <- dat[kmer25==1 & dpLofFlag==1 & dpCanFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & !is.na(pLof071114),sum(pCan<1e-3)]
ncomb <- dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')),sum(pCan<1e-3 & pLof071114<1e-3)]
nexp <- nLof/ntot * nCan/ntot * ntot # expected number in both
nLof
nCan
ncomb # number in both
nexp

dat[kmer25==1 & dpCanFlag==1 & dpLofFlag==1 & !(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced')) & pCan<1e-3 & pLof071114<1e-3,.(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, pCan, pLof071114)]

binom.test(x=ncomb, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations
binom.test(x=5, n=ntot, p=nLof/ntot * nCan/ntot) # statistical test relative to expectations: each cluster separately



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
