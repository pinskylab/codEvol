## Calculate correlation in allele frequency and frequency change across datasets
## New comments!

# load functions
require(data.table)

# read in data on all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
dat <- fread("gunzip -c analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz")
	dat


###################################
# correlations in allele frequency among datasets
###################################
# old to new within a site
dat[!is.na(Freq_07) & !is.na(Freq_11),cor.test(Freq_07, Freq_11, method='pearson')]
dat[!is.na(Freq_07) & !is.na(Freq_14),cor.test(Freq_07, Freq_14, method='pearson')]
dat[!is.na(Freq_11) & !is.na(Freq_14),cor.test(Freq_11, Freq_14, method='pearson')]
dat[!is.na(Freq_Can40) & !is.na(Freq_CanMod),cor.test(Freq_Can40, Freq_CanMod, method='pearson')]

# across sites
dat[!is.na(Freq_07) & !is.na(Freq_Can40),cor.test(Freq_07, Freq_Can40, method='pearson')]
dat[!is.na(Freq_11) & !is.na(Freq_CanMod),cor.test(Freq_11, Freq_CanMod, method='pearson')]
dat[!is.na(Freq_14) & !is.na(Freq_CanMod),cor.test(Freq_14, Freq_CanMod, method='pearson')]


###############################################################
# plot correlations in allele frequency change among regions
###############################################################

	# 07-11 and 07-14
dat07.11_07.14[,plot(ABS_DIFF_0711, ABS_DIFF_0714, col=rgb(0,0,0,0.5), pch=16, cex=0.3)]
	dat07.11_07.14[,cor.test(ABS_DIFF_0711, ABS_DIFF_0714)]

	# 07-14 and 11-14
dat07.14_11.14[,plot(ABS_DIFF_0714, ABS_DIFF_1114, col=rgb(0,0,0,0.5), pch=16, cex=0.3)]
	dat07.14_11.14[,cor.test(ABS_DIFF_0714, ABS_DIFF_1114)]

	# 07-11 and 11-14
dat07.11_11.14[,plot(ABS_DIFF_0711, ABS_DIFF_1114, col=rgb(0,0,0,0.5), pch=16, cex=0.3)]
	dat07.11_11.14[,cor.test(ABS_DIFF_0711, ABS_DIFF_1114)]

	# CAN and 07-14
#datCAN_14[,plot(ABS_DIFF_40TGA, ABS_DIFF_0714, col=rgb(0,0,0,0.5), pch=16, cex=0.3)]
#	datCAN_14[,cor.test(ABS_DIFF_40TGA, ABS_DIFF_0714)]
#
#	datCAN_14[ABS_DIFF_40TGA>0.3 & ABS_DIFF_0714>0.3,]
#	
#	# all 3 comparisons
#	datCAN_11_14[ABS_DIFF_40TGA>0.3 & ABS_DIFF_0711>0.3 & ABS_DIFF_0714>0.3,]

