## Calculate correlation in allele frequency change across datasets
## Add some test comments

# load functions
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){ # not cod node
	require(data.table)
	require(plyr)
	require(parallel)
	ncores=2
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	ncores=20
}

# read in data
dat14_25k <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(dat14_25k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_0714')) # for 1907 vs. 2014
dat11_25k <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof11_25k.txt', header=TRUE); setnames(dat11_25k, 3:7, c('N_CHR_07', 'Freq_07', 'N_CHR_11', 'Freq_11', 'ABS_DIFF_0711')) # for 1907 vs. 2011
dat1114_25k <- fread('data/data_29.06.17/Frequency_table_Lof11_Lof14_25k.txt', header=TRUE); setnames(dat1114_25k, 3:7, c('N_CHR_11', 'Freq_11', 'N_CHR_14', 'Freq_14', 'ABS_DIFF_1114')) # for 2011 vs. 2014
#datCAN <- fread('analysis/Frequency_table_Can_40_Can_TGA.txt', header=TRUE); setnames(datCAN, 3:7, c('N_CHR_40', 'Freq_40', 'N_CHR_TGA', 'Freq_TGA', 'ABS_DIFF_40TGA')) # for 1907 vs. 2011


# make a nucleotide position for the whole genome
chrmax <- dat14_25k[,.(len=max(POS)), by=CHROM]
#chrmax <- dat14orig[,.(len=max(POS)), by=CHROM]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, CHROM)

setkey(dat14_25k, CHROM)
dat14_25k <- dat14_25k[chrmax[,.(CHROM, start)], ]
dat14_25k[,POSgen:=POS+start]

setkey(dat11_25k, CHROM)
dat11_25k <- dat11_25k[chrmax[,.(CHROM, start)], ]
dat11_25k[,POSgen:=POS+start]

setkey(dat1114_25k, CHROM)
dat1114_25k <- dat1114_25k[chrmax[,.(CHROM, start)], ]
dat1114_25k[,POSgen:=POS+start]

#setkey(datCAN, CHROM)
#datCAN <- datCAN[chrmax[,.(CHROM, start)], ]
#datCAN[,POSgen:=POS+start]

#setkey(dat14orig, CHROM)
#dat14orig <- dat14orig[chrmax[,.(CHROM, start)], ]
#dat14orig[,POSgen:=POS+start]


# combine the datasets to look at correlations across them
setkey(dat11_25k, CHROM, POS)
setkey(dat14_25k, CHROM, POS)
setkey(dat1114_25k, CHROM, POS)
#setkey(datCAN, CHROM, POS)

dat07.11_07.14 <- dat11_25k[dat14_25k, .(CHROM, POS, POSgen, N_CHR_07, N_CHR_11, N_CHR_14, Freq_07, Freq_11, Freq_14, ABS_DIFF_0711, ABS_DIFF_0714)]
dat07.14_11.14 <- dat14_25k[dat1114_25k, .(CHROM, POS, POSgen, N_CHR_07, N_CHR_11, N_CHR_14, Freq_07, Freq_11, Freq_14, ABS_DIFF_0714, ABS_DIFF_1114)]
dat07.11_11.14 <- dat11_25k[dat1114_25k, .(CHROM, POS, POSgen, N_CHR_07, N_CHR_11, N_CHR_14, Freq_07, Freq_11, Freq_14, ABS_DIFF_0711, ABS_DIFF_1114)]

#datCAN_14 <- datCAN[dat14, .(CHROM, POS, POSgen, N_CHR_40, N_CHR_TGA, N_CHR_07, N_CHR_14, Freq_40, Freq_TGA, Freq_07, Freq_14, ABS_DIFF_40TGA, ABS_DIFF_0714)]

#datCAN_11_14 <- datCAN_14[dat11, .(CHROM, POS, POSgen, N_CHR_40, N_CHR_TGA, N_CHR_07, N_CHR_11, N_CHR_14, Freq_40, Freq_TGA, Freq_07, Freq_11, Freq_14, ABS_DIFF_40TGA, ABS_DIFF_0711, ABS_DIFF_0714)]


###################################
# plot correlations among regions
###################################
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