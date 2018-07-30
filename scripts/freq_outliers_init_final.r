# compare initial and final allele frequencies near outlier loci and near non-outlier loci

#outliertype <- 'bypop' # use q3.Lof071114 and q3.Can < 0.2 to define outlier loci
outliertype <- 'combinedpop' # use q3.comb071114Can < 0.2


######################
# calculate initial and final allele frequencies
######################
#require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node
require(data.table) # not on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
#outl <- fread("zcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz") # on cod node 270MB file
outl <- fread("gzcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz") # on mac 108MB file
#	outl

# remove inversions and unplaced (also remove unmappable?)
outl <- outl[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]

# remove low-quality rows (SNPs that fail a filter)
# dont' filter on coverage since that is population-specific. apply that later (and already included in outlier determination column)
outl <- outl[kmer25==1,]

nrow(outl) # 415238

# how many outlier comparisons?
outl[,table(outlierLof071114_q3)] # 40
outl[,table(outlierCan_q3)] # 4
outl[,table(outlierLof071114_Can_q3)] # 23

# set up 1000 not-outlier loci
	outl[,nearoutlierLof:=FALSE]
	outl[,nearoutlierCAN:=FALSE]
	outl[,nearoutlierComb:=FALSE]

	# for Norway
	outl[outlierLof071114_q3==1, nearoutlierLof:=TRUE] # mark loci that are outliers
	for(i in outl[,which(outlierLof071114_q3==1)]) { # add loci that are <= 5000 bp from outliers
		outl[CHROM==outl$CHROM[i] & abs(POS - outl$POS[i]) <= 5000, nearoutlierLof:=TRUE]
	}
	allloci <- outl[nearoutlierLof==FALSE & !is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_14), .(CHROM, POS)] # list of all loci that aren't near outliers
	allloci$dist = NA
	allloci$dist[2:nrow(allloci)] <- allloci$POS[2:nrow(allloci)] - allloci$POS[1:(nrow(allloci)-1)] # calculate distance to next locus to the left
		dim(allloci)
	allloci <- allloci[dist > 5000 & POS > 5000,] # trim out loci <5000 bp from start of CHROM or from another locus (doesn't check end of CHROM...)
		dim(allloci)
	notoutl <- allloci[sample(.N, 1000)] # sample 1000 loci to be the "not outliers"
	notoutl[,notoutlierLof:=TRUE] # label these as the 1000 not outliers
	notoutl[,dist:=NULL] # remove dist column
	outl <- merge(outl, notoutl, by=c('CHROM', 'POS'), all.x=TRUE) # label the not-outliers
	outl[is.na(notoutlierLof), notoutlierLof:=FALSE]


	# for Canada
	outl[outlierCan_q3==1, nearoutlierCAN:=TRUE]
	for(i in outl[,which(outlierCan_q3==1)]) {
		outl[CHROM==outl$CHROM[i] & abs(POS - outl$POS[i]) <= 5000, nearoutlierCAN:=TRUE]
	}
	allloci <- outl[nearoutlierCAN==FALSE & !is.na(Freq_Can40) & !is.na(Freq_CanMod), .(CHROM, POS)]
	allloci$dist = NA
	allloci$dist[2:nrow(allloci)] <- allloci$POS[2:nrow(allloci)] - allloci$POS[1:(nrow(allloci)-1)]
		dim(allloci)
	allloci <- allloci[dist > 5000 & POS > 5000,]
		dim(allloci)
	notoutl <- allloci[sample(.N, 1000)]
	notoutl[,notoutlierCAN:=TRUE]
	notoutl[,dist:=NULL]
	outl <- merge(outl, notoutl, by=c('CHROM', 'POS'), all.x=TRUE)
	outl[is.na(notoutlierCAN), notoutlierCAN:=FALSE]

	# for Combined
	outl[outlierLof071114_Can_q3==1, nearoutlierCAN:=TRUE]
	for(i in outl[,which(outlierLof071114_Can_q3==1)]) {
		outl[CHROM==outl$CHROM[i] & abs(POS - outl$POS[i]) <= 5000, nearoutlierCAN:=TRUE]
	}
	allloci <- outl[nearoutlierComb==FALSE & !is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !is.na(Freq_CanMod), .(CHROM, POS)]
	allloci$dist = NA
	allloci$dist[2:nrow(allloci)] <- allloci$POS[2:nrow(allloci)] - allloci$POS[1:(nrow(allloci)-1)]
		dim(allloci)
	allloci <- allloci[dist > 5000 & POS > 5000,]
		dim(allloci)
	notoutl <- allloci[sample(.N, 1000)]
	notoutl[,notoutlierComb:=TRUE]
	notoutl[,dist:=NULL]
	outl <- merge(outl, notoutl, by=c('CHROM', 'POS'), all.x=TRUE)
	outl[is.na(notoutlierComb), notoutlierComb:=FALSE]

	outl[,summary(notoutlierLof)] # 1000
	outl[,summary(notoutlierCAN)]
	outl[,summary(notoutlierComb)]


# calculate initial frequency of the allele that increases in frequency (use 2014 as modern for Laf)
outl[Freq_07 < Freq_14, ':=' (Freq_07low=Freq_07, Freq_11low=Freq_11, Freq_14low=Freq_14)]
outl[Freq_07 >= Freq_14, ':=' (Freq_07low=1-Freq_07, Freq_11low=1-Freq_11, Freq_14low=1-Freq_14)]
outl[Freq_Can40 < Freq_CanMod, ':=' (Freq_Can40low=Freq_Can40, Freq_CanModlow=Freq_CanMod)]
outl[Freq_Can40 >= Freq_CanMod, ':=' (Freq_Can40low=1-Freq_Can40, Freq_CanModlow=1-Freq_CanMod)]


###############
# plot
###############
require(RColorBrewer)
require(data.table)


# plotting parameters
cols <- brewer.pal(5, 'Set1')
colstrans <- paste(cols, '55', sep='')
cex=0.5
lwd=2
ylims=c(0,12)

### plot density of initial and final frequencies
quartz(width=6, height=3)
# pdf(width=6, height=3, file='figures/freq_init_final_outliers.pdf')
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Initial
		# outlier single pop
outl[outlierLof071114_q3==1,plot(density(Freq_07low, from=0, to=1), ylim=ylims, xlim=c(0,1), type='l', xlab='Initial allele frequency', ylab='Density', main='', log='', col=cols[1], lwd=2)]
outl[outlierCan_q3==1,lines(density(Freq_Can40low, from=0, to=1, adjust=0.1), col=cols[2], lwd=2)]
		# outlier combined pop
outl[outlierLof071114_Can_q3==1,lines(density(Freq_07low, from=0, to=1), col=cols[1], lwd=2, lty=2)]
outl[outlierLof071114_Can_q3==1,lines(density(Freq_Can40low, from=0, to=1, adjust=1), col=cols[2], lwd=2, lty=2)]
		# notoutlier
outl[notoutlierLof==TRUE,lines(density(Freq_07low, from=0, to=1), col=cols[1], lty=3)]
outl[notoutlierCAN==TRUE,lines(density(Freq_Can40low, from=0, to=1), col=cols[2], lty=3)]

legend('topright', legend=c('LOF_07', 'CAN40', 'Outlier single pop', 'Outlier combined', 'Not outlier'), col=c(cols[c(1,2)], 'black', 'black', 'black'), lty=c(rep(1,3), 2, 3), bty='n', cex=0.5)

	# Final
		# outlier single pop
outl[outlierLof071114_q3==1,plot(density(Freq_14low, from=0, to=1), ylim=ylims, xlim=c(0,1), type='l', xlab='Final allele frequency', ylab='Density', main='', log='', col=cols[1], lwd=2)]
outl[outlierLof071114_q3==1,lines(density(Freq_11low, from=0, to=1), col=cols[3], lwd=2)]
outl[outlierCan_q3==1,lines(density(Freq_CanModlow, from=0, to=1, adjust=50), col=cols[2], lwd=2)]
		# outlier combined pop
outl[outlierLof071114_Can_q3==1,lines(density(Freq_14low, from=0, to=1), col=cols[1], lwd=2, lty=2)]
outl[outlierLof071114_Can_q3==1,lines(density(Freq_11low, from=0, to=1), col=cols[3], lwd=2, lty=2)]
outl[outlierLof071114_Can_q3==1,lines(density(Freq_CanModlow, from=0, to=1, adjust=1), col=cols[2], lwd=2, lty=2)]
		# notoutlier
outl[notoutlierLof==TRUE,lines(density(Freq_14low, from=0, to=1), col=cols[1], lty=3)]
outl[notoutlierLof==TRUE,lines(density(Freq_11low, from=0, to=1), col=cols[3], lty=3)]
outl[notoutlierCAN==TRUE,lines(density(Freq_CanModlow, from=0, to=1), col=cols[2], lty=3)]

legend('topright', legend=c('LOF_S_14', 'LOF_S_11', 'CANMod', 'Outlier single pop', 'Outlier combined', 'Not outlier'), col=c(cols[c(1,3,2)], 'black', 'black', 'black'), lty=c(rep(1,4), 2, 3), bty='n', cex=0.5)

dev.off()