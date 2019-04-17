# compare initial and final allele frequencies near outlier loci and near non-outlier loci

######################
# calculate initial and final allele frequencies
######################
#require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node
require(data.table) # not on cod node

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
#outl <- fread("zcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz") # on cod node
outl <- fread("gzcat analysis/wfs_nullmodel_outliers_07-11-14_Can.csv.gz") # on mac
#	outl

# remove inversions and unplaced
outl <- outl[!(CHROM %in% c('LG01', 'LG02', 'LG07', 'LG12', 'Unplaced'))]

nrow(outl) # 661091

# how many outlier comparisons?
outl[,table(q3.Lof0711<0.2)] # 11
outl[,table(q3.Lof0714<0.2)] # 6
outl[,table(q3.Can<0.2)] # 13
outl[,table(q3.comb0711Can<0.2)] # 8
outl[,table(q3.comb0714Can<0.2)] # 10
outl[,table(q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2)] # 23 Norway union
outl[,table(q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2)] # 26 Canada union

# set up 1000 not-outlier loci
	outl[,nearoutlier:=FALSE] # to mark loci near outliers in any population
	outl[q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2, nearoutlier:=TRUE] # first mark loci that are outliers
	for(i in outl[,which(q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2)]) { # add loci that are <= 5000 bp from outliers
		outl[CHROM==outl$CHROM[i] & abs(POS - outl$POS[i]) <= 5000, nearoutlier:=TRUE]
	}

	# non-outliers for Norway
	allloci <- outl[nearoutlier==FALSE & !is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_14), .(CHROM, POS)] # list of all loci that aren't near outliers
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
	allloci <- outl[nearoutlier==FALSE & !is.na(Freq_Can40) & !is.na(Freq_CanMod), .(CHROM, POS)]
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
	allloci <- outl[nearoutlier==FALSE & !is.na(Freq_07) & !is.na(Freq_11) & !is.na(Freq_14) & !is.na(Freq_Can40) & !is.na(Freq_CanMod), .(CHROM, POS)]
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


# write out
out <- outl[q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2 | notoutlierLof==TRUE | notoutlierCAN==TRUE, .(CHROM, POS, q3.Lof0711, q3.Lof0714, q3.Can, q3.comb0711Can, q3.comb0714Can, notoutlierLof, notoutlierCAN, Freq_07low, Freq_Can40low, Freq_11low, Freq_14low, Freq_CanModlow)]
write.csv(out, file='analysis/freq_init_final_outliers.csv')

###############
# plot
###############
require(RColorBrewer)
require(data.table)

out <- fread('analysis/freq_init_final_outliers.csv')

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
out[q3.Lof0711<0.2 | q3.Lof0714<0.2,plot(density(Freq_07low, from=0, to=1), ylim=ylims, xlim=c(0,1), type='l', xlab='Initial allele frequency', ylab='Density', main='', log='', col=cols[1], lwd=2)]
out[q3.Can<0.2,lines(density(Freq_Can40low, from=0, to=1), col=cols[2], lwd=2)]
		# outlier combined pop
out[q3.comb0711Can<0.2  | q3.comb0714Can<0.2,lines(density(Freq_07low, from=0, to=1), col=cols[1], lwd=2, lty=2)]
out[q3.comb0711Can<0.2  | q3.comb0714Can<0.2,lines(density(Freq_Can40low, from=0, to=1, adjust=1), col=cols[2], lwd=2, lty=2)]
		# notoutlier
out[notoutlierLof==TRUE,lines(density(Freq_07low, from=0, to=1), col=cols[1], lty=3)]
out[notoutlierCAN==TRUE,lines(density(Freq_Can40low, from=0, to=1), col=cols[2], lty=3)]

legend('topright', legend=c('LOF_07', 'CAN40', 'Outlier single pop', 'Outlier combined', 'Not outlier'), col=c(cols[c(1,2)], 'black', 'black', 'black'), lty=c(rep(1,3), 2, 3), bty='n', cex=0.5)

	# Final
		# outlier single pop
out[q3.Lof0711<0.2 | q3.Lof0714<0.2,plot(density(Freq_14low, from=0, to=1), ylim=ylims, xlim=c(0,1), type='l', xlab='Final allele frequency', ylab='Density', main='', log='', col=cols[1], lwd=2)]
out[q3.Lof0711<0.2 | q3.Lof0714<0.2,lines(density(Freq_11low, from=0, to=1), col=cols[3], lwd=2)]
out[q3.Can<0.2,lines(density(Freq_CanModlow, from=0, to=1), col=cols[2], lwd=2)]
		# outlier combined pop
out[q3.comb0711Can<0.2 | q3.comb0714Can<0.2,lines(density(Freq_14low, from=0, to=1), col=cols[1], lwd=2, lty=2)]
out[q3.comb0711Can<0.2 | q3.comb0714Can<0.2,lines(density(Freq_11low, from=0, to=1), col=cols[3], lwd=2, lty=2)]
out[q3.comb0711Can<0.2 | q3.comb0714Can<0.2,lines(density(Freq_CanModlow, from=0, to=1, adjust=1), col=cols[2], lwd=2, lty=2)]
		# notoutlier
out[notoutlierLof==TRUE,lines(density(Freq_14low, from=0, to=1), col=cols[1], lty=3)]
out[notoutlierLof==TRUE,lines(density(Freq_11low, from=0, to=1), col=cols[3], lty=3)]
out[notoutlierCAN==TRUE,lines(density(Freq_CanModlow, from=0, to=1), col=cols[2], lty=3)]

legend('topright', legend=c('LOF_S_14', 'LOF_S_11', 'CANMod', 'Outlier single pop', 'Outlier combined', 'Not outlier'), col=c(cols[c(1,3,2)], 'black', 'black', 'black'), lty=c(rep(1,4), 2, 3), bty='n', cex=0.5)

dev.off()



### plot density of initial and final frequencies: Union of single-population and combined outliers
quartz(width=6, height=3)
# pdf(width=6, height=3, file='figures/freq_init_final_outliers_union.pdf')
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Initial
		# outlier
out[q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2,plot(density(Freq_07low, from=0, to=1), ylim=ylims, xlim=c(0,1), type='l', xlab='Initial allele frequency', ylab='Density', main='', log='', col=cols[1], lwd=2)]
out[q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2,lines(density(Freq_Can40low, from=0, to=1), col=cols[2], lwd=2)]
		# notoutlier
out[notoutlierLof==TRUE,lines(density(Freq_07low, from=0, to=1), col=cols[1], lty=3)]
out[notoutlierCAN==TRUE,lines(density(Freq_Can40low, from=0, to=1), col=cols[2], lty=3)]

legend('topright', legend=c('LOF_07', 'CAN40', 'Outlier', 'Not outlier'), col=c(cols[c(1,2)], 'black', 'black'), lty=c(rep(1,3), 3), bty='n', cex=0.5)

	# Final
		# outlier
out[q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2,plot(density(Freq_14low, from=0, to=1), ylim=ylims, xlim=c(0,1), type='l', xlab='Final allele frequency', ylab='Density', main='', log='', col=cols[1], lwd=2)]
out[q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2,lines(density(Freq_11low, from=0, to=1), col=cols[3], lwd=2)]
out[q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2,lines(density(Freq_CanModlow, from=0, to=1), col=cols[2], lwd=2)]
		# notoutlier
out[notoutlierLof==TRUE,lines(density(Freq_14low, from=0, to=1), col=cols[1], lty=3)]
out[notoutlierLof==TRUE,lines(density(Freq_11low, from=0, to=1), col=cols[3], lty=3)]
out[notoutlierCAN==TRUE,lines(density(Freq_CanModlow, from=0, to=1), col=cols[2], lty=3)]

legend('topright', legend=c('LOF_S_14', 'LOF_S_11', 'CANMod', 'Outlier', 'Not outlier'), col=c(cols[c(1,3,2)], 'black', 'black'), lty=c(rep(1,4), 3), bty='n', cex=0.5)

dev.off()


### orca plot of initial and final frequencies
out[,freqinitmean:=mean(c(Freq_07low, Freq_Can40low), na.rm=TRUE), by = seq_len(nrow(out))]
setkey(out, freqinitmean)
n <- out[q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2,.N]
locinds <- out[, which(q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2)]

quartz(width=6, height=3)
# pdf(width=6, height=3, file='figures/freq_init_final_outliers_union_bylocus.pdf')
par(mfrow=c(1,2), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.7, las=1, mgp=c(1.5, 0.3, 0), tcl=-0.15)

	# Canada
plot(-1,-1, xlim=c(0,1), ylim=c(0,n), xlab='Frequency', ylab='Loci', main='Canada')
for(i in 1:length(locinds)){
	thislty <- 2
	if(out[locinds[i], (q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2) & !is.na(q3.Can)]) thislty <- 1
	if(out[locinds[i], (q3.Can>=0.2 & q3.comb0711Can>=0.2 & q3.comb0714Can>=0.2) | is.na(q3.Can)]) thislty <- 3

	out[locinds[i], lines(x=c(Freq_Can40low, Freq_CanModlow), y=c(i,i), type='l', lty=thislty)]
	out[locinds[i], points(x=Freq_Can40low, y=i, type='p', cex=0.5, col=cols[1])]
	out[locinds[i], points(x=Freq_CanModlow, y=i, type='p', cex=0.5, col=cols[3])]
}
legend('topleft', legend=c('CAN40', 'CANMod', 'outlier', 'not outlier'), col=c(cols[c(1,3)], 'black', 'black'), pch=c(1,1,NA,NA), lty=c(NA, NA, 1,3), bty='n', cex=0.5)

	# Norway
plot(-1,-1, xlim=c(0,1), ylim=c(0,n), xlab='Frequency', ylab='Loci', main='Norway')
for(i in 1:length(locinds)){
	thislty <- 2
	if(out[locinds[i], (q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2) & (!is.na(q3.Lof0711) | !is.na(q3.Lof0714))]) thislty <- 1
	if(out[locinds[i], (q3.Lof0711>=0.2 & q3.Lof0714>=0.2 & q3.comb0711Can>=0.2 & q3.comb0714Can>=0.2) | (is.na(q3.Lof0711) & is.na(q3.Lof0714))]) thislty <- 3

	out[locinds[i], lines(x=c(Freq_07low, Freq_11low), y=c(i,i), type='l', lty=thislty)]
	out[locinds[i], lines(x=c(Freq_07low, Freq_14low), y=c(i,i), type='l', lty=thislty)]
	out[locinds[i], points(x=Freq_07low, y=i, type='p', cex=0.5, col=cols[1])]
	out[locinds[i], points(x=Freq_11low, y=i, type='p', cex=0.5, col=cols[2])]
	out[locinds[i], points(x=Freq_14low, y=i, type='p', cex=0.5, col=cols[3])]
	
}

legend('topleft', legend=c('LOF_S_07', 'LOF_S_11', 'LOF_S_14', 'outlier', 'not outlier'), col=c(cols[1:3], 'black', 'black'), pch=c(1,1,1,NA,NA), lty=c(NA,NA,NA,1,3), bty='n', cex=0.5)

dev.off()