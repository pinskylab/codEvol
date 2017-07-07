# compare WFS null model simulations (s=0) to observations
## examine how much drift to expect for large Ne
# to run on my mac

require(data.table)
require(bit)
require(ff)
require(RColorBrewer)

########################################################################
# plot histogram of diff vs. initial observed frequency (25k or 150k)
########################################################################
	# observations
	dat <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); kmer='25k'; yrs='07-14' # the name and observed frequencies of all the loci, from output by Bastiaan Star
	dat <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE); setnames(dat, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')); kmer='150k'; yrs='07-14' # the name and observed frequencies of all the loci, from output by Bastiaan Star
	dat[,summary(N_CHR_1)] # 44
	dat[,summary(N_CHR_2)] # 44

	# simulations
	ffload('analysis/temp/wfs_simsnull_ff44,44', overwrite=TRUE) # loads thisout.ff. quick on my mac


	xlims <- range(c(dat[,Freq_2-Freq_1], thisout.ff['f2samp',]-thisout.ff['f1samp',]))
	cols <- c('black', '#7fcdbb', '#edf8b1', '#2c7fb8') # sims, obs 44, low obs, high obs

	quartz(width=7, height=6)
	# png(width=7, height=7, filename=paste('analysis/figures/wfs_nullmodel_hist_obsdiff_by_obsf1_with_data', yrs, ',', kmer, '.png', sep=''), units='in', res=300)
	f1s <- c(seq(0,1,by=0.1), 1.01)
	bks <- seq(-1,1,by=0.05)
	par(mfrow=c(3,4), las=1, tcl=-0.3, mgp=c(2.4,0.6,0), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.8)
	for(i in 2:length(f1s)){
		inds <- (thisout.ff['f1samp',] >= f1s[i-1]) & (thisout.ff['f1samp',] < f1s[i])
		hst <- hist(thisout.ff['f2samp',inds] - thisout.ff['f1samp',inds], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j <- which(hst$density>0)
		plot(hst$mids[j], hst$density[j], type='o', xlab='Sample ∆freq', ylab='log(density)', main=paste('initial ', f1s[i-1], '-', f1s[i], sep=''), pch=16, xlim=xlims, log='y', ylim=c(0.0001,20), col=cols[1])
		abline(v=c(-0.06, 0.06), col='grey', lty=3) # genome-wide allele frequency change

			# low sample size
#		inds3 <- dat[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1<44 | N_CHR_2<44)]
#		hst3 <- hist(dat[inds3,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
#		j3 <- which(hst3$density>0)
#		lines(hst3$mids[j3], hst3$density[j3], type='o', pch=16, col=cols[3])

			# high sample size
#		inds4 <- dat[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1>44 & N_CHR_2>44)]
#		hst4 <- hist(dat[inds4,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
#		j4 <- which(hst4$density>0)
#		lines(hst4$mids[j4], hst4$density[j4], type='o', pch=16, col=cols[4])

			# sample size 44
		inds2 <- dat[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & N_CHR_1==44 & N_CHR_2==44]
		hst2 <- hist(dat[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[2])

	}
	
	legend('bottomright', col=cols[1:2], lty=1, legend=c('Sims n=44', 'Obs n=44'), cex=0.7, bty='n')
#	legend('bottomright', col=cols, lty=1, legend=c('Sims n=44', 'Obs n=44', 'Obs low n', 'Obs high n'), cex=0.7, bty='n')

	dev.off()


#########################################################################
# plot histogram of diff vs. initial observed frequency 
# 25k AND 150k
#########################################################################
	# load data
		# for 1907-2014
	dat25 <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); setnames(dat25, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	dat150 <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE); setnames(dat150, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	yrs='07-14'

		# for 1907-2011
	dat25 <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof11_25k.txt', header=TRUE); setnames(dat25, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	dat150 <- fread('data/data_29.06.17/Frequency_table_Lof07_Lof11_150k.txt', header=TRUE); setnames(dat150, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	yrs='07-11'


		# for 2011-2014
	dat25 <- fread('data/data_29.06.17/Frequency_table_Lof11_Lof14_25k.txt', header=TRUE); setnames(dat25, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	dat150 <- fread('data/data_29.06.17/Frequency_table_Lof11_Lof14_150k.txt', header=TRUE); setnames(dat150, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
	yrs='11-14'


	xlims <- range(c(dat25[,Freq_2-Freq_1], dat150[,Freq_2-Freq_1], thisout.ff['f2samp',]-thisout.ff['f1samp',]))
	cols <- c('black', brewer.pal(6, 'PRGn')[c(3:1,4:6)]) # sims, obs 25k n<44, n=44, n>44, obs150k n<44, n=44, n>44

	quartz(width=7, height=6)
	# png(width=7, height=7, filename=paste('analysis/figures/wfs_nullmodel_hist_obsdiff_by_obsf1_with_data', yrs, ',25k&150k.png', sep=''), units='in', res=300)
	f1s <- c(seq(0,1,by=0.1), 1.01)
	bks <- seq(-1,1,by=0.05)
	par(mfrow=c(3,4), las=1, tcl=-0.3, mgp=c(2.4,0.6,0), mai=c(0.5, 0.5, 0.3, 0.05), cex.axis=0.8)

	# add panel not divided by initial freq
		hst <- hist(thisout.ff['f2samp',] - thisout.ff['f1samp',], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j <- which(hst$density>0)
		plot(hst$mids[j], hst$density[j], type='o', xlab='Sample ∆freq', ylab='log(density)', main=paste('All loci', yrs), pch=16, xlim=xlims, log='y', ylim=c(0.00001,20), col=cols[1])
		abline(v=c(-0.06, 0.06), col='grey', lty=3) # genome-wide allele frequency change

			# 25kmer n<44
		inds2 <- dat25[, N_CHR_1<44 | N_CHR_2<44]
		hst2 <- hist(dat25[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[2])

			# 25kmer n=44
		inds2 <- dat25[, N_CHR_1==44 & N_CHR_2==44]
		hst2 <- hist(dat25[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[3])

			# 25kmer n>44
		inds2 <- dat25[, N_CHR_1>44 | N_CHR_2>44]
		hst2 <- hist(dat25[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[4])

			# 150kmer n<44
		inds2 <- dat150[, N_CHR_1<44 | N_CHR_2<44]
		hst2 <- hist(dat150[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[5])

			# 150kmer n=44
		inds2 <- dat150[, N_CHR_1==44 & N_CHR_2==44]
		hst2 <- hist(dat150[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[6])

			# 150kmer n>44
		inds2 <- dat150[, N_CHR_1>44 | N_CHR_2>44]
		hst2 <- hist(dat150[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[7])

	# panels for each bin of initial frequencies
	for(i in 2:length(f1s)){
		print(i)
		
		inds <- (thisout.ff['f1samp',] >= f1s[i-1]) & (thisout.ff['f1samp',] < f1s[i])
		hst <- hist(thisout.ff['f2samp',inds] - thisout.ff['f1samp',inds], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j <- which(hst$density>0)
		plot(hst$mids[j], hst$density[j], type='o', xlab='Sample ∆freq', ylab='log(density)', main=paste('initial ', f1s[i-1], '-', f1s[i], sep=''), pch=16, xlim=xlims, log='y', ylim=c(0.0001,20), col=cols[1])
		abline(v=c(-0.06, 0.06), col='grey', lty=3) # genome-wide allele frequency change

			# 25kmer n<44
		inds2 <- dat25[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1<44 | N_CHR_2<44)]
		hst2 <- hist(dat25[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[2])

			# 25kmer n=44
		inds2 <- dat25[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & N_CHR_1==44 & N_CHR_2==44]
		hst2 <- hist(dat25[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[3])

			# 25kmer n>44
		inds2 <- dat25[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1>44 | N_CHR_2>44)]
		hst2 <- hist(dat25[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[4])

			# 150kmer n<44
		inds2 <- dat150[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1<44 | N_CHR_2<44)]
		hst2 <- hist(dat150[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[5])

			# 150kmer n=44
		inds2 <- dat150[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & N_CHR_1==44 & N_CHR_2==44]
		hst2 <- hist(dat150[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[6])

			# 150kmer n>44
		inds2 <- dat150[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1>44 | N_CHR_2>44)]
		hst2 <- hist(dat150[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[7])

	}

	legend('bottomright', col=cols, lty=1, legend=c('Sims n=44', 'Obs 25k n<44', 'Obs 25k n=44', 'Obs 25k n>44', 'Obs 150k n<44', 'Obs 150k n=44', 'Obs 150k n>44'), cex=0.7, bty='n')
#	legend('bottomright', col=cols, lty=1, legend=c('Sims n=44', 'Obs n=44', 'Obs low n', 'Obs high n'), cex=0.7, bty='n')

	dev.off()

