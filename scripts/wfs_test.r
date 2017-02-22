# test wfs simulations

require(foreach)
require(doParallel)
source('scripts/wfs_simp.r')

# make cluster
cl <- makeCluster(3)
registerDoParallel(cl)


# probability of fixation of a neutral allele should equal the initial frequency
	ne <- 100
	gen=1000 # run for a long time: likely to fixaton (expected time is 4N)
	smin <- smax <- 0
	fs <- seq(1/100, 91/100, by=5/100) # initial allele frequencies
	nsims <- 1000


	sims <- data.frame(f1=rep(fs, rep(nsims, length(fs))))
	sims$f2 <- sims$f1samp <- sims$f2samp <- NA
	for(j in 1:length(fs)){
		print(j)
		thesesims <- foreach(i=1:nsims, .combine=rbind) %dopar% {
			wfs_simp(i, f1=fs[j], s=0, c1=50, c2=50, gen=gen, ne=ne)
		}
		sims[((j-1)*1000+1):(j*1000),c('f2', 'f1samp', 'f2samp')] <- as.data.frame(thesesims)
	}

	propfix <- aggregate(list(prop=sims$f2), by=list(f1=sims$f1), FUN=function(x) sum(x==1)/nsims)

	# plot of proportion fixed vs. initial frequency
	quartz(width=5, height=4)
	# png(width=5, height=4, units='in', res=300, file='analysis/figures/wfs_sim_test_s=0.png')
	par(las=1)
	with(propfix, plot(f1, prop, xlab='Initial allele frequency', ylab='Proportion fixed', main='Ne=100, s=0, nsims=1000'))
	abline(0,1)
	
	dev.off()




# probability of fixation of a non-neutral allele should equal s if initial frequency is 1/ne
# (usually written as 2s, but that's when fitness of waa is 1+s2. here it is 1+s)
# for large N and small s
	ne <- 100
	gen=1000 # run for a long time: likely to fixation in 4N generations with drift
	ss <- seq(0.05,1,by=0.05)
	f1 <- 1/ne # starting allele freq
	nsims <- 1000

	sims2 <- data.frame(s=rep(ss, rep(nsims, length(ss))))
	sims2$f2 <- sims2$f1samp <- sims2$f2samp <- NA
	length(ss)
	for(j in 1:length(ss)){
		print(j)
		thesesims <- foreach(i=1:nsims, .combine=rbind) %dopar% {
			if(i %% 100 == 0) cat(i)
			wfs_simp(i, f1=f1, s=ss[j], c1=50, c2=50, gen=gen, ne=ne)
		}
		sims2[((j-1)*1000+1):(j*1000),c('f2', 'f1samp', 'f2samp')] <- as.data.frame(thesesims)
	}

	sum(sims2$f2 %in% c(0,1))/nsims/length(ss) # fraction that reached fixation: make gen larger until this=1
	propfix2 <- aggregate(list(prop=sims2$f2), by=list(s=sims2$s), FUN=function(x) sum(x==1)/nsims)
	propfix2$s_exp1 <- (1-exp(-propfix2$s))/(1-exp(-2*ne*propfix2$s)) # Eq. 2.21 in Graur & Li
	propfix2$s_exp2 <- propfix2$s/(1-exp(-2*ne*propfix2$s)) # Eq. 2.22 in Graur & Li

	# plot of proportion fixed vs. s
	quartz(width=5, height=4)
	# png(width=5, height=4, units='in', res=300, file='analysis/figures/wfs_sim_test_s>0.png')
	par(las=1)
	with(propfix2, plot(s, prop, xlab='Selection coefficient (s)', ylab='Proportion fixed at p=1', main='Ne=100, p0=1/Ne, nsims=1000', ylim=c(0,0.7)))
	with(propfix2, lines(s, s_exp1, col='green'))
	legend('topleft', col='green', lty=1, legend='(1-exp(-s))/(1-exp(-2*ne*s))', bty='n', cex=0.7)
	
	dev.off()



## examine how much drift vs. sampling variance for large ne
	ne <- 10000 # like NEA cod
	gen=11 # like NEA cod
	smin <- smax <- 0
	f1min <- f1max <- ne/2/ne # starting allele freq
	nsims <- 1000


	sims <- foreach(i=1:nsims, .combine=rbind) %dopar% {
		if(i %% 100 == 0) cat(i)
		wfs(i, f1min=f1min, f1max=f1max, smin=smin, smax=smax, c1=50, c2=50, gen=gen, ne=ne)
	}
	sims <- as.data.frame(sims)

	# examine drift variance and sampling variance
	par(mfrow=c(3,1))
	bks <- seq(0,1,length.out=20)
	hist(sims$f1samp, main='f1samp', breaks=bks, col='grey'); abline(v=f1min, col='red', lwd=2)
	hist(sims$f2, main='f2', breaks=bks, col='grey'); abline(v=f1min, col='red', lwd=2)
	hist(sims$f2samp-sims$f2, main='delta f2samp', breaks=seq(-0.5,0.5,length.out=20), col='grey'); abline(v=0, col='red', lwd=2)


## examine how much drift to expect for large Ne
	# observations
	dat <- locnms <- fread('analysis/Frequency_table_Lof07_Lof14.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star


	# simulations
	ne <- 10000
	gen=11 # run for a long time: likely to fixaton (expected time is 4N)
	fs <- seq(1/100, 91/100, by=5/100) # initial allele frequencies
	nsims <- 100000 # takes 20 min or so

	sims <- data.frame(f1=rep(fs, rep(nsims, length(fs))))
	sims$f2 <- sims$f1samp <- sims$f2samp <- NA
	for(j in 1:length(fs)){
		print(j)
		thesesims <- foreach(i=1:nsims, .combine=rbind) %dopar% {
			wfs_simp(i, f1=fs[j], s=0, c1=44, c2=44, gen=gen, ne=ne)
		}
		sims[((j-1)*nsims+1):(j*nsims),c('f2', 'f1samp', 'f2samp')] <- as.data.frame(thesesims)
	}
	
	xlims <- range(c(dat[,Freq_2-Freq_1], sims$f2samp-sims$f1samp))
	cols <- c('black', '#7fcdbb', '#edf8b1', '#2c7fb8') # sims, obs 44, low obs, high obs

	# plot histogram of diff vs. initial observed frequency
	quartz(width=7, height=6)
	# png(width=7, height=7, filename='analysis/figures/wfs_sim_test_hist_obsdiff_by_obsf1_with_data.png', units='in', res=300)
	f1s <- seq(0,1.1,by=0.1)
	bks <- seq(-1,1,by=0.05)
	par(mfrow=c(3,4), las=1, tcl=-0.3, mgp=c(2.4,0.6,0), mai=c(0.5, 0.5, 0.3, 0.05))
	for(i in 2:length(f1s)){
		inds <- (sims$f1samp >= f1s[i-1]) & (sims$f1samp < f1s[i])
		hst <- hist(sims$f2samp[inds] - sims$f1samp[inds], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j <- which(hst$density>0)
		plot(hst$mids[j], hst$density[j], type='o', xlab='Sample ∆freq', ylab='log(density)', main=paste('initial ', f1s[i-1], '-', f1s[i], sep=''), pch=16, xlim=xlims, log='y', ylim=c(0.0001,20), col=cols[1])
		abline(v=c(-0.06, 0.06), col='grey', lty=3) # genome-wide allele frequency change

			# low sample size
		inds3 <- dat[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1<44 | N_CHR_2<44)]
		hst3 <- hist(dat[inds3,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j3 <- which(hst3$density>0)
		lines(hst3$mids[j3], hst3$density[j3], type='o', pch=16, col=cols[3])

			# high sample size
		inds4 <- dat[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & (N_CHR_1>44 & N_CHR_2>44)]
		hst4 <- hist(dat[inds4,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j4 <- which(hst4$density>0)
		lines(hst4$mids[j4], hst4$density[j4], type='o', pch=16, col=cols[4])

			# sample size 44
		inds2 <- dat[,Freq_1 >= f1s[i-1] & Freq_1 < f1s[i] & N_CHR_1==44 & N_CHR_2==44]
		hst2 <- hist(dat[inds2,Freq_2 - Freq_1], plot=FALSE, breaks=bks) # returns some 0 values that get eliminated, but appear as holes in the histogram: doesn't change our interpretation
		j2 <- which(hst2$density>0)
		lines(hst2$mids[j2], hst2$density[j2], type='o', pch=16, col=cols[2])

	}
	
	legend('bottomright', col=cols, lty=1, legend=c('Sims n=44', 'Obs n=44', 'Obs low n', 'Obs high n'), cex=0.7, bty='n')

	dev.off()


#######################
# clean up cluster
#######################
stopCluster(cl)

