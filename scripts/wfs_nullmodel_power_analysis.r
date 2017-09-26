# Calculate probability of the null model analysis detecting selection when it exists
# run after wfs_make_sims_null_power.r, wfs_nullmodel_function.r and wfs_nullmodel_combine.r
# Script analyzes the output

#############################
## Read in and summarize Lof
#############################
require(data.table)

# Read in data
	# Lof-like simulations
locnms <- fread('analysis/Frequency_table_PowerSims_Lof_Ne46000_cnt46_44.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from simulations
load('analysis/wfs_nullmodel_pvals_Power_Lof_Ne46000_cnt46_44.rdata') # dat
suffix <- '_Lof_Ne46000_cnt46_44'

# Then continue here
	# make sure the rows match
	nrow(dat)
	nrow(locnms)


# merge
locnms[,locusnum:=1:nrow(locnms)] # add a locusnumber for plotting and merging
dat <- as.data.table(dat)
dat <- merge(dat, locnms[,.(locusnum, CHROM, POS, Freq_1, Freq_2, ABS_DIFF)], by='locusnum')

# re-order columns
setcolorder(dat, c('CHROM', 'POS', 'locusnum', 'n', 'cnt1', 'cnt2', 'Freq_1', 'Freq_2', 'ABS_DIFF', 'p', 'pmax', 'p.adj'))

# extract s and f1 values from CHROM field
dat[,s := as.numeric(gsub('s=', '', unlist(strsplit(CHROM,split=','))[seq(1,2*nrow(dat),by=2)]))]
dat[,f1 := as.numeric(gsub('f1=', '', unlist(strsplit(CHROM,split=','))[seq(2,2*nrow(dat),by=2)]))]

# summarize fraction p<0.05, 0.005, 0.0005 for each combination of s and f1
power <- dat[,.(pow05=sum(p<0.05)/.N, pow005=sum(p<0.005)/.N, pow0005=sum(p<0.0005)/.N),by=.(s,f1)]


# f1 vs. s with pow
	require(RColorBrewer)
#	colrmp <- colorRamp(brewer.pal(9, 'OrRd'))
	colrmp <- colorRamp(terrain.colors(9))
	colrmpfun <- function(x){
		rgbs <- colrmp(x)
		return(rgb(rgbs, maxColorValue=256))
	}

	quartz(width=6,height=6)
	layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,1))
	power[, plot(s, f1, type='p', cex=2, pch=15, xlab='s', ylab='f1', col=colrmpfun(pow), main=paste('null model power analysis\n', suffix), xlim=c(0,1), ylim=c(0,1))]

	# add a legend
	par(mai=c(1,0,1,0))
	plot(rep(0,11), seq(0,0.8,length.out=11), col=colrmpfun(seq(0,1,length.out=11)), pch=15, cex=2, type='p', xaxt='n', yaxt='n', bty='n', xlim=c(0,3), ylim=c(0,0.9), xlab='', ylab='') 
	text(rep(1,6), seq(0,0.8,length.out=6), labels=seq(0,1,length.out=6), cex=0.8)
	text(-0.05, 0.85, 'power', cex=0.8, pos=4)
		
	dev.off()

# f1 vs. s with pow (with lattice contourplot)
	require(lattice)
	require(RColorBrewer)
	require(gridExtra)
	
	colrmppal <- colorRampPalette(rev(brewer.pal(9, 'PuBu')))

	plot05 <- contourplot(pow05 ~ s*f1, data=power, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.05')
	plot005 <- contourplot(pow005 ~ s*f1, data=power, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main=paste(suffix, '\nPower analysis at p<0.005'))
	plot0005 <- contourplot(pow0005 ~ s*f1, data=power, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=TRUE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.0005')

	
	quartz(width=9, height=3)
	# png(width=9, height=3, filename=paste('analysis/figures/wfs_nullpower_f1_vs_s_with_pow', suffix, '.png', sep=''), units='in', res=300)
	print(plot05, position=c(0,0,0.315,1), more=TRUE)
	print(plot005, position=c(0.315,0,0.63,1), more=TRUE)
	print(plot0005, position=c(0.63,0,1,1), more=FALSE)

	dev.off()
	
	
	
##################################################
## Read in and summarize Lof and Can together
#################################################
require(data.table)

# Read in both
	# Lof-like simulations
locnms <- fread('analysis/Frequency_table_PowerSims_Lof_Ne46000_cnt46_44.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from simulations
	locnms1 <- locnms
load('analysis/wfs_nullmodel_pvals_Power_Lof_Ne46000_cnt46_44.rdata') # dat
	dat1 <- dat

	# Can-like simulations
locnms <- fread('analysis/Frequency_table_PowerSims_Can_Ne5900_cnt32_40.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from simulations
	locnms2 <- locnms
load('analysis/wfs_nullmodel_pvals_Power_Can_Ne5900_cnt32_40.rdata') # dat
	dat2 <- dat

suffix <- '_Lof_Can_comb'

# Then continue here
	# make sure the rows match
	nrow(dat1)
	nrow(locnms1)

	nrow(dat2)
	nrow(locnms2)

# merge
locnms1[,locusnum:=1:nrow(locnms1)] # add a locusnumber for plotting and merging
locnms2[,locusnum:=1:nrow(locnms2)] # add a locusnumber for plotting and merging
dat1 <- as.data.table(dat1)
dat2 <- as.data.table(dat2)
dat1 <- merge(dat1, locnms1[,.(locusnum, CHROM, POS, Freq_1, Freq_2, ABS_DIFF)], by='locusnum')
dat2 <- merge(dat2, locnms2[,.(locusnum, CHROM, POS, Freq_1, Freq_2, ABS_DIFF)], by='locusnum')

# re-order columns
setcolorder(dat1, c('CHROM', 'POS', 'locusnum', 'n', 'cnt1', 'cnt2', 'Freq_1', 'Freq_2', 'ABS_DIFF', 'p', 'pmax', 'p.adj'))
setcolorder(dat2, c('CHROM', 'POS', 'locusnum', 'n', 'cnt1', 'cnt2', 'Freq_1', 'Freq_2', 'ABS_DIFF', 'p', 'pmax', 'p.adj'))

# extract s and f1 values from CHROM field
dat1[,s := as.numeric(gsub('s=', '', unlist(strsplit(CHROM,split=','))[seq(1,2*nrow(dat1),by=2)]))]
dat1[,f1 := as.numeric(gsub('f1=', '', unlist(strsplit(CHROM,split=','))[seq(2,2*nrow(dat1),by=2)]))]

dat2[,s := as.numeric(gsub('s=', '', unlist(strsplit(CHROM,split=','))[seq(1,2*nrow(dat2),by=2)]))]
dat2[,f1 := as.numeric(gsub('f1=', '', unlist(strsplit(CHROM,split=','))[seq(2,2*nrow(dat2),by=2)]))]

# summarize fraction p<0.05, 0.005, 0.0005 for each combination of s and f1
power1 <- dat1[,.(pow05=sum(p<0.05)/.N, pow005=sum(p<0.005)/.N, pow0005=sum(p<0.0005)/.N),by=.(s,f1)]
power2 <- dat2[,.(pow05=sum(p<0.05)/.N, pow005=sum(p<0.005)/.N, pow0005=sum(p<0.0005)/.N),by=.(s,f1)]

# Examine false-positive rate
power1[s==0,]
power2[s==0,]

power1[s==0,summary(pow05)]
power2[s==0,summary(pow05)]

power1[s==0, plot(f1, pow05)]
power1[s==0, lines(f1, predict(loess(pow05~f1)))]

# combined p-values
	# make sure rows match
	all(dat1$s == dat2$s)
	all(dat1$f1 == dat2$f1)
	all(dat1$locusnum == dat2$locusnum)

setkey(dat1, locusnum)
setkey(dat2, locusnum)
dat <- dat1[dat2, ] # column names of dat2 renamed to i._____
dat[, p.comb := pchisq(-2 * sum(log(c(pmax, i.pmax))),df=2*2,lower=FALSE), by=1:nrow(dat)]

power <- dat[,.(pow05=sum(p.comb<0.05)/.N, pow005=sum(p.comb<0.005)/.N, pow0005=sum(p.comb<0.0005)/.N),by=.(s,f1)]


### Plots ###

# Lof and Can separately: f1 vs. s with pow (with lattice contourplot)
	require(lattice)
	require(RColorBrewer)
	require(gridExtra)
	
	colrmppal <- colorRampPalette(rev(brewer.pal(9, 'PuBu')))

	plot05_1 <- contourplot(pow05 ~ s*f1, data=power1, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.05')
	plot005_1 <- contourplot(pow005 ~ s*f1, data=power1, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main=paste('Lofoten', '\nPower analysis at p<0.005'))
	plot0005_1 <- contourplot(pow0005 ~ s*f1, data=power1, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=TRUE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.0005')
	plot05_2 <- contourplot(pow05 ~ s*f1, data=power2, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.05')
	plot005_2 <- contourplot(pow005 ~ s*f1, data=power2, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main=paste('Canada', '\nPower analysis at p<0.005'))
	plot0005_2 <- contourplot(pow0005 ~ s*f1, data=power2, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=TRUE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.0005')

	
	quartz(width=9, height=6)
	# png(width=9, height=6, filename='analysis/figures/wfs_nullpower_f1_vs_s_with_pow_Lof&Can.png', units='in', res=300)
	print(plot05_1, position=c(0,0.5,0.315,1), more=TRUE)
	print(plot005_1, position=c(0.315,0.5,0.63,1), more=TRUE)
	print(plot0005_1, position=c(0.63,0.5,1,1), more=TRUE)
	print(plot05_2, position=c(0,0,0.315,0.5), more=TRUE)
	print(plot005_2, position=c(0.315,0,0.63,0.5), more=TRUE)
	print(plot0005_2, position=c(0.63,0,1,0.5), more=FALSE)

	dev.off()
	
# Lof & Can combined p-value f1 vs. s with pow (with lattice contourplot)
	require(lattice)
	require(RColorBrewer)
	require(gridExtra)
	
	colrmppal <- colorRampPalette(rev(brewer.pal(9, 'PuBu')))

	plot05 <- contourplot(pow05 ~ s*f1, data=power, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.05')
	plot005 <- contourplot(pow005 ~ s*f1, data=power, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main=paste(suffix, '\nPower analysis at p<0.005'))
	plot0005 <- contourplot(pow0005 ~ s*f1, data=power, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=TRUE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.0005')

	
	quartz(width=9, height=3)
	# png(width=9, height=3, filename=paste('analysis/figures/wfs_nullpower_f1_vs_s_with_pow', suffix, '.png', sep=''), units='in', res=300)
	print(plot05, position=c(0,0,0.315,1), more=TRUE)
	print(plot005, position=c(0.315,0,0.63,1), more=TRUE)
	print(plot0005, position=c(0.63,0,1,1), more=FALSE)

	dev.off()
	
	
##################################################
## Read in and summarize Lof and Can together
## Simulate Lof 07-11 and 07-14 by reshuffling p-values
#################################################
require(data.table)

# Read in both
	# Lof-like simulations
locnms <- fread('analysis/Frequency_table_PowerSims_Lof_Ne46000_cnt46_44.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from simulations
	locnms1 <- locnms
load('analysis/wfs_nullmodel_pvals_Power_Lof_Ne46000_cnt46_44.rdata') # dat
	dat1 <- dat

	# Can-like simulations
locnms <- fread('analysis/Frequency_table_PowerSims_Can_Ne5900_cnt32_40.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from simulations
	locnms2 <- locnms
load('analysis/wfs_nullmodel_pvals_Power_Can_Ne5900_cnt32_40.rdata') # dat
	dat2 <- dat

suffix <- '_Lof07-11_Lof-07-14_Can_comb'

# Then continue here
	# make sure the rows match
	nrow(dat1)
	nrow(locnms1)

	nrow(dat2)
	nrow(locnms2)

# merge
locnms1[,locusnum:=1:nrow(locnms1)] # add a locusnumber for plotting and merging
locnms2[,locusnum:=1:nrow(locnms2)] # add a locusnumber for plotting and merging
dat1 <- as.data.table(dat1)
dat2 <- as.data.table(dat2)
dat1 <- merge(dat1, locnms1[,.(locusnum, CHROM, POS, Freq_1, Freq_2, ABS_DIFF)], by='locusnum')
dat2 <- merge(dat2, locnms2[,.(locusnum, CHROM, POS, Freq_1, Freq_2, ABS_DIFF)], by='locusnum')

# re-order columns
setcolorder(dat1, c('CHROM', 'POS', 'locusnum', 'n', 'cnt1', 'cnt2', 'Freq_1', 'Freq_2', 'ABS_DIFF', 'p', 'pmax', 'p.adj'))
setcolorder(dat2, c('CHROM', 'POS', 'locusnum', 'n', 'cnt1', 'cnt2', 'Freq_1', 'Freq_2', 'ABS_DIFF', 'p', 'pmax', 'p.adj'))

# extract s and f1 values from CHROM field
dat1[,s := as.numeric(gsub('s=', '', unlist(strsplit(CHROM,split=','))[seq(1,2*nrow(dat1),by=2)]))]
dat1[,f1 := as.numeric(gsub('f1=', '', unlist(strsplit(CHROM,split=','))[seq(2,2*nrow(dat1),by=2)]))]

dat2[,s := as.numeric(gsub('s=', '', unlist(strsplit(CHROM,split=','))[seq(1,2*nrow(dat2),by=2)]))]
dat2[,f1 := as.numeric(gsub('f1=', '', unlist(strsplit(CHROM,split=','))[seq(2,2*nrow(dat2),by=2)]))]

# simulate Lof 07-11 by reshuffling p-values
dat3 <- dat1[,.(locusnum, s, f1, pmax)]
dat3[, pmax2 := sample(pmax, .N, replace=TRUE), by=.(s,f1)]
setkey(dat3, locusnum)

# combined p-values
	# make sure rows match
	all(dat1$s == dat2$s)
	all(dat1$f1 == dat2$f1)
	all(dat1$s == dat3$s)
	all(dat1$f1 == dat3$f1)

setkey(dat1, locusnum)
setkey(dat2, locusnum)
setkey(dat3, locusnum)
dat4 <- dat1[dat2, ] # column names of dat2 renamed to i._____
dat4 <- dat4[dat3[.(locusnum, s, f1, pmax2)],]
dat4[, p.comb := pchisq(-2 * sum(log(c(pmax, i.pmax, pmax2))),df=2*3,lower=FALSE), by=1:nrow(dat4)]

power4 <- dat4[,.(pow05=sum(p.comb<0.05)/.N, pow005=sum(p.comb<0.005)/.N, pow0005=sum(p.comb<0.0005)/.N),by=.(s,f1)]


### Plots ###
	
# Lof07-11, Lof07-14 & Can combined p-value f1 vs. s with pow (with lattice contourplot)
	require(lattice)
	require(RColorBrewer)
	require(gridExtra)
	
	colrmppal <- colorRampPalette(rev(brewer.pal(9, 'PuBu')))

	plot05 <- contourplot(pow05 ~ s*f1, data=power4, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.05')
	plot005 <- contourplot(pow005 ~ s*f1, data=power4, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=FALSE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main=paste(suffix, '\nPower analysis at p<0.005'))
	plot0005 <- contourplot(pow0005 ~ s*f1, data=power4, at=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01), colorkey=TRUE, col.regions=colrmppal, labels=list(cex=0.7), region=TRUE, main='\nPower analysis at p<0.0005')

	
	quartz(width=9, height=3)
	# png(width=9, height=3, filename=paste('analysis/figures/wfs_nullpower_f1_vs_s_with_pow', suffix, '.png', sep=''), units='in', res=300)
	print(plot05, position=c(0,0,0.315,1), more=TRUE)
	print(plot005, position=c(0.315,0,0.63,1), more=TRUE)
	print(plot0005, position=c(0.63,0,1,1), more=FALSE)

	dev.off()