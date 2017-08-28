# Calculate probability of null model producing results as extreme as our observations
# run after wfs_nullmodel_function.r and wfs_nullmodel_combine.r
# Script analyzes the output

#############################
## Read in and summarize
#############################
require(data.table)

# Choose one set to read in
locnms <- fread('analysis/Frequency_table_PowerSims_Lof_Ne46000_cnt46_44.txt', header=TRUE); setnames(locnms, 3:7, c('N_CHR_1', 'Freq_1', 'N_CHR_2', 'Freq_2', 'ABS_DIFF')) # the name and observed frequencies of all the loci, from output by Bastiaan Star
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