# Compare freq change, LD change, Tajima's D change, pi change

require(data.table)
require(vioplot)

################################
# read in data
################################

# wd <- 1e4; width='10k'
wd <- 3e4; wdnm='3e4'; width='30k'


freq <- readRDS(paste('analysis/Frequency_table_ABS_DIFF_runmean', wdnm, '.rds', sep=''))
pi <- readRDS(paste('analysis/pi_change_region_', wdnm, '.rds', sep=''))
tajd <- readRDS(paste('analysis/tajimasD_change_region_', wdnm, '.rds', sep=''))
ld <- readRDS(paste('analysis/ld_change_region_', wdnm, '.rds', sep=''))

# adjust column names
setnames(freq, c('perc0711', 'perc0714', 'percCAN'), c('freq_region_perc0711', 'freq_region_perc0714', 'freq_region_percCAN'))
setnames(freq, c('ABS_DIFF_0711', 'ABS_DIFF_0714', 'ABS_DIFF_Can'), c('freq_diff_0711', 'freq_diff_0714', 'freq_diff_CAN'))
setnames(freq, c('cluster0711', 'cluster0714', 'clusterCAN'), c('freq_region_cluster0711', 'freq_region_cluster0714', 'freq_region_clusterCAN'))
setnames(ld, 'CHR', 'CHROM')

# adjust pi bins (they are +1 compared to others)
pi[,BIN_START:= BIN_START-1]

# merge
bins <- merge(
	freq[, .(CHROM, BIN_START, freq_diff_0711, freq_diff_0714, freq_diff_CAN, freq_region_perc0711, freq_region_perc0714, freq_region_percCAN, freq_region_cluster0711, freq_region_cluster0714, freq_region_clusterCAN)], 
	pi[, .(CHROM, BIN_START, pi_diff_0711, pi_diff_0714, pi_diff_CAN, pi_region_perc0711, pi_region_perc0714, pi_region_percCAN, pi_region_cluster0711, pi_region_cluster0714, pi_region_clusterCAN)], 
	by=c('CHROM', 'BIN_START'), all=TRUE)
	
bins <- merge(bins, 
	tajd[, .(CHROM, BIN_START, D_diff_0711, D_diff_0714, D_diff_CAN, D_region_perc0711, D_region_perc0714, D_region_percCAN, D_region_cluster0711, D_region_cluster0714, D_region_clusterCAN)], 
	by=c('CHROM', 'BIN_START'), all=TRUE)
	
bins <- merge(bins, 
	ld[, .(CHROM, BIN_START, ld_diff_0711, ld_diff_0714, ld_diff_CAN, ld_region_perc0711, ld_region_perc0714, ld_region_percCAN, ld_region_cluster0711, ld_region_cluster0714, ld_region_clusterCAN)], 
	by=c('CHROM', 'BIN_START'), all=TRUE)
	
######################
## Biplots of change
## Within same kind of metric, across populations
######################

cols <- c('grey', 'purple', 'red')

quartz(height=6, width=6)
outfile=paste('figures/region_change_across_pops_', wdnm, '.png', sep='')
outfile
# png(height=6, width=6, units='in', res=300, file=outfile)
par(mfrow=c(4, 3), las=1, mai=c(0.5, 0.6, 0.1, 0.1), tcl=-0.2, mgp=c(2.8,0.4,0), cex.axis=0.7)
# frequency change
bins[freq_region_perc0711<=0.99 & freq_region_percCAN<=0.99, plot(freq_diff_CAN, freq_diff_0711, cex=0.2, xlim=c(0,1), ylim=c(0,1), col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_perc0711>0.99 & freq_region_percCAN<=0.99, points(freq_diff_CAN, freq_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & freq_region_percCAN>0.99,points(freq_diff_CAN, freq_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & freq_region_percCAN>0.99,points(freq_diff_CAN, freq_diff_0711, cex=0.2, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one pop', '1% outlier in both pops'), cex=0.5)

bins[freq_region_perc0714<=0.99 & freq_region_percCAN<=0.99, plot(freq_diff_CAN, freq_diff_0714, cex=0.2, xlim=c(0,1), ylim=c(0,1), col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_perc0714>0.99 & freq_region_percCAN<=0.99, points(freq_diff_CAN, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714<=0.99 & freq_region_percCAN>0.99,points(freq_diff_CAN, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714>0.99 & freq_region_percCAN>0.99,points(freq_diff_CAN, freq_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[freq_region_perc0711<=0.99 & freq_region_perc0714<=0.99, plot(freq_diff_0711, freq_diff_0714, cex=0.2, xlim=c(0,1), ylim=c(0,1), col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region_perc0711>0.99 & freq_region_perc0714<=0.99, points(freq_diff_0711, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & freq_region_perc0714>0.99,points(freq_diff_0711, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99,points(freq_diff_0711, freq_diff_0714, cex=0.2, col=cols[3], pch=16)]

# pi change
xlims <- ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, plot(pi_diff_CAN, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='pi_diff_CAN', line=1.5)
bins[(pi_region_perc0711>0.995 | pi_region_perc0711<0.005) & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, points(pi_diff_CAN, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(pi_diff_CAN, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[(pi_region_perc0711>0.995 | pi_region_perc0711<0.005) & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(pi_diff_CAN, pi_diff_0711, cex=0.2, col=cols[3], pch=16)]

bins[pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, plot(pi_diff_CAN, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='pi_diff_CAN', line=1.5)
bins[(pi_region_perc0714>0.995 | pi_region_perc0714<0.005) & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, points(pi_diff_CAN, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005),points(pi_diff_CAN, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(pi_region_perc0714>0.995 | pi_region_perc0714<0.005) & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005),points(pi_diff_CAN, pi_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, plot(pi_diff_0711, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='pi_diff_0711', line=1.5)
bins[(pi_region_perc0711>0.995 | pi_region_perc0711<0.005) & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, points(pi_diff_0711, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(pi_diff_0711, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(pi_region_perc0711>0.995 | pi_region_perc0711<0.005) & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(pi_diff_0711, pi_diff_0714, cex=0.2, col=cols[3], pch=16)]

# tajd change
xlims <- ylims <- bins[,range(D_diff_CAN, D_diff_0711, D_diff_0714, na.rm=TRUE)]
bins[D_region_perc0711<=0.995 & D_region_perc0711>=0.005 & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, plot(D_diff_CAN, D_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='D_diff_CAN', line=1.5)
bins[(D_region_perc0711>0.995 | D_region_perc0711<0.005) & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, points(D_diff_CAN, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[D_region_perc0711<=0.995 & D_region_perc0711>=0.005 & (D_region_percCAN>0.995 | D_region_percCAN<0.005), points(D_diff_CAN, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[(D_region_perc0711>0.995 | D_region_perc0711<0.005) & (D_region_percCAN>0.995 | D_region_percCAN<0.005), points(D_diff_CAN, D_diff_0711, cex=0.2, col=cols[3], pch=16)]

bins[D_region_perc0714<=0.995 & D_region_perc0714>=0.005 & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, plot(D_diff_CAN, D_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='D_diff_CAN', line=1.5)
bins[(D_region_perc0714>0.995 | D_region_perc0714<0.005) & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, points(D_diff_CAN, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[D_region_perc0714<=0.995 & D_region_perc0714>=0.005 & (D_region_percCAN>0.995 | D_region_percCAN<0.005),points(D_diff_CAN, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(D_region_perc0714>0.995 | D_region_perc0714<0.005) & (D_region_percCAN>0.995 | D_region_percCAN<0.005),points(D_diff_CAN, D_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[D_region_perc0711<=0.995 & D_region_perc0711>=0.005 & D_region_perc0714<=0.995 & D_region_perc0714>=0.005, plot(D_diff_0711, D_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='D_diff_0711', line=1.5)
bins[(D_region_perc0711>0.995 | D_region_perc0711<0.005) & D_region_perc0714<=0.995 & D_region_perc0714>=0.005, points(D_diff_0711, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[D_region_perc0711<=0.995 & D_region_perc0711>=0.005 & (D_region_perc0714>0.995 | D_region_perc0714<0.005), points(D_diff_0711, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(D_region_perc0711>0.995 | D_region_perc0711<0.005) & (D_region_perc0714>0.995 | D_region_perc0714<0.005), points(D_diff_0711, D_diff_0714, cex=0.2, col=cols[3], pch=16)]

# LD change
xlims <- ylims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
bins[ld_region_perc0711<=0.99 & ld_region_percCAN<=0.99, plot(ld_diff_CAN, ld_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_CAN', line=1.5)
bins[ld_region_perc0711>0.99 & ld_region_percCAN<=0.99, points(ld_diff_CAN, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711<=0.99 & ld_region_percCAN>0.99,points(ld_diff_CAN, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711>0.99 & ld_region_percCAN>0.99,points(ld_diff_CAN, ld_diff_0711, cex=0.2, col=cols[3], pch=16)]

bins[ld_region_perc0714<=0.99 & ld_region_percCAN<=0.99, plot(ld_diff_CAN, ld_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_CAN', line=1.5)
bins[ld_region_perc0714>0.99 & ld_region_percCAN<=0.99, points(ld_diff_CAN, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0714<=0.99 & ld_region_percCAN>0.99,points(ld_diff_CAN, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0714>0.99 & ld_region_percCAN>0.99,points(ld_diff_CAN, ld_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[ld_region_perc0711<=0.99 & ld_region_perc0714<=0.99, plot(ld_diff_0711, ld_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0711', line=1.5)
bins[ld_region_perc0711>0.99 & ld_region_perc0714<=0.99, points(ld_diff_0711, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711<=0.99 & ld_region_perc0714>0.99,points(ld_diff_0711, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711>0.99 & ld_region_perc0714>0.99,points(ld_diff_0711, ld_diff_0714, cex=0.2, col=cols[3], pch=16)]

dev.off()


######################
## Biplots of change
## Within same population, across metrics
######################

cols <- c('grey', 'purple', 'red')

quartz(height=6, width=6)
outfile=paste('figures/region_change_across_metrics_', wdnm, '.png', sep='')
outfile
# png(height=6, width=6, units='in', res=300, file=outfile)
par(mfrow=c(4, 3), las=1, mai=c(0.5, 0.6, 0.1, 0.1), tcl=-0.2, mgp=c(2.8,0.4,0), cex.axis=0.7)
# frequency change vs. pi change
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[freq_region_perc0711<=0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, plot(freq_diff_0711, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region_perc0711>0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, points(freq_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(freq_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(freq_diff_0711, pi_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region_perc0714<=0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, plot(freq_diff_0714, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region_perc0714>0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, points(freq_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714<=0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(freq_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714>0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(freq_diff_0714, pi_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region_percCAN<=0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, plot(freq_diff_CAN, pi_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_percCAN>0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, points(freq_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN<=0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(freq_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN>0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(freq_diff_CAN, pi_diff_CAN, cex=0.4, col=cols[3], pch=16)]

# freq vs. D
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(D_diff_CAN, D_diff_0711, D_diff_0714, na.rm=TRUE)]
bins[freq_region_perc0711<=0.99 & D_region_perc0711<=0.995 & D_region_perc0711>=0.005, plot(freq_diff_0711, D_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region_perc0711>0.99 & D_region_perc0711<=0.995 & D_region_perc0711>=0.005, points(freq_diff_0711, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & (D_region_perc0711>0.995 | D_region_perc0711<0.005), points(freq_diff_0711, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & (D_region_perc0711>0.995 | D_region_perc0711<0.005), points(freq_diff_0711, D_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region_perc0714<=0.99 & D_region_perc0714<=0.995 & D_region_perc0714>=0.005, plot(freq_diff_0714, D_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region_perc0714>0.99 & D_region_perc0714<=0.995 & D_region_perc0714>=0.005, points(freq_diff_0714, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714<=0.99 & (D_region_perc0714>0.995 | D_region_perc0714<0.005), points(freq_diff_0714, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714>0.99 & (D_region_perc0714>0.995 | D_region_perc0714<0.005), points(freq_diff_0714, D_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region_percCAN<=0.99 & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, plot(freq_diff_CAN, D_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_percCAN>0.99 & D_region_percCAN<=0.995 & D_region_percCAN>=0.005, points(freq_diff_CAN, D_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN<=0.99 & (D_region_percCAN>0.995 | D_region_percCAN<0.005), points(freq_diff_CAN, D_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN>0.99 & (D_region_percCAN>0.995 | D_region_percCAN<0.005), points(freq_diff_CAN, D_diff_CAN, cex=0.4, col=cols[3], pch=16)]

# freq vs. ld
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
bins[freq_region_perc0711<=0.99 & ld_region_perc0711<=0.99, plot(freq_diff_0711, ld_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region_perc0711>0.99 & ld_region_perc0711<=0.99, points(freq_diff_0711, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711<=0.99 & ld_region_perc0711>0.99, points(freq_diff_0711, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0711>0.99 & ld_region_perc0711>0.99, points(freq_diff_0711, ld_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region_perc0714<=0.99 & ld_region_perc0714<=0.99, plot(freq_diff_0714, ld_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region_perc0714>0.99 & ld_region_perc0714<=0.99, points(freq_diff_0714, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714<=0.99 & ld_region_perc0714>0.99, points(freq_diff_0714, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_perc0714>0.99 & ld_region_perc0714>0.99, points(freq_diff_0714, ld_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region_percCAN<=0.99 & ld_region_percCAN<=0.99, plot(freq_diff_CAN, ld_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region_percCAN>0.99 & ld_region_percCAN<=0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN<=0.99 & ld_region_percCAN>0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region_percCAN>0.99 & ld_region_percCAN>0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.4, col=cols[3], pch=16)]


# ld vs. pi
xlims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[ld_region_perc0711<=0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, plot(ld_diff_0711, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0711', line=1.5)
bins[ld_region_perc0711>0.99 & pi_region_perc0711<=0.995 & pi_region_perc0711>=0.005, points(ld_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711<=0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(ld_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0711>0.99 & (pi_region_perc0711>0.995 | pi_region_perc0711<0.005), points(ld_diff_0711, pi_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[ld_region_perc0714<=0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, plot(ld_diff_0714, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0714', line=1.5)
bins[ld_region_perc0714>0.99 & pi_region_perc0714<=0.995 & pi_region_perc0714>=0.005, points(ld_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0714<=0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(ld_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_perc0714>0.99 & (pi_region_perc0714>0.995 | pi_region_perc0714<0.005), points(ld_diff_0714, pi_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[ld_region_percCAN<=0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, plot(ld_diff_CAN, pi_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_CAN', line=1.5)
bins[ld_region_percCAN>0.99 & pi_region_percCAN<=0.995 & pi_region_percCAN>=0.005, points(ld_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_percCAN<=0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(ld_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[ld_region_percCAN>0.99 & (pi_region_percCAN>0.995 | pi_region_percCAN<0.005), points(ld_diff_CAN, pi_diff_CAN, cex=0.4, col=cols[3], pch=16)]

dev.off()


#####################
# Find outlier regions
#####################

#### regions with large frequency change and large increase in LD
	bins[,freqLD0711 := freq_region_perc0711 + ld_region_perc0711]
	bins[,freqLD0714 := freq_region_perc0714 + ld_region_perc0714]
	bins[,freqLDCAN := freq_region_percCAN + ld_region_percCAN]
	
	# 07-11
	setkey(bins, freqLD0711); bins[freq_region_perc0711>0.99 & ld_region_perc0711>0.99,.(CHROM, BIN_START, freq_diff_0711, ld_diff_0711, freqLD0711, freq_region_perc0711, ld_region_perc0711)]
	# 07-14
	setkey(bins, freqLD0714); bins[freq_region_perc0714>0.99 & ld_region_perc0714>0.99,.(CHROM, BIN_START, freq_diff_0714, ld_diff_0714, freqLD0714, freq_region_perc0714, ld_region_perc0714)]
	# 07-11 and 07-14
	setkey(bins, freqLD0711, freqLD0714); bins[freq_region_perc0711>0.95 & ld_region_perc0711>0.95 & freq_region_perc0714>0.95 & ld_region_perc0714>0.95,.(CHROM, BIN_START, freq_diff_0711, freq_diff_0714, ld_diff_0711, ld_diff_0714, freqLD0711, freqLD0714)]
	# CAN
	setkey(bins, freqLDCAN); bins[freq_region_percCAN>0.99 & ld_region_percCAN>0.99,.(CHROM, BIN_START, freq_diff_CAN, ld_diff_CAN, freqLDCAN, freq_region_percCAN, ld_region_percCAN)]
	# 07-11 and 07-14 and CAN
	setkey(bins, freqLD0711, freqLD0714, freqLDCAN); bins[freq_region_perc0711>0.95 & ld_region_perc0711>0.95 & freq_region_perc0714>0.95 & ld_region_perc0714>0.95 & freq_region_percCAN>0.95 & ld_region_percCAN>0.95,.(CHROM, BIN_START, freq_diff_0711, freq_diff_0714, freq_diff_CAN, ld_diff_0711, ld_diff_0714, ld_diff_CAN)]

	par(mfrow=c(1,3)); bins[,plot(freqLD0711, freqLD0714, cex=0.2, col='#00000033')]; bins[,plot(freqLD0711, freqLDCAN, cex=0.2, col='#00000033')]; bins[,plot(freqLD0714, freqLDCAN, cex=0.2, col='#00000033')]
	bins[,cor.test(freqLD0711, freqLD0714)]
	bins[,cor.test(freqLD0711, freqLDCAN)]
	bins[,cor.test(freqLDCAN, freqLD0714)]
	
	
	# write out a selection of highly ranked loci
	setkey(bins, CHROM, BIN_START)

	write.csv(bins[(freq_region_perc0711>0.99 & ld_region_perc0711>0.99) | (freq_region_perc0714>0.99 & ld_region_perc0714>0.99) | (freq_region_percCAN>0.99 & ld_region_percCAN>0.99), .(CHROM, BIN_START, BIN_END=BIN_START+wd, outlier0711=(freq_region_perc0711>0.99 & ld_region_perc0711>0.99), outlier0714=(freq_region_perc0714>0.99 & ld_region_perc0714>0.99), outlierCAN=(freq_region_percCAN>0.99 & ld_region_percCAN>0.99), freq_diff_0711, freq_diff_0714, freq_diff_CAN, ld_diff_0711, ld_diff_0714, ld_diff_CAN)], file='analysis/outlier_10kregions_freqLD_07-11-14_Can.csv', row.names=FALSE)
	
	
#### rank the outlier regions based in clustering the frequency changes together (see Oziolor et al. 2019 Science)
	outlreg0711 = bins[,.(CHROM=unique(CHROM), BIN_START=min(BIN_START), BIN_END=max(BIN_START+wd), WIDTH=max(BIN_START+wd)-min(BIN_START), score=sum(freq_diff_0711)), by=freq_region_cluster0711]
	outlreg0714 = bins[,.(CHROM=unique(CHROM), BIN_START=min(BIN_START), BIN_END=max(BIN_START+wd), WIDTH=max(BIN_START+wd)-min(BIN_START), score=sum(freq_diff_0714)), by=freq_region_cluster0714]
	outlregCAN = bins[,.(CHROM=unique(CHROM), BIN_START=min(BIN_START), BIN_END=max(BIN_START+wd), WIDTH=max(BIN_START+wd)-min(BIN_START), score=sum(freq_diff_CAN)), by=freq_region_clusterCAN]

	setkey(outlreg0711, score)
	setkey(outlreg0714, score)
	setkey(outlregCAN, score)

	outlreg0711[,sort(unique(WIDTH))] # cluster width always 10kb
	outlreg0714[,sort(unique(WIDTH))]
	outlregCAN[,sort(unique(WIDTH))]

	tail(outlreg0711, 10)
	tail(outlreg0714, 10)
	tail(outlregCAN, 10)

#### rank outlier regions based on shared change in allele frequency across populations
	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(CHROM, BIN_START, BIN_END=BIN_START+wd, freq_diff_0711, freq_diff_0714, freq_diff_CAN)]

	# write out a selection of highly ranked loci
	setkey(bins, CHROM, BIN_START)

	outfile <- paste('analysis/outlier_', width, 'regions_freqshared_07-11-14_Can.csv', sep='')
	outfile
	write.csv(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(CHROM, BIN_START, BIN_END=BIN_START+wd)], file=outfile, row.names=FALSE)


		
	# LD change in freq outlier regions
	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(CHROM, BIN_START, BIN_END=BIN_START+wd, freq_diff_0711, freq_diff_0714, freq_diff_CAN, ld_diff_0711, ld_diff_0714, ld_diff_CAN)]

	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(mean0711=mean(ld_diff_0711, na.rm=TRUE), se0711=sd(ld_diff_0711, na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0711))), mean0714=mean(ld_diff_0714, na.rm=TRUE), se0714=sd(ld_diff_0714, na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0714))), meanCAN=mean(ld_diff_CAN, na.rm=TRUE), seCAN=sd(ld_diff_CAN, na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_CAN))))] # in outlier regions
	bins[, .(mean0711=mean(ld_diff_0711, na.rm=TRUE), se0711=sd(ld_diff_0711, na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0711))), mean0714=mean(ld_diff_0714, na.rm=TRUE), se0714=sd(ld_diff_0714, na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0714))), meanCAN=mean(ld_diff_CAN, na.rm=TRUE), seCAN=sd(ld_diff_CAN, na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_CAN))))] # genome-wide for comparison
	
		# simple plot
		quartz(height=5, width=9)
		vioplot(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_0711], 
			bins[,ld_diff_0711],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_0714], 
			bins[,ld_diff_0714],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_CAN], 
			bins[,ld_diff_CAN],
			col=c('purple', 'grey', 'purple', 'grey', 'purple', 'grey'),
			names=c('0711 outliers', '0711 genome', '0714 outliers', '0714 genome', 'CAN outliers', 'CAN genome'),
			ylab="Change in LD",
			main='Regions that are 99% frequency change outliers in all three comparisons')
		abline(h=0, lty=2, col='grey')
		

	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_0711], bins[,ld_diff_0711])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_0714], bins[,ld_diff_0714])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_CAN], bins[,ld_diff_CAN])

		# one-sample version
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_0711], bins[,ld_diff_0711])
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_0714], bins[,ld_diff_0714])
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, ld_diff_CAN], bins[,ld_diff_CAN])


	# abs(LD change)
	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(mean0711=mean(abs(ld_diff_0711), na.rm=TRUE), se0711=sd(abs(ld_diff_0711), na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0711))), mean0711=mean(abs(ld_diff_0714), na.rm=TRUE), se0714=sd(abs(ld_diff_0714), na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0714))), meanCAN=mean(abs(ld_diff_CAN), na.rm=TRUE), seCAN=sd(abs(ld_diff_CAN), na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_CAN))))] # in outlier regions
	bins[, .(mean0711=mean(abs(ld_diff_0711), na.rm=TRUE), se0711=sd(abs(ld_diff_0711), na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0711))), mean0711=mean(abs(ld_diff_0714), na.rm=TRUE), se0714=sd(abs(ld_diff_0714), na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_0714))), meanCAN=mean(abs(ld_diff_CAN), na.rm=TRUE), seCAN=sd(abs(ld_diff_CAN), na.rm=TRUE)/sqrt(sum(!is.na(ld_diff_CAN))))] # genome-wide for comparison
	
		# simple plot
		quartz(height=5, width=9)
		vioplot(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(ld_diff_0711)], 
			bins[,abs(ld_diff_0711)],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(ld_diff_0714)], 
			bins[,abs(ld_diff_0714)],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(ld_diff_CAN)], 
			bins[,abs(ld_diff_CAN)],
			col=c('purple', 'grey', 'purple', 'grey', 'purple', 'grey'),
			names=c('0711 outliers', '0711 genome', '0714 outliers', '0714 genome', 'CAN outliers', 'CAN genome'),
			ylab="abs(Change in LD)",
			main='Regions that are 99% frequency change outliers in all three comparisons')
		abline(h=0, lty=2, col='grey')

	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(ld_diff_0711)], bins[,abs(ld_diff_0711)])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(ld_diff_0714)], bins[,abs(ld_diff_0714)])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(ld_diff_CAN)], bins[,abs(ld_diff_CAN)])

		
	# Tajima's D change
	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(CHROM, BIN_START, BIN_END=BIN_START+wd, freq_diff_0711, freq_diff_0714, freq_diff_CAN, D_diff_0711, D_diff_0714, D_diff_CAN)]

	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(mean0711=mean(D_diff_0711, na.rm=TRUE), se0711=sd(D_diff_0711, na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0711))), mean0711=mean(D_diff_0714, na.rm=TRUE), se0714=sd(D_diff_0714, na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0714))), meanCAN=mean(D_diff_CAN, na.rm=TRUE), seCAN=sd(D_diff_CAN, na.rm=TRUE)/sqrt(sum(!is.na(D_diff_CAN))))] # in outlier regions
	bins[, .(mean0711=mean(D_diff_0711, na.rm=TRUE), se0711=sd(D_diff_0711, na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0711))), mean0711=mean(D_diff_0714, na.rm=TRUE), se0714=sd(D_diff_0714, na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0714))), meanCAN=mean(D_diff_CAN, na.rm=TRUE), seCAN=sd(D_diff_CAN, na.rm=TRUE)/sqrt(sum(!is.na(D_diff_CAN))))] # genome-wide for comparison
	
		# simple plot
		quartz(height=5, width=9)
		outfile=paste('figures/region_change_shared_freq_outliers_Dchange_', wdnm, '.pdf', sep='')
		outfile
		# pdf(height=5, width=9, file=outfile)
		vioplot(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_0711], 
			bins[,D_diff_0711],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_0714], 
			bins[,D_diff_0714],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_CAN], 
			bins[,D_diff_CAN],
			col=c('purple', 'grey', 'purple', 'grey', 'purple', 'grey'),
			names=c('0711 outliers', '0711 genome', '0714 outliers', '0714 genome', 'CAN outliers', 'CAN genome'),
			ylab="Change in Tajima's D",
			main='Regions that are 99% frequency change outliers in all three comparisons')
		abline(h=0, lty=2, col='grey')
		
		dev.off()

	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_0711], bins[,D_diff_0711])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_0714], bins[,D_diff_0714])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_CAN], bins[,D_diff_CAN])

		# one-sample version
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_0711], bins[,D_diff_0711])
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_0714], bins[,D_diff_0714])
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, D_diff_CAN], bins[,D_diff_CAN])

	

	# abs(Tajima's D change)
	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(mean0711=mean(abs(D_diff_0711), na.rm=TRUE), se0711=sd(abs(D_diff_0711), na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0711))), mean0711=mean(abs(D_diff_0714), na.rm=TRUE), se0714=sd(abs(D_diff_0714), na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0714))), meanCAN=mean(abs(D_diff_CAN), na.rm=TRUE), seCAN=sd(abs(D_diff_CAN), na.rm=TRUE)/sqrt(sum(!is.na(D_diff_CAN))))] # in outlier regions
	bins[, .(mean0711=mean(abs(D_diff_0711), na.rm=TRUE), se0711=sd(abs(D_diff_0711), na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0711))), mean0711=mean(abs(D_diff_0714), na.rm=TRUE), se0714=sd(abs(D_diff_0714), na.rm=TRUE)/sqrt(sum(!is.na(D_diff_0714))), meanCAN=mean(abs(D_diff_CAN), na.rm=TRUE), seCAN=sd(abs(D_diff_CAN), na.rm=TRUE)/sqrt(sum(!is.na(D_diff_CAN))))] # genome-wide for comparison
	
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(D_diff_0711)], bins[,abs(D_diff_0711)])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(D_diff_0714)], bins[,abs(D_diff_0714)])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, abs(D_diff_CAN)], bins[,abs(D_diff_CAN)])


	# Pi change
	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(CHROM, BIN_START, BIN_END=BIN_START+wd, freq_diff_0711, freq_diff_0714, freq_diff_CAN, pi_diff_0711, pi_diff_0714, pi_diff_CAN)]

	bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, .(mean0711=mean(pi_diff_0711, na.rm=TRUE), se0711=sd(pi_diff_0711, na.rm=TRUE)/sqrt(sum(!is.na(pi_diff_0711))), mean0711=mean(pi_diff_0714, na.rm=TRUE), se0714=sd(pi_diff_0714, na.rm=TRUE)/sqrt(sum(!is.na(pi_diff_0714))), meanCAN=mean(pi_diff_CAN, na.rm=TRUE), seCAN=sd(pi_diff_CAN, na.rm=TRUE)/sqrt(sum(!is.na(pi_diff_CAN))))] # in outlier regions
	bins[, .(mean0711=mean(pi_diff_0711, na.rm=TRUE), se0711=sd(pi_diff_0711, na.rm=TRUE)/sqrt(sum(!is.na(pi_diff_0711))), mean0711=mean(pi_diff_0714, na.rm=TRUE), se0714=sd(pi_diff_0714, na.rm=TRUE)/sqrt(sum(!is.na(pi_diff_0714))), meanCAN=mean(pi_diff_CAN, na.rm=TRUE), seCAN=sd(pi_diff_CAN, na.rm=TRUE)/sqrt(sum(!is.na(pi_diff_CAN))))] # genome-wide for comparison
	
		# simple plot
		quartz(height=5, width=9)
		outfile=paste('figures/region_change_shared_freq_outliers_pichange_', wdnm, '.pdf', sep='')
		outfile
		# pdf(height=5, width=9, file=outfile)
		vioplot(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_0711], 
			bins[,pi_diff_0711],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_0714], 
			bins[,pi_diff_0714],
			bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_CAN], 
			bins[,pi_diff_CAN],
			col=c('purple', 'grey', 'purple', 'grey', 'purple', 'grey'),
			names=c('0711 outliers', '0711 genome', '0714 outliers', '0714 genome', 'CAN outliers', 'CAN genome'),
			ylab="Change in pi",
			main='Regions that are 99% frequency change outliers in all three comparisons')
		abline(h=0, lty=2, col='grey')
		
		dev.off()

	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_0711], bins[,pi_diff_0711])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_0714], bins[,pi_diff_0714])
	t.test(bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_CAN], bins[,pi_diff_CAN])

		# one-sample
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_0711], bins[,pi_diff_0711])
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_0714], bins[,pi_diff_0714])
	t.test(mu=bins[freq_region_perc0711>0.99 & freq_region_perc0714>0.99 & freq_region_percCAN>0.99, pi_diff_CAN], bins[,pi_diff_CAN])
