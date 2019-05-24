# Compare freq change, LD change, Tajima's D change, pi change

# read in data
width='1e4'
freq <- readRDS('analysis/Frequency_table_ABS_DIFF_runmean1e4.rds')
pi <- readRDS('analysis/pi_change_region_1e4.rds')
tajd <- readRDS('analysis/tajimasD_change_region_1e4.rds')
ld <- readRDS('analysis/ld_change_region_1e4.rds')

# adjust column names
setnames(freq, c('perc0711', 'perc0714', 'percCAN'), c('freq_region10kb_perc0711', 'freq_region10kb_perc0714', 'freq_region10kb_percCAN'))
setnames(freq, c('ABS_DIFF_0711', 'ABS_DIFF_0714', 'ABS_DIFF_Can'), c('freq_diff_0711', 'freq_diff_0714', 'freq_diff_CAN'))
setnames(freq, c('cluster0711', 'cluster0714', 'clusterCAN'), c('freq_region10kb_cluster0711', 'freq_region10kb_cluster0714', 'freq_region10kb_clusterCAN'))
setnames(ld, 'CHR', 'CHROM')

# adjust pi bins (they are +1 compared to others)
pi[,BIN_START:= BIN_START-1]

# merge
bins <- merge(
	freq[, .(CHROM, BIN_START, freq_diff_0711, freq_diff_0714, freq_diff_CAN, freq_region10kb_perc0711, freq_region10kb_perc0714, freq_region10kb_percCAN, freq_region10kb_cluster0711, freq_region10kb_cluster0714, freq_region10kb_clusterCAN)], 
	pi[, .(CHROM, BIN_START, pi_diff_0711, pi_diff_0714, pi_diff_CAN, pi_region10kb_perc0711, pi_region10kb_perc0714, pi_region10kb_percCAN, pi_region10kb_cluster0711, pi_region10kb_cluster0714, pi_region10kb_clusterCAN)], 
	by=c('CHROM', 'BIN_START'), all=TRUE)
	
bins <- merge(bins, 
	tajd[, .(CHROM, BIN_START, D_diff_0711, D_diff_0714, D_diff_CAN, D_region10kb_perc0711, D_region10kb_perc0714, D_region10kb_percCAN, D_region10kb_cluster0711, D_region10kb_cluster0714, D_region10kb_clusterCAN)], 
	by=c('CHROM', 'BIN_START'), all=TRUE)
	
bins <- merge(bins, 
	ld[, .(CHROM, BIN_START, ld_diff_0711, ld_diff_0714, ld_diff_CAN, ld_region10kb_perc0711, ld_region10kb_perc0714, ld_region10kb_percCAN, ld_region10kb_cluster0711, ld_region10kb_cluster0714, ld_region10kb_clusterCAN)], 
	by=c('CHROM', 'BIN_START'), all=TRUE)
	
######################
## Biplots of change
## Within same kind of metric, across populations
######################

cols <- c('grey', 'purple', 'red')

quartz(height=6, width=6)
# png(height=6, width=6, units='in', res=300, file='figures/region_change_across_pops.png')
par(mfrow=c(4, 3), las=1, mai=c(0.5, 0.6, 0.1, 0.1), tcl=-0.2, mgp=c(2.8,0.4,0), cex.axis=0.7)
# frequency change
bins[freq_region10kb_perc0711<=0.99 & freq_region10kb_percCAN<=0.99, plot(freq_diff_CAN, freq_diff_0711, cex=0.2, xlim=c(0,1), ylim=c(0,1), col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region10kb_perc0711>0.99 & freq_region10kb_percCAN<=0.99, points(freq_diff_CAN, freq_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711<=0.99 & freq_region10kb_percCAN>0.99,points(freq_diff_CAN, freq_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711>0.99 & freq_region10kb_percCAN>0.99,points(freq_diff_CAN, freq_diff_0711, cex=0.2, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one pop', '1% outlier in both pops'), cex=0.5)

bins[freq_region10kb_perc0714<=0.99 & freq_region10kb_percCAN<=0.99, plot(freq_diff_CAN, freq_diff_0714, cex=0.2, xlim=c(0,1), ylim=c(0,1), col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region10kb_perc0714>0.99 & freq_region10kb_percCAN<=0.99, points(freq_diff_CAN, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714<=0.99 & freq_region10kb_percCAN>0.99,points(freq_diff_CAN, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714>0.99 & freq_region10kb_percCAN>0.99,points(freq_diff_CAN, freq_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[freq_region10kb_perc0711<=0.99 & freq_region10kb_perc0714<=0.99, plot(freq_diff_0711, freq_diff_0714, cex=0.2, xlim=c(0,1), ylim=c(0,1), col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region10kb_perc0711>0.99 & freq_region10kb_perc0714<=0.99, points(freq_diff_0711, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711<=0.99 & freq_region10kb_perc0714>0.99,points(freq_diff_0711, freq_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711>0.99 & freq_region10kb_perc0714>0.99,points(freq_diff_0711, freq_diff_0714, cex=0.2, col=cols[3], pch=16)]

# pi change
xlims <- ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005 & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, plot(pi_diff_CAN, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='pi_diff_CAN', line=1.5)
bins[(pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005) & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, points(pi_diff_CAN, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005 & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005), points(pi_diff_CAN, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[(pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005) & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005), points(pi_diff_CAN, pi_diff_0711, cex=0.2, col=cols[3], pch=16)]

bins[pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005 & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, plot(pi_diff_CAN, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='pi_diff_CAN', line=1.5)
bins[(pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005) & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, points(pi_diff_CAN, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005 & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005),points(pi_diff_CAN, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005) & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005),points(pi_diff_CAN, pi_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005 & pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005, plot(pi_diff_0711, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='pi_diff_0711', line=1.5)
bins[(pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005) & pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005, points(pi_diff_0711, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005 & (pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005), points(pi_diff_0711, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005) & (pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005), points(pi_diff_0711, pi_diff_0714, cex=0.2, col=cols[3], pch=16)]

# tajd change
xlims <- ylims <- bins[,range(D_diff_CAN, D_diff_0711, D_diff_0714, na.rm=TRUE)]
bins[D_region10kb_perc0711<=0.995 & D_region10kb_perc0711>=0.005 & D_region10kb_percCAN<=0.995 & D_region10kb_percCAN>=0.005, plot(D_diff_CAN, D_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='D_diff_CAN', line=1.5)
bins[(D_region10kb_perc0711>0.995 | D_region10kb_perc0711<0.005) & D_region10kb_percCAN<=0.995 & D_region10kb_percCAN>=0.005, points(D_diff_CAN, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[D_region10kb_perc0711<=0.995 & D_region10kb_perc0711>=0.005 & (D_region10kb_percCAN>0.995 | D_region10kb_percCAN<0.005), points(D_diff_CAN, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[(D_region10kb_perc0711>0.995 | D_region10kb_perc0711<0.005) & (D_region10kb_percCAN>0.995 | D_region10kb_percCAN<0.005), points(D_diff_CAN, D_diff_0711, cex=0.2, col=cols[3], pch=16)]

bins[D_region10kb_perc0714<=0.995 & D_region10kb_perc0714>=0.005 & D_region10kb_percCAN<=0.995 & D_region10kb_percCAN>=0.005, plot(D_diff_CAN, D_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='D_diff_CAN', line=1.5)
bins[(D_region10kb_perc0714>0.995 | D_region10kb_perc0714<0.005) & D_region10kb_percCAN<=0.995 & D_region10kb_percCAN>=0.005, points(D_diff_CAN, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[D_region10kb_perc0714<=0.995 & D_region10kb_perc0714>=0.005 & (D_region10kb_percCAN>0.995 | D_region10kb_percCAN<0.005),points(D_diff_CAN, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(D_region10kb_perc0714>0.995 | D_region10kb_perc0714<0.005) & (D_region10kb_percCAN>0.995 | D_region10kb_percCAN<0.005),points(D_diff_CAN, D_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[D_region10kb_perc0711<=0.995 & D_region10kb_perc0711>=0.005 & D_region10kb_perc0714<=0.995 & D_region10kb_perc0714>=0.005, plot(D_diff_0711, D_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='D_diff_0711', line=1.5)
bins[(D_region10kb_perc0711>0.995 | D_region10kb_perc0711<0.005) & D_region10kb_perc0714<=0.995 & D_region10kb_perc0714>=0.005, points(D_diff_0711, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[D_region10kb_perc0711<=0.995 & D_region10kb_perc0711>=0.005 & (D_region10kb_perc0714>0.995 | D_region10kb_perc0714<0.005), points(D_diff_0711, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[(D_region10kb_perc0711>0.995 | D_region10kb_perc0711<0.005) & (D_region10kb_perc0714>0.995 | D_region10kb_perc0714<0.005), points(D_diff_0711, D_diff_0714, cex=0.2, col=cols[3], pch=16)]

# LD change
xlims <- ylims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
bins[ld_region10kb_perc0711<=0.99 & ld_region10kb_percCAN<=0.99, plot(ld_diff_CAN, ld_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_CAN', line=1.5)
bins[ld_region10kb_perc0711>0.99 & ld_region10kb_percCAN<=0.99, points(ld_diff_CAN, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0711<=0.99 & ld_region10kb_percCAN>0.99,points(ld_diff_CAN, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0711>0.99 & ld_region10kb_percCAN>0.99,points(ld_diff_CAN, ld_diff_0711, cex=0.2, col=cols[3], pch=16)]

bins[ld_region10kb_perc0714<=0.99 & ld_region10kb_percCAN<=0.99, plot(ld_diff_CAN, ld_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_CAN', line=1.5)
bins[ld_region10kb_perc0714>0.99 & ld_region10kb_percCAN<=0.99, points(ld_diff_CAN, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0714<=0.99 & ld_region10kb_percCAN>0.99,points(ld_diff_CAN, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0714>0.99 & ld_region10kb_percCAN>0.99,points(ld_diff_CAN, ld_diff_0714, cex=0.2, col=cols[3], pch=16)]

bins[ld_region10kb_perc0711<=0.99 & ld_region10kb_perc0714<=0.99, plot(ld_diff_0711, ld_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0711', line=1.5)
bins[ld_region10kb_perc0711>0.99 & ld_region10kb_perc0714<=0.99, points(ld_diff_0711, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0711<=0.99 & ld_region10kb_perc0714>0.99,points(ld_diff_0711, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0711>0.99 & ld_region10kb_perc0714>0.99,points(ld_diff_0711, ld_diff_0714, cex=0.2, col=cols[3], pch=16)]

dev.off()


######################
## Biplots of change
## Within same population, across metrics
######################

cols <- c('grey', 'purple', 'red')

quartz(height=6, width=6)
# png(height=6, width=6, units='in', res=300, file='figures/region_change_across_metrics.png')
par(mfrow=c(4, 3), las=1, mai=c(0.5, 0.6, 0.1, 0.1), tcl=-0.2, mgp=c(2.8,0.4,0), cex.axis=0.7)
# frequency change vs. pi change
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[freq_region10kb_perc0711<=0.99 & pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005, plot(freq_diff_0711, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region10kb_perc0711>0.99 & pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005, points(freq_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711<=0.99 & (pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005), points(freq_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711>0.99 & (pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005), points(freq_diff_0711, pi_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region10kb_perc0714<=0.99 & pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005, plot(freq_diff_0714, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region10kb_perc0714>0.99 & pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005, points(freq_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714<=0.99 & (pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005), points(freq_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714>0.99 & (pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005), points(freq_diff_0714, pi_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region10kb_percCAN<=0.99 & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, plot(freq_diff_CAN, pi_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region10kb_percCAN>0.99 & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, points(freq_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_percCAN<=0.99 & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005), points(freq_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_percCAN>0.99 & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005), points(freq_diff_CAN, pi_diff_CAN, cex=0.4, col=cols[3], pch=16)]

# freq vs. D
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(D_diff_CAN, D_diff_0711, D_diff_0714, na.rm=TRUE)]
bins[freq_region10kb_perc0711<=0.99 & D_region10kb_perc0711<=0.995 & D_region10kb_perc0711>=0.005, plot(freq_diff_0711, D_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region10kb_perc0711>0.99 & D_region10kb_perc0711<=0.995 & D_region10kb_perc0711>=0.005, points(freq_diff_0711, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711<=0.99 & (D_region10kb_perc0711>0.995 | D_region10kb_perc0711<0.005), points(freq_diff_0711, D_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711>0.99 & (D_region10kb_perc0711>0.995 | D_region10kb_perc0711<0.005), points(freq_diff_0711, D_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region10kb_perc0714<=0.99 & D_region10kb_perc0714<=0.995 & D_region10kb_perc0714>=0.005, plot(freq_diff_0714, D_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region10kb_perc0714>0.99 & D_region10kb_perc0714<=0.995 & D_region10kb_perc0714>=0.005, points(freq_diff_0714, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714<=0.99 & (D_region10kb_perc0714>0.995 | D_region10kb_perc0714<0.005), points(freq_diff_0714, D_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714>0.99 & (D_region10kb_perc0714>0.995 | D_region10kb_perc0714<0.005), points(freq_diff_0714, D_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region10kb_percCAN<=0.99 & D_region10kb_percCAN<=0.995 & D_region10kb_percCAN>=0.005, plot(freq_diff_CAN, D_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region10kb_percCAN>0.99 & D_region10kb_percCAN<=0.995 & D_region10kb_percCAN>=0.005, points(freq_diff_CAN, D_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_percCAN<=0.99 & (D_region10kb_percCAN>0.995 | D_region10kb_percCAN<0.005), points(freq_diff_CAN, D_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_percCAN>0.99 & (D_region10kb_percCAN>0.995 | D_region10kb_percCAN<0.005), points(freq_diff_CAN, D_diff_CAN, cex=0.4, col=cols[3], pch=16)]

# freq vs. ld
xlims <- bins[,range(freq_diff_CAN, freq_diff_0711, freq_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
bins[freq_region10kb_perc0711<=0.99 & ld_region10kb_perc0711<=0.99, plot(freq_diff_0711, ld_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0711', line=1.5)
bins[freq_region10kb_perc0711>0.99 & ld_region10kb_perc0711<=0.99, points(freq_diff_0711, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711<=0.99 & ld_region10kb_perc0711>0.99, points(freq_diff_0711, ld_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0711>0.99 & ld_region10kb_perc0711>0.99, points(freq_diff_0711, ld_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[freq_region10kb_perc0714<=0.99 & ld_region10kb_perc0714<=0.99, plot(freq_diff_0714, ld_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_0714', line=1.5)
bins[freq_region10kb_perc0714>0.99 & ld_region10kb_perc0714<=0.99, points(freq_diff_0714, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714<=0.99 & ld_region10kb_perc0714>0.99, points(freq_diff_0714, ld_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_perc0714>0.99 & ld_region10kb_perc0714>0.99, points(freq_diff_0714, ld_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[freq_region10kb_percCAN<=0.99 & ld_region10kb_percCAN<=0.99, plot(freq_diff_CAN, ld_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='freq_diff_CAN', line=1.5)
bins[freq_region10kb_percCAN>0.99 & ld_region10kb_percCAN<=0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_percCAN<=0.99 & ld_region10kb_percCAN>0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[freq_region10kb_percCAN>0.99 & ld_region10kb_percCAN>0.99, points(freq_diff_CAN, ld_diff_CAN, cex=0.4, col=cols[3], pch=16)]


# ld vs. pi
xlims <- bins[,range(ld_diff_CAN, ld_diff_0711, ld_diff_0714, na.rm=TRUE)]
ylims <- bins[,range(pi_diff_CAN, pi_diff_0711, pi_diff_0714, na.rm=TRUE)]
bins[ld_region10kb_perc0711<=0.99 & pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005, plot(ld_diff_0711, pi_diff_0711, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0711', line=1.5)
bins[ld_region10kb_perc0711>0.99 & pi_region10kb_perc0711<=0.995 & pi_region10kb_perc0711>=0.005, points(ld_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0711<=0.99 & (pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005), points(ld_diff_0711, pi_diff_0711, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0711>0.99 & (pi_region10kb_perc0711>0.995 | pi_region10kb_perc0711<0.005), points(ld_diff_0711, pi_diff_0711, cex=0.4, col=cols[3], pch=16)]

legend('topright', col=cols, pch=16, legend=c('Locus', '1% outlier in one metric', '1% outlier in both metrics'), cex=0.5)

bins[ld_region10kb_perc0714<=0.99 & pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005, plot(ld_diff_0714, pi_diff_0714, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_0714', line=1.5)
bins[ld_region10kb_perc0714>0.99 & pi_region10kb_perc0714<=0.995 & pi_region10kb_perc0714>=0.005, points(ld_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0714<=0.99 & (pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005), points(ld_diff_0714, pi_diff_0714, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_perc0714>0.99 & (pi_region10kb_perc0714>0.995 | pi_region10kb_perc0714<0.005), points(ld_diff_0714, pi_diff_0714, cex=0.4, col=cols[3], pch=16)]

bins[ld_region10kb_percCAN<=0.99 & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, plot(ld_diff_CAN, pi_diff_CAN, cex=0.2, xlim=xlims, ylim=ylims, col=cols[1], pch=16, xlab='')]; title(xlab='ld_diff_CAN', line=1.5)
bins[ld_region10kb_percCAN>0.99 & pi_region10kb_percCAN<=0.995 & pi_region10kb_percCAN>=0.005, points(ld_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_percCAN<=0.99 & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005), points(ld_diff_CAN, pi_diff_CAN, cex=0.2, col=cols[2], pch=16)]
bins[ld_region10kb_percCAN>0.99 & (pi_region10kb_percCAN>0.995 | pi_region10kb_percCAN<0.005), points(ld_diff_CAN, pi_diff_CAN, cex=0.4, col=cols[3], pch=16)]

dev.off()