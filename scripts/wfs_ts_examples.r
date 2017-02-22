source('scripts/wfs_ts.r')

# plot drift vs. selection
e1 <- wfs_ts(f1=0.5)
e2 <- wfs_ts(f1=0.5)
e3 <- wfs_ts(f1=0.5)

s1 <- wfs_ts(f1=0.5, s=0.5)
s2 <- wfs_ts(f1=0.5, s=0.5)
s3 <- wfs_ts(f1=0.5, s=0.5)

with(e1, plot(gen, p, type='l', ylim=c(0,1)))
with(e2, lines(gen, p))
with(e3, lines(gen, p))
with(s1, lines(gen, p, col='red'))
with(s2, lines(gen, p, col='red'))
with(s3, lines(gen, p, col='red'))



# plot drift example and a selection example for ne=10000, gen=11
e <- list(); for(i in 1:1000) e[[i]] <- wfs_ts(f1=0.5, s=0, gen=11, ne=10000, cnt=44)
sel <- list(); for(i in 1:1000) sel[[i]] <- wfs_ts(f1=0.5, s=0.2, gen=11, ne=10000, cnt=44)

col=rgb(0,0,0,0.05)
col2=rgb(0.5,0,0,0.05)
quartz(width=8, height=4)
# png(width=8, height=4, units='in', res=300, file='analysis/figures/wfs_sim_examples.png')
par(mai=c(0.7, 0.7, 0.15, 0.2), mfrow=c(2,1), mgp=c(1.5,0.5,0), tcl=-0.3)
with(e[[1]], plot(gen, fsamp-fsamp[1], type='l', ylim=c(0,0.7), col=col, xlab='Generation', ylab='∆Obs allele freq', main='s=0'))
for(i in 2:length(e)) with(e[[i]], lines(gen, abs(fsamp-fsamp[1]), col=col))
for(i in 1:length(e)) with(e[[i]], lines(gen, abs(p-p[1]), col=col2))
	abline(h=0.06, lty=2, col='grey')

with(sel[[1]], plot(gen, fsamp-fsamp[1], type='l', ylim=c(0,0.7), col=col, xlab='Generation', ylab='∆Obs allele freq', main='s=0.2'))
for(i in 2:length(sel)) with(sel[[i]], lines(gen, abs(fsamp-fsamp[1]), col=col))
for(i in 1:length(sel)) with(sel[[i]], lines(gen, abs(p-p[1]), col=col2))
	abline(h=0.06, lty=2, col='grey')

dev.off()