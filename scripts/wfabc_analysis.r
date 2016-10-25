require(boa) # for HPD

# Plot posterior for Ne
post_N=read.table("analysis/LOF_07_LG03_to_LOF_S_14_LG03.w_Ne_bootstrap.txt") 

quartz(width=4,height=4)
# pdf(width=4,height=4, file='analysis/figures/wfabc_Ne.pdf')
par(mai=c(0.7,1,0.3, 0.1), cex.axis=0.7, las=1, tcl=-0.2, mgp=c(3.5,0.6,0))
plot(density(post_N[,1]),lwd=2,main="Ne Posterior",xlab="", xlim=c(1,1e5), log='x')
mtext("Ne",1,line=2)

dev.off()


# Analyze posterior for s for selected loci
post_s=read.table("analysis/LOF_07_LG03_to_LOF_S_14_LG03.w_posterior_s.txt") 
s_hpd <- apply(post_s,1,boa.hpd,1-0.95) # calc 95% HPD for s at each locus
posinds <- which(s_hpd[1,]>0 & s_hpd[2,]>0) # evidence for positive selection: 113
neginds <- which(s_hpd[1,]<0 & s_hpd[2,]<0) # evidence for positive selection: 27

# Plot posterior for s for selected loci
post_s=read.table("analysis/LOF_07_LG03_to_LOF_S_14_LG03.w_posterior_s.txt") 

inds = c(5,1)
xlims=c(-0.6, 0.6)
quartz(height=6, width=3)
# pdf(height=6, width=3, file='analysis/figures/wfabc_s_examples.pdf')
par(mfrow=c(2,1), mai=c(0.7, 0.7, 0.2, 0.1), las=1, mgp=c(2.5, 0.7, 0), tcl=-0.2)
plot(density(t(post_s[inds[1],])),lwd=2,main=paste("s Posterior (locus ", inds[1], ")", sep=''), xlab="s", xlim=xlims)
	abline(v=0, col='grey')
plot(density(t(post_s[inds[2],])),lwd=2,main=paste("s Posterior (locus ", inds[2], ")", sep=''), xlab="s", xlim=xlims) 
	abline(v=0, col='grey')

dev.off()

# Plot posterior for Ne for selected loci: are they different?
post_Nloc=read.table("analysis/LOF_07_LG03_to_LOF_S_14_LG03.w_posterior_N.txt") 

inds = c(5,1)
xlims=c(2000, 5e4)
quartz(height=6, width=3)
# pdf(height=6, width=3, file='analysis/figures/wfabc_s_examples.pdf')
par(mfrow=c(2,1), mai=c(0.7, 1, 0.2, 0.1), las=1, mgp=c(3.5, 0.7, 0), tcl=-0.2, cex.axis=0.7)
plot(density(t(post_Nloc[inds[1],])),lwd=2,main=paste("N Posterior (locus ", inds[1], ")", sep=''), xlab='', xlim=xlims, log='x')
	mtext("Ne",1,line=2)
	abline(v=0, col='grey')
plot(density(t(post_Nloc[inds[2],])),lwd=2,main=paste("N Posterior (locus ", inds[1], ")", sep=''), xlab='', xlim=xlims, log='x')
	mtext("Ne",1,line=2)
	abline(v=0, col='grey')

dev.off()
