require(boa) # for HPD

# Plot posterior for Ne
post_N=read.table("analysis/LOF_07_LG03_to_LOF_S_14_LG03.w_Ne_bootstrap.txt") 

quartz(width=4,height=4)
# pdf(width=4,height=4, file='analysis/figures/wfabc_Ne.pdf')
par(mai=c(0.7,1,0.3, 0.1), cex.axis=0.7, las=1, tcl=-0.2, mgp=c(3.5,0.6,0))
plot(density(post_N[,1]),lwd=2,main="Ne Posterior",xlab="", xlim=c(1,1e5), log='x')
mtext("Ne",1,line=2)

dev.off()


