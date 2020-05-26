# Plot PCA from PCAngsd
# input files created on saga with angsd_pcangsd_pca.sh
# srun --ntasks=1 --mem-per-cpu=100G --time=00-02:00:00 --qos=devel --account=nn9244k --pty bash -i	
# module load PCAngsd/200115-foss-2019a-Python-2.7.15	
# pcangsd.py -beagle data_31_01_20/All_ind_beagle.GATK.gz -minMaf 0.05 -threads 16 -o analysis/pcangsd_pca_gatk # gatk loci a little surprised it ran with 16 threads. or maybe it didn't use them
# pcangsd.py -beagle data_31_01_20/All_ind_beagle.GATK_no_dam.gz -minMaf 0.05 -threads 16 -o analysis/pcangsd_pca_gatk_no_dam # gatk loci no adna damage
# pcangsd.py -beagle data_31_01_20/All_ind_beagle.GATK_ex_inv.gz -minMaf 0.05 -threads 16 -o analysis/pcangsd_pca_gatk_ex_inv # gatk loci no inversion LGs


library(data.table)
library(RColorBrewer)

pcagt <- fread('analysis/pcangsd_pca_gatk.cov') # only GATK snps
pcagtnd <- fread('analysis/pcangsd_pca_gatk_no_dam.cov') # only GATK no aDNA damage snps
pcagtei <- fread('analysis/pcangsd_pca_gatk_ex_inv.cov') # only GATK w/out inversions
pcagtul <- fread('analysis/pcangsd_pca_gatk_unlink.cov') # only GATK snps unlinked
pcagtndul <- fread('analysis/pcangsd_pca_gatk_no_dam_unlink.cov') # only GATK no aDNA damage snps unlinked
pcagteiul <- fread('analysis/pcangsd_pca_gatk_ex_inv_unlink.cov') # only GATK w/out inversions unlinked
pops <- fread('data_31_01_20/List_to_malin.tab')

# eigenvector decomposition of the covariance matrix
egt <- eigen(pcagt)
egtnd <- eigen(pcagtnd)
egtei <- eigen(pcagtei)
egtul <- eigen(pcagtul)
egtndul <- eigen(pcagtndul)
egteiul <- eigen(pcagteiul)

# make colors for populations
cols <- brewer.pal(6, 'RdYlBu')
pops[, unique(Pop)]
pops[Pop == 'Historic_Canada', col := cols[1]]
pops[Pop == 'Canada_2013', col := cols[2]]
pops[Pop == 'Historic_Lofoten', col := cols[6]]
pops[Pop == 'Lofoten_2011', col := cols[5]]
pops[Pop == 'Lofoten_2014', col := cols[4]]

# plot PCAs
# red is early, black is late
pdf('figures/pcangsd_pca_all.pdf', width = 6, height = 4)
par(mfrow=c(2,3), mai = c(0.6, 0.7, 0.3, 0.1), las = 1, mgp = c(2, 0.7, 0))

plot(egt$vectors[, 1:2], xlab = 'PC1', ylab = 'PC2', col = pops$col,
     main = 'GATK')

plot(egtnd$vectors[, 1:2], xlab = 'PC1', ylab = 'PC2', col = pops$col,
     main = 'GATK no dam')

plot(egtei$vectors[, 1:2], xlab = 'PC1', ylab = 'PC2', col = pops$col,
     main = 'GATK ex inv')

legend('top', pch = 1, col = pops[!duplicated(Pop), col], legend = pops[!duplicated(Pop), Pop], cex = 0.5)

plot(egtul$vectors[, 1:2], xlab = 'PC1', ylab = 'PC2', col = pops$col,
     main = 'GATK unlinked')

plot(egtndul$vectors[, 1:2], xlab = 'PC1', ylab = 'PC2', col = pops$col,
     main = 'GATK no dam unlinked')

plot(egteiul$vectors[, 1:2], xlab = 'PC1', ylab = 'PC2', col = pops$col,
     main = 'GATK ex inv unlinked')


dev.off()
