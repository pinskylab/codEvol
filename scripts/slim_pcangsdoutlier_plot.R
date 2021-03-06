# Plots pcangsd outlier calculations from SLiM soft sweep simulations
# see end of script for plotting PCAs

require(data.table)
library(RcppCNPy)
require(ggplot2)

# need to have run code in slim_calcpcangsdoutlier.sh

##############
# read in data
##############

files <- list.files('analysis/slim_sim/', pattern = 'slim_sim_n[[:alnum:][:punct:]]*.selection.npy') # get the pcangsd file names
length(files)

for(i in 1:length(files)){
  if(i %% 100 == 0) cat(paste0(i, ' '))
  # read in selection calculations and sites
  temp <- as.data.table(npyLoad(paste0('analysis/slim_sim/',  files[i]))) # selection statistics along a tested PC; χ²-distributed with 1 degree of freedom
  sitestemp <- fread(paste0('analysis/slim_sim/',  gsub('.selection.npy', '.sites', files[i])), header = FALSE) # sites
  
  # number of PCs selected by MAP test in pcangsd
  npctemp <- data.table(npc = ncol(temp))
  
  # trim to first pc
  temp <- temp[, 1] # trim to first pc
  setnames(temp, 1, 'stat')

  # add site information
  sitestemp[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)] # split chromosome and position
  temp <- cbind(temp, sitestemp[, .(CHROM, POS)])

  # calculate p-value and FDR correct
  temp[, p := pchisq(q = stat, df = 1, lower.tail = FALSE)]
  temp[, pfdr := p.adjust(p, method = 'fdr')]

  # calc num sites, trim to p<0.001 (or at least one site with the min FDR)
  temp[, npos := .N]
  keepinds <- temp[, which(p < 0.001)]
  if(length(keepinds)<1) keepinds <- temp[, which.min(pfdr)] # at least keep the min fdr
  temp <- temp[keepinds, ]
  sitestemp <- sitestemp[keepinds,]

  
  # add sim metadata
  temp[, ne := as.numeric(gsub('slim_sim_n|_s.*', '', files[i]))] # get Ne from file name
  temp[, s := as.numeric(gsub('slim_sim_n[[:digit:]]*_s|_f.*', '', files[i]))] # get s
  temp[, f := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f|_i.*', '', files[i]))] # get min f
  temp[, i := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\.selection.npy|\\_comb.selection.npy', '', files[i]))] # get iteration
  combval <- ifelse(grepl('comb', files[i]), 1, 0) # not sure why, but this fails if inside the data.table expression
  temp[, comb := combval] # get whether chromsome has been combined with 20 s=0 chromosomes
  
  npctemp[, ne := temp$ne[1]]
  npctemp[, s := temp$s[1]]
  npctemp[, f := temp$f[1]]
  npctemp[, i := temp$i[1]]
  npctemp[, comb := temp$comb[1]]
  
  if(i == 1){
    outl <- temp
    npc <- npctemp
  } 
  if(i > 1){
    outl <- rbind(outl, temp)
    npc <- rbind(npc, npctemp)
  } 
}
dim(outl)



##############################
# summarize across iterations
##############################

# stats on the region under selection
sum1 <- outl[, .(minp = min(p), minpfdr = min(pfdr), npfdrlow = sum(pfdr < 0.05), npos = unique(npos)),  # min p and pfdr
              by = .(ne, s, f, i, comb)]

sum2 <- sum1[, .(prop = sum(minpfdr <= 0.05)/.N, npos = mean(npos)), by = .(ne, s, f, comb)] # proportion of sims that detected at least one SNP under selection


# write out
write.csv(sum2, file = gzfile('analysis/slim_pcangsd.summary.csv.gz'), row.names = FALSE)

##########################
# plot outlier detection
##########################

# plot of each sim. now too many to visualize well. VERY SLOW (hours)
# p <- ggplot(outl, aes(x=POS, y= -log10(pfdr))) + 
#   geom_point(data = outl[-log10(pfdr) >= -log10(0.05),], size = 0.2, shape = 1, alpha = 0.3, color = 'red') +
#   geom_point(data = outl[-log10(pfdr) < -log10(0.05),], size = 0.2, shape = 1, alpha = 0.3, color = 'black') +
#   facet_grid(ne + i ~ s + f + comb)
#   #geom_hline(yintercept= -log10(0.05), linetype="dashed")
# ggsave('figures/slim_pcangsdoutlier.png', plot = p, width = 16, height = 32, dpi = 300)

# number of SNPs
ggplot(sum1, aes(x = ne, y = npos, group = as.factor(f), color = as.factor(f))) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  facet_grid(~ comb, scales = 'free') +
  labs(color = 'Initial frequency') +
  scale_y_log10() +
  scale_x_log10() +
  ylab('# NSPs')
ggsave('figures/slim_nsnps.png', width = 6, height = 3, dpi = 150)

# plot min p vs. s and ne and comb (group by f)
ggplot(sum1, aes(x = s, y = -log10(minp), group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(comb ~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Min -log10 p-value')
ggsave('figures/slim_pcangsdoutlier_minp.png', width = 7, height = 4, dpi = 150)

# plot min p by burn-in file
ggplot(sum1, aes(x = as.factor(i), y = -log10(minp))) +
  geom_boxplot() +
  facet_grid(~ne)
  labs(color = 'Burn-in file') +
  ylab('Min -log10 p-value')

# plot min pfdr vs. s and ne and comb (group by f)
ggplot(sum1, aes(x = s, y = -log10(minpfdr), group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(comb ~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Min -log10\nFDR-corrected p-value') + 
  geom_hline(yintercept= -log10(0.05), linetype="dashed")
ggsave('figures/slim_pcangsdoutlier_minpfdr.png', width = 7, height = 4, dpi = 150)


# plot number of pfdr<0.05 vs. s and ne (group by f)
ggplot(sum1, aes(x = s, y = npfdrlow+1, group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
#  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(comb ~ ne) +
  labs(color = 'Initial frequency') +
  scale_y_log10() +
  ylab('# of pfdr<0.05 + 1')
ggsave('figures/slim_pcangsdoutlier_npfdrlow.png', width = 7, height = 4, dpi = 150)

# plot proportion of sims that have at least one pfdr<0.05 vs. s and ne (group by f)
ggplot(sum2, aes(x = s, y = prop, group = f, color = as.factor(f))) +
  geom_point(data = sum2[f == 0.05,], size = 3) +
  geom_point(data = sum2[f == 0.2,], size = 1) +
#  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(comb ~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Proportion min FDR-corrected p < 0.05')
ggsave('figures/slim_pcangsdoutlier_propminpfdrlow.png', width = 7, height = 4, dpi = 150)

# plot proportion of sims that have at least one pfdr<0.05 vs. npos (group by f)
ggplot(sum2, aes(x = npos, y = prop, group = f, color = s)) +
  geom_point() +
  #facet_grid(~s) +
  labs(color = 's') +
  scale_x_log10() +
  ylab('Proportion min FDR-corrected p < 0.05')


#########################
# Plot PCA from SLiM sims
#########################

# set up colors. first 22 are 1907, latter 24 are 2011
cols <- c(rep('#af8dc3', 22), rep('#7fbf7b', 24))

# pick which sim to plot
ne=300
s=0
f=0.05
i=1
comb='_comb' # or ''

cov <- fread(paste0('analysis/slim_sim/slim_sim_n', ne, '_s', s, '_f', f, '_i', i, comb, '.cov'))
eig <- eigen(cov)

png(filename = paste0('figures/slim_pca_n', ne, '_s', s, '_f', f, '_i', i, 
                      ifelse(comb=='_comb', '_wholegen.png', 'onechrom.png')),
    width = 4, height = 4, units = 'in', res = 300)
plot(eig$vectors[, 1:2], xlab = 'PC1', ylab = 'PC2', col = cols,
     main = paste0('ne=', ne, ' s=', s, ' init freq=', f, '\niter#', i, ifelse(comb=='_comb', ' whole genome', ' one chromosome')))
legend('topright', legend = c('1', '11'), col = cols[c(1,23)], pch = 1, title = 'generation')
dev.off()

