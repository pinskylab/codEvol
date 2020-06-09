# Plots from SLiM soft sweep simulations
# Plot pcangsd outlier calculations

require(data.table)
library(RcppCNPy)
require(ggplot2)

# need to have run code in slim_calcpcangsdoutlier.sh

##############
# read in data
##############

files <- list.files('analysis/slim_sim/', pattern = 'slim_sim_n[[:alnum:][:punct:]]*.selection.npy') # get the fst file names
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

  # calculate p-value and FDR correct (whole genome)
  temp[, p := pchisq(q = stat, df = 1, lower.tail = FALSE)]
  temp[, pfdr := p.adjust(p, method = 'fdr')]

  # FDR correct (only chrom with selection)
  temp[CHROM==1, pfdr_chr1 := p.adjust(p, method = 'fdr')]

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
  temp[, i := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\.selection.npy', '', files[i]))] # get iteration
  
  npctemp[, ne := temp$ne[1]]
  npctemp[, s := temp$s[1]]
  npctemp[, f := temp$f[1]]
  npctemp[, i := temp$i[1]]
  
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
              by = .(ne, s, f, i)]

sum2 <- sum1[, .(prop = sum(minpfdr <= 0.05)/.N), by = .(ne, s, f)] # proportion of sims that detected at least one SNP under selection
############
# plot
############

# quick plot of each sim. now too many to visualize well.
p <- ggplot(outl, aes(x=POS, y= -log10(pfdr))) + 
  geom_point(data = outl[-log10(pfdr) >= -log10(0.05),], size = 0.2, shape = 1, alpha = 0.3, color = 'red') +
  geom_point(data = outl[-log10(pfdr) < -log10(0.05),], size = 0.2, shape = 1, alpha = 0.3, color = 'black') +
  facet_grid(ne + i ~ s + f)
  #geom_hline(yintercept= -log10(0.05), linetype="dashed")
ggsave('figures/slim_pcangsdoutlier.png', plot = p, width = 8, height = 32, dpi = 300)



# plot min p vs. s and ne (group by f)
ggplot(sum1, aes(x = s, y = -log10(minp), group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Min -log10 p-value')
ggsave('figures/slim_pcangsdoutlier_minp.png', width = 7, height = 2, dpi = 150)

# plot min pfdr vs. s and ne (group by f)
ggplot(sum1, aes(x = s, y = -log10(minpfdr), group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Min -log10\nFDR-corrected p-value') + 
  geom_hline(yintercept= -log10(0.05), linetype="dashed")
ggsave('figures/slim_pcangsdoutlier_minpfdr.png', width = 7, height = 2, dpi = 150)


# plot number of pfdr<0.05 vs. s and ne (group by f)
ggplot(sum1, aes(x = s, y = log10(npfdrlow+1), group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
#  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('log10(# of pfdr<0.05 + 1')
ggsave('figures/slim_pcangsdoutlier_npfdrlow.png', width = 7, height = 2, dpi = 150)

# plot proportion of sims that have at least one pfdr<0.05 vs. s and ne (group by f)
ggplot(sum2, aes(x = s, y = prop, group = f, color = as.factor(f))) +
  geom_point(data = sum2[f == 0.05,], size = 3) +
  geom_point(data = sum2[f == 0.2,], size = 1) +
#  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Proportion min FDR-corrected p < 0.05')
ggsave('figures/slim_pcangsdoutlier_propminpfdrlow.png', width = 7, height = 2, dpi = 150)