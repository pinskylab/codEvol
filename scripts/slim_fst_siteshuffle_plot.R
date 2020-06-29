## examine results from site reshuffled FST from SLiM genotypes
## run after slim_fst_siteshuffle_null.sh

#################
# parameters
#################
minloci <- 2 # should match slim_fst_siteshuffle_null.r
winsz <- 50000 # window size
winstp <- 10000 # window step

###########################
# load functions
###########################
require(data.table)
#require(plyr)
require(ggplot2)
require(RColorBrewer)

calcp <- function(fst, null) return((sum(null > fst)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen

#############################
# read in and summarize data
#############################

files <- list.files('analysis/slim_sim/', pattern = 'slim_sim_n[[:alnum:][:punct:]]*.fst.siteshuffle.csv') # get the max fst file names from site shuffling
length(files)

for(i in 1:length(files)){
  if(i %% 100 == 0) cat(paste0(i, ' '))
  
  # max FST per genome from reshuffling
  null <- fread(paste0('analysis/slim_sim/', files[i]))

  # sliding window FST A and B components from ANGSD, after collapsing to unlinked loci
  # header is missing the fst column, so have to skip and make our own
  # need to make the all loci AB files
  fstAB <- fread(paste0('analysis/slim_sim/', gsub('.siteshuffle|_comb', '', files[i]))) # fst components

  # create new columns as indices for windows
  for(j in 1:(winsz/winstp)){
    fstAB[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
  }
  
  # calc fst and # snps per window
  for(j in 1:(winsz/winstp)){
    if(j ==1){
      fstwin <- fstAB[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
    } 
    if(j > 1){
      fstwin <- rbind(fstwin, fstAB[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
    } 
  }
  
  # calc p-values per window where >=minloci
  fstwin[nloci >= minloci, p := calcp(fst, null$x), by = .(CHROM, midPos)]

  # summarize
  tempsum <- fstwin[, .(minp = min(p, na.rm = TRUE), npl05 = sum(p < 0.05, na.rm = TRUE), nloci = sum(nloci), nwin = .N)]
  tempsum[, ne := as.numeric(gsub('slim_sim_n|_s.*', '', files[i]))] # get Ne from file name
  tempsum[, s := as.numeric(gsub('slim_sim_n[[:digit:]]*_s|_f.*', '', files[i]))] # get s
  tempsum[, f := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f|_i.*', '', files[i]))] # get min f
  tempsum[, i := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\.fst.siteshuffle.csv.gz|\\_comb.fst.siteshuffle.csv.gz', '', files[i]))] # get iteration
  combval <- ifelse(grepl('comb', files[i]), 1, 0) # not sure why, but this fails if inside the data.table expression
  tempsum[, comb := combval] # get whether chromsome has been combined with 20 s=0 chromosomes
  
  
  if(i == 1){
    sum <- tempsum
  } else {
    sum <- rbind(sum, tempsum)
  }
}

# summarize across iterations
# fpl05: fraction of sims with at least one window p<0.05
# f2pl05: fraction of sims with at least two windows p<0.05
sum2 <- sum[, .(fpl05 = sum(minp < 0.05)/.N, 
                f2pl05 = sum(npl05 > 1)/.N,
                f3pl05 = sum(npl05 > 2)/.N), 
            by = .(ne, s, f, comb)]

# write out
write.csv(sum2[, .(ne, s, f, comb, fpl05)], file = gzfile('analysis/slim_fst_siteshuffle.summary.csv.gz'), row.names = FALSE)

#########################################
# plots of FST calculation summaries
# see next section for a Manhattan plot
#########################################

# min p
ggplot(sum, aes(s, -log10(minp), group = f, color = as.factor(f))) +
  geom_point() +
  geom_smooth(method ='lm') +
  facet_grid(comb ~ ne)
ggsave('figures/slim_fst_siteshuffle_minp.png', width = 7, height = 4, dpi = 150)

# number of windows p<0.05
ggplot(sum, aes(s, npl05, group = f, color = as.factor(f))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(comb ~ ne)
ggsave('figures/slim_fst_siteshuffle_npl05.png', width = 7, height = 4, dpi = 150)

# fraction of sims with at least one window p<0.05
ggplot(sum2, aes(s, fpl05, group = f, color = as.factor(f))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(comb ~ ne) +
  coord_cartesian(ylim = c(0,1))
ggsave('figures/slim_fst_siteshuffle_fpl05.png', width = 7, height = 4, dpi = 150)

# fraction of sims with at least one window p<0.05
ggplot(sum2, aes(s, f2pl05, group = f, color = as.factor(f))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(comb ~ ne) +
  coord_cartesian(ylim = c(0,1))

# fraction of sims with at least two windows p<0.05
ggplot(sum2, aes(s, f3pl05, group = f, color = as.factor(f))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(comb ~ ne) +
  coord_cartesian(ylim = c(0,1))

################
# Manhattan plot 
################
# pick which sim to plot
ne=3000
s=1.5
f=0.05
i=4

thisfstAB <- fread(paste0('analysis/slim_sim/slim_sim_n', ne, '_s', s, '_f', f, '_i', i, '.fst.csv.gz'))

# create new columns as indices for windows
for(j in 1:(winsz/winstp)){
  thisfstAB[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
}

# calc fst and # snps per window
for(j in 1:(winsz/winstp)){
  if(j ==1){
    thisfstwin <- thisfstAB[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
  } 
  if(j > 1){
    thisfstwin <- rbind(thisfstwin, thisfstAB[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
  } 
}

thisfstwin[, plot(midPos, fst, xlab = 'POS', ylab = 'FST', main = paste0('ne=', ne, ' s=', s, ' init freq=', f, ' iter#', i))]




