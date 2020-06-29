# Plots from SLiM soft sweep simulations
# Plot per-locus Reynolds Fst as output by slim_calcfst_reynolds.sh

require(data.table)
require(ggplot2)


##############
# read in data
##############

fstfiles <- list.files('analysis/slim_sim/', pattern = 'slim_sim_n[[:alnum:][:punct:]]*_i[[:digit:]]+.fst.csv.gz') # get the fst file names
length(fstfiles)

# read in fsts
for(i in 1:length(fstfiles)){
  if(i %% 100 == 0) cat(paste0(i, ' '))
  temp <- fread(paste0('analysis/slim_sim/',fstfiles[i]))
  
  temp[, ne := as.numeric(gsub('slim_sim_n|_s.*', '', fstfiles[i]))] # get Ne from file name
  temp[, s := as.numeric(gsub('slim_sim_n[[:digit:]]*_s|_f.*', '', fstfiles[i]))] # get s
  temp[, f := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f|_i.*', '', fstfiles[i]))] # get min f
  temp[, i := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\.fst.csv.gz', '', fstfiles[i]))] # get iteration
  
  if(i == 1) fst <- temp
  if(i > 1) fst <- rbind(fst, temp)
}
dim(fst)


##############################
# summarize across iterations
##############################

# mean fst per genome
fst[B>0, .(fst = sum(A)/sum(B)), by = .(ne, s, f, i)][, .(fstmean = mean(fst)), by = .(ne)]

# stats
fstsum1 <- fst[, .(fstmax = max(A/B, na.rm = TRUE), fstmean = sum(A)/sum(B)),  # min and max fst by simulation
              by = .(ne, s, f, i)]
fst0s <- fst[s==0 & B>0, .(fstmaxs0 = max(A/B)), by = .(ne, f)] # max FST for s=0 sims
fst2 <- merge(fst, fst0s) # add max fst from s=0 to the simulated fsts
fstsum2 <- fst2[(A/B) > fstmaxs0, .(selposmin = min(POS), 
                                    selposmax = max(POS),
                                    selwidth = (diff(range(POS)) + 1)), # width of the selected region as region where fst > maxfst from s=0
                by = .(ne, s, f, i)]
fstsum2[selposmin < 7/16*3e7 | selposmax > 9/16*3e7, ] # check for any ranges that extend outside the region in which selection was simulated


# average fst for s= 0
fst[s==0, .(fstmean = mean(WEIGHTED_FST)), by = .(ne)]

############
# plot
############

# plot FST Manhattan for two simulations
ggplot(fst[ne == 3000 & f == 0.05 & s %in% c(0, 1.5) & B>0, ], aes(x = POS, y = A/B)) +
  geom_point(size = 0.2) +
  geom_smooth() +
  facet_grid(s ~ .)
ggsave('tmp/slim_fst_persite_examples.png', width = 7, height = 4, dpi = 300)


# plot max fst vs. s, n, f
ggplot(fstsum1, aes(x = s, y = fstmax, group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Max Fst')


# plot selected width vs. s and n
ggplot(fstsum2, aes(x = s, y = selwidth/1e6, group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Width of selected region (Mb)')
ggsave('figures/slim_selwidth_persite.png', width = 7, height = 2, dpi = 150)
