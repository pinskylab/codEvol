# Plots from SLiM soft sweep simulations
# Plot Fst from sliding-window

require(data.table)
require(ggplot2)

# need to have run code in slim_run.R and followed instructions there to calculat sliding window Fst as well

##############
# read in data
##############

fstfiles <- list.files('analysis/slim_sim/', pattern = 'slim_sim_n[[:alnum:][:punct:]]*.fst') # get the fst file names
length(fstfiles)

# read in fsts
for(i in 1:length(fstfiles)){
  if(i %% 100 == 0) cat(paste0(i, ' '))
  temp <- fread(paste0('analysis/slim_sim/',fstfiles[i]))
  
  temp[, ne := as.numeric(gsub('slim_sim_n|_s.*', '', fstfiles[i]))] # get Ne from file name
  temp[, s := as.numeric(gsub('slim_sim_n[[:digit:]]*_s|_f.*', '', fstfiles[i]))] # get s
  temp[, f := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f|_i.*', '', fstfiles[i]))] # get min f
  temp[, i := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\.windowed.weir.fst', '', fstfiles[i]))] # get iteration
  
  if(i == 1) fst <- temp
  if(i > 1) fst <- rbind(fst, temp)
}
dim(fst)

fst[, mid := rowMeans(cbind(BIN_START, BIN_END))]

##############################
# summarize across iterations
##############################

# an average Manhattan plot
fstave <- fst[, .(fst = mean(WEIGHTED_FST), fstsd = sd(WEIGHTED_FST)), by = .(mid, ne, s, f)]


# stats on the region under selection
fst[, selection_region := 0]
fst[mid > 7/16*3e7 & mid < 9/16*3e7, selection_region := 1] # label the region in which selection was simulated
fstsum1 <- fst[, .(fstmax = max(WEIGHTED_FST), fstmean = mean(WEIGHTED_FST)),  # min and max fst in and outside selected region
              by = .(selection_region, ne, s, f, i)]
fst2 <- merge(fst, fst[s==0, .(fstmaxs0 = max(WEIGHTED_FST)),
                       by = .(ne, f)]) # add max fst from s=0 to the simulated fsts
fstsum2 <- fst2[WEIGHTED_FST > fstmaxs0, .(selposmin = min(mid), 
                                         selposmax = max(mid),
                                         selwidth = (diff(range(mid)) + 1)/1e6), # width of the selected region as region where fst > maxfst from s=0
              by = .(ne, s, f, i)]
fstsum2[selposmin < 7/16*3e7 | selposmax > 9/16*3e7, ] # check for any ranges that extend outside the region in which selection was simulated
fstsum3 <- dcast(fstsum1, ne + s + f + i ~ selection_region, value.var = c('fstmax', 'fstmean')) # prep to calculat in/out ratios
fstsum3[, ':='(fstmean_ratio = fstmean_1/fstmean_0, fstmax_ratio = fstmax_1/fstmax_0), by = .(ne, s, f, i)]
  
fstsum <- merge(fstsum2[, .(ne, s, f, i, selwidth)], 
                fstsum3[, .(ne, s, f, i, fstmax_sel = fstmax_1, fstmax_nosel = fstmax_0, fstmean_sel = fstmean_1, fstmean_nosel = fstmean_0,
                            fstmean_ratio, fstmax_ratio)],
                all.y = TRUE, by = c('ne', 's', 'f', 'i'))
fstsum[is.na(selwidth), selwidth := 0] # set selwidth to 0 where no region existed with fst>max in nonselected region



############
# plot
############

# quick plot of each sim
p <- ggplot(fst, aes(x=mid, y=WEIGHTED_FST)) +
  geom_point(size = 0.2, shape = 1, alpha = 0.3) +
  facet_grid(ne + i ~ s + f, scales = 'free')
ggsave('figures/slim_fst.png', plot = p, width = 8, height = 24, dpi = 300)


# plot of average pattern
p2 <- ggplot(fstave, aes(x = mid, y = fst, ymin = fst - fstsd, ymax = fst + fstsd)) +
  annotate('rect', xmin = 7/16*3e7, xmax = 9/16*3e7, ymin = 0, ymax = Inf, fill = 'red', color = NA, alpha = 0.1) + # for the region under selection
  geom_ribbon(fill = 'grey', color = NA) +
  geom_point(size = 0.2, shape = 1, alpha = 0.3) +
  facet_grid(ne ~ s + f, scales = 'free')
ggsave('figures/slim_fstave.png', plot = p2, width = 8, height = 4, dpi = 300)

# plot max fst vs. ne and s
ggplot(fst[, .(maxfst = max(WEIGHTED_FST)), by = .(ne,s,i,f)], aes(s, y = maxfst, group = as.factor(f), color = as.factor(f))) +
  geom_point() +
  facet_grid(~ne)

# plot fst ratio vs. s and n
ggplot(fstsum, aes(x = s, y = fstmax_ratio, group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Ratio of max Fst inside vs. \noutside selected region')
ggsave('figures/slim_fstmaxratio.png', width = 7, height = 2, dpi = 150)


# plot selected width vs. s and n
ggplot(fstsum, aes(x = s, y = selwidth, group = as.factor(f), color = as.factor(f))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', size = 0.5) +
  facet_grid(~ ne) +
  labs(color = 'Initial frequency') +
  ylab('Width of selected region (Mb)')
ggsave('figures/slim_selwidth.png', width = 7, height = 2, dpi = 150)
