# Plots from SLiM soft sweep simulations
# Plot Fst from sliding-window

require(data.table)

# need to have run code in slim_run.R and followed instructions there to calculat sliding window Fst as well

############
# plot
############
require(ggplot2)

fstfiles <- list.files('analysis/', pattern = 'slim_sim_n[[:alnum:][:punct:]]*.fst') # get the fst file names
length(fstfiles)

# read in fsts
for(i in 1:length(fstfiles)){
  temp <- fread(paste0('analysis/',fstfiles[i]))
  
  temp[, ne := as.numeric(gsub('slim_sim_n|_s.*', '', fstfiles[i]))] # get Ne from file name
  temp[, s := as.numeric(gsub('slim_sim_n[[:digit:]]*_s|_f.*', '', fstfiles[i]))] # get s
  temp[, f := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f|_i.*', '', fstfiles[i]))] # get min f
  temp[, i := as.numeric(gsub('slim_sim_n[[:digit:]]*_s[[:digit:].]*_f0\\.[[:digit:]]*_i|\\.windowed.weir.fst', '', fstfiles[i]))] # get iteration
  
  if(i == 1) fst <- temp
  if(i > 1) fst <- rbind(fst, temp)
}
dim(fst)

fst[, mid := rowMeans(cbind(BIN_START, BIN_END))]


# quick plot
p <- ggplot(fst, aes(x=mid, y=WEIGHTED_FST)) +
  geom_line(alpha = 0.5, size = 0.2) +
  facet_grid(ne + i ~ s)
ggsave('figures/slim_fst.png', plot = p, width = 8, height = 16, dpi = 300)


