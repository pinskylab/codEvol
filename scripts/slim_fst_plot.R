# Make SLiM soft sweep simulations
# then plot Fst from sliding-window

require(data.table)


###########
# functions
###########

runslim <- function(L=30000000, ne=50, f=0.1, n1=22, n2=22, s=1.0, r=3.11e-8, g=11, o='tmp/test', i='tmp/test.vcf'){
  # notice complicated quoting needed to pass a string to slim
  # use -long to enable verbose output for troubleshooting
  system2("slim", c(paste0("-d ", c(paste0("L=", L), paste0("ne=", ne), paste0("f=", f), paste0("n1=", n1), 
                                    paste0("n2=", n2), paste0("s=", s), paste0("r=", r), 
                                    paste0("g=", g), paste0("\"o='", o, "'\""), paste0("\"i='", i, "'\""))), 
                    '-long', " scripts/slim_softsweep.slim"))
}

###########
# run slim
###########

bifiles <- list.files('analysis/', pattern = 'slim_burnin_n') # get the burn-in files
nes <- as.numeric(gsub('slim_burnin_n|_[12345]\\.vcf', '', bifiles)) # find the ne values from the file names
iter <- as.numeric(gsub('slim_burnin_n[[:digit:]]*_|\\.vcf', '', bifiles)) # find the iteraction #s from the file names

# Lofoten 07-11 parameters
paramtable <- expand.grid(iter = sort(unique(iter)), ne = sort(unique(nes)), f = c(0.1, 0.01), n1 = 22, n2 = 23, s = c(0.01, 0.03, 0.1, 0.3, 1), 
                          r = 3.11e-8, g = 11, L = 3e7)
dim(paramtable)

# run slim for Lofoten
for(i in 1:nrow(paramtable)){
  infile = with(paramtable[i,], paste0('analysis/slim_burnin_n', ne, '_', iter, '.vcf'))
  outfile = with(paramtable[i,], paste0('analysis/slim_sim_n', ne, '_s', s, '_f', f, '_i', iter))
  with(paramtable[i,], runslim(L = L, ne = ne, f = f, n1 = n1, n2 = n2, s = s, r = r, g = g, o = outfile, i = infile))
  
}



#############
# calculate fst
#############


#rsync -tz analysis/slim_sim_n*.vcf USER@saga.sigma2.no:~/historic_malin/analysis/ # copy to saga
# on saga:
#sbatch scripts/slim_calcfst.sh # calc fst with vcftools
#rsync -tz USER@saga.sigma2.no:~/historic_malin/analysis/slim_sim_n*.fst analysis/ # copy back from saga


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


