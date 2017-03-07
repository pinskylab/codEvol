# compare rejection vs. regression ABC

sampsize <- '50,48' # choose the set of simulations

filesregr <- list.files(path='analysis/temp', pattern=paste('wfs_abcregr_sampsize', sampsize, sep=''), full.names=TRUE)

posts <- read.csv(gzfile(filesregr[1]))
