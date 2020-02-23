# plot pFst output

require(ggplot2)
require(data.table)

# read in file
dat<-fread('gunzip -c analysis/pfst_can.out.gz', header=FALSE ); nm <- 'can'

# prep for plotting
nrow(dat)
dat[, indx := 1:length(dat$V2)] # add a snp index (could make this genome position instead)
dat<-dat[dat$V3 < 0.9,] # trim out high p-values
nrow(dat)

# make a nucleotide position for the whole genome (start position for each chr)
chrmax <- dat[,.(len=max(V2)), by=V1]
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, V1)
setkey(dat, V1)
dat <- merge(dat, chrmax[,.(V1, start)], by='V1')
dat[,POSgen := V2 + start]
dat[,start := NULL]

# fdr-adjust the p-values
dat[, pfdr := p.adjust(V3, method = 'fdr')]

# add a vector for color by LG
cols <- c('#a6cee333', '#1f78b433') # light blue, blue, partially transparent: for alternating LGs
dat[, lgcol := cols[1]]
dat[V1 %in% chrmax$V1[seq(2, nrow(chrmax),by=2)], lgcol := cols[2]]

# plot and save out
png(filename = paste0('figures/pFst_', nm, '.png'), width = 20, height = 4, units = "in", res = 300)
dat[, plot(POSgen/1e6, -log10(pfdr), xlab = 'Position (Mb)', ylab = "-log10(pFst FDR-adjusted)",
           col = lgcol)]

dev.off()
