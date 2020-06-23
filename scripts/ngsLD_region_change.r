# average LD within windows and calculate change through time
# for ngsLD output (ANGSD-called loci)


######################
# calculate LD change
# run on saga
######################
#srun --ntasks=1 --mem-per-cpu=30G --time=00-02:00:00 --qos=devel --account=nn9244k --pty bash -i # need more memory
#module load R/3.6.2-foss-2019b
#R

require(data.table)

width <- 5e4; stp <- 1e4; windsz='5e4'; windnm='50kb' # window parameters. first two are numeric, second two for naming

# read in nodam2 list of loci
nodam2 <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab') # list of loci that pass nodam2 filter
setnames(nodam2, c('CHROM', 'POS', 'REF', 'ALT'))
nodam2[, posnm := paste0(CHROM, ':', POS), by = 1:nrow(nodam2)]

# read in data (gatk loci only). files from ngsLD_bypop.sh
datCan40 <- fread('analysis/ld.Can_40.gatk.gz')
datCan14 <- fread('analysis/ld.Can_14.gatk.gz')
dat07 <- fread('analysis/ld.Lof_07.gatk.gz')
dat11 <- fread('analysis/ld.Lof_11.gatk.gz')
dat14 <- fread('analysis/ld.Lof_14.gatk.gz')

# add column names
nms <- c('pos1nm', 'pos2nm', 'dist', 'r2', 'D', 'Dprime', 'r2em')
setnames(datCan14, nms)
setnames(datCan40, nms)
setnames(dat14, nms)
setnames(dat11, nms)
setnames(dat07, nms)

# trim to nodam2
nrow(datCan40)
datCan40 <- datCan40[pos1nm %in% nodam2$posnm & pos2nm %in% nodam2$posnm, ]
nrow(datCan40)

nrow(datCan14)
datCan14 <- datCan14[pos1nm %in% nodam2$posnm & pos2nm %in% nodam2$posnm, ]
nrow(datCan14)

nrow(dat07)
dat07 <- dat07[pos1nm %in% nodam2$posnm & pos2nm %in% nodam2$posnm, ]
nrow(dat07)

nrow(dat11)
dat11 <- dat11[pos1nm %in% nodam2$posnm & pos2nm %in% nodam2$posnm, ]
nrow(dat11)

nrow(dat14)
dat14 <- dat14[pos1nm %in% nodam2$posnm & pos2nm %in% nodam2$posnm, ]
nrow(dat14)

# make a chromosome column
datCan40[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
datCan14[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
dat07[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
dat11[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
dat14[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]

# make position columns
datCan40[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
datCan40[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
datCan14[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
datCan14[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
dat07[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
dat07[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
dat11[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
dat11[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
dat14[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
dat14[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]



# round pos1 and pos2 to nearest window midpoint, once for each step in width
chroms <- datCan40[, sort(unique(chr))]
for(j in 1:(width/stp)){
	nm <- paste('pos1mid', j, sep='') # create column name
	dat14[, eval(nm) := floor((pos1 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	dat11[, eval(nm) := floor((pos1 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	dat07[, eval(nm) := floor((pos1 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	datCan40[, eval(nm) := floor((pos1 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	datCan14[, eval(nm) := floor((pos1 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]

	nm <- paste('pos2mid', j, sep='') # create column name
	dat14[, eval(nm) := floor((pos2 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	dat11[, eval(nm) := floor((pos2 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	dat07[, eval(nm) := floor((pos2 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	datCan40[, eval(nm) := floor((pos2 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
	datCan14[, eval(nm) := floor((pos2 + (j-1)*stp)/width)*width + (j-1)*stp + width/2]
}
posnms <- grep('pos[[:digit:]]mid', colnames(dat14), value=TRUE) # get column names created, to check
	posnms


# calculate average LD within regions
bins14 <- dat14[!is.na(r2) & pos1mid1 == pos2mid1, .(r2ave_14 = mean(r2), r2n_14 = .N), by = .(chr, pos1mid1)]
bins11 <- dat11[!is.na(r2) & pos1mid1 == pos2mid1, .(r2ave_11 = mean(r2), r2n_11 = .N), by = .(chr, pos1mid1)]
bins07 <- dat07[!is.na(r2) & pos1mid1 == pos2mid1, .(r2ave_07 = mean(r2), r2n_07 = .N), by = .(chr, pos1mid1)] 
binsCan40 <- datCan40[!is.na(r2) & pos1mid1 == pos2mid1, .(r2ave_Can40 = mean(r2), r2n_Can40 = .N), by = .(chr, pos1mid1)] 
binsCan14 <- datCan14[!is.na(r2) & pos1mid1 == pos2mid1, .(r2ave_Can14 = mean(r2), r2n_Can14 = .N), by = .(chr, pos1mid1)] 

setkey(bins14, chr, pos1mid1)
setkey(bins11, chr, pos1mid1)
setkey(bins07, chr, pos1mid1)
setkey(binsCan40, chr, pos1mid1)
setkey(binsCan14, chr, pos1mid1)

# merge
bins <- merge(bins14, bins11, by=c('chr', 'pos1mid1'), all=TRUE)
bins <- merge(bins, bins07, by=c('chr', 'pos1mid1'), all=TRUE)
bins <- merge(bins, binsCan40, by=c('chr', 'pos1mid1'), all=TRUE)
bins <- merge(bins, binsCan14, by=c('chr', 'pos1mid1'), all=TRUE)

nrow(bins)

# add BIN_START
bins[,BIN_START:=pos1mid1 - width/2]

# calculate change
bins[, ld_diff_Lof0711 := r2ave_11 - r2ave_07]
bins[, ld_diff_Lof0714 := r2ave_14 - r2ave_07]
bins[, ld_diff_Lof1114 := r2ave_14 - r2ave_11]
bins[, ld_diff_Can := r2ave_Can14 - r2ave_Can40]



# save region means data
filenm <- paste('analysis/ld_change_region_', windsz, '_ngsLD.gatk.csv.gz', sep='')
filenm
write.csv(bins, file = gzfile(filenm), row.names = FALSE)





###########################################
# plot LD change from regions
# run on laptop
###########################################
require(data.table)
require(RColorBrewer)
require(ggplot2)

# read in
bins <- fread('analysis/ld_change_region_5e4_ngsLD.gatk.csv.gz'); width='5e4'

# add genome position
chrmax <- fread('data/lg_length.csv')
chrmax[, start := c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])]
bins <- merge(bins, chrmax[, .(chr, start)], by = 'chr')
bins[, POSgen := BIN_START + start + as.numeric(width)/2]

# long format for plotting
binsl <- melt(bins[, .(chr, POSgen, ld_diff_Can, ld_diff_Lof0711, ld_diff_Lof0714, ld_diff_Lof1114)], id.vars = c('chr', 'POSgen'), variable.name = 'pop', value.name = 'ldchange')
binsl[, pop := gsub('ld_diff_', '', pop)]

# plot
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p1 <- ggplot(binsl, aes(POSgen, ldchange, color = chr)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1, scales = 'free') +
  scale_color_manual(values = cols)
p1
ggsave(plot = p1, device = 'png', filename = paste0('figures/ld_change_vs_pos_NEA_CAN_runmean', width, '_ngsLD.gatk.png'), 
       width = 7.5, height = 6, units = 'in', dpi = 300)
