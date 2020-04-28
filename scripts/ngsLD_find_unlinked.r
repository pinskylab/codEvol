# identify unlinked loci from ngsLD
# run after nsdLD_find_blocks.sh
# currently only set up for GATK loci


# functions
require(data.table)

findmid <- function(POS){ # function to return a value near the middle of a vector of positions
	mn <- mean(range(POS))
	return(POS[which.min(abs(POS - mn))])
}


# read in list of linkage blocks
ld <- fread('analysis/ld.blocks.gatk.csv.gz', drop = 1) # linkage blocks from ngsLD_find_blocks.r

# find a locus near the middle of each block to keep
canclust <- ld[!is.na(cluster_can), .(CHROM = unique(CHROM), POS = findmid(POS), nloci = length(POS)), by = .(cluster = cluster_can)] # find the cluster midpoints
canclust <- rbind(canclust, ld[is.na(cluster_can), .(CHROM, POS, nloci = 1, cluster = NA_integer_)]) # add the loci not in blocks
setkey(canclust, CHROM, POS)
setcolorder(canclust, c('CHROM', 'POS', 'cluster', 'nloci'))

lofclust <- ld[!is.na(cluster_lof), .(CHROM = unique(CHROM), POS = findmid(POS), nloci = length(POS)), by = .(cluster = cluster_lof)] # find the cluster midpoints
lofclust <- rbind(lofclust, ld[is.na(cluster_lof), .(CHROM, POS, nloci = 1, cluster = NA_integer_)]) # add the loci not in blocks
setkey(lofclust, CHROM, POS)
setcolorder(lofclust, c('CHROM', 'POS', 'cluster', 'nloci'))


# compare
nrow(ld)
nrow(canclust)
nrow(lofclust)

# write out loci trimmed to unlinked blocks. no quotes so can process easily on command line
write.csv(canclust, gzfile('analysis/ld.unlinked.Can.gatk.csv.gz'), quote = FALSE)
write.csv(lofclust, gzfile('analysis/ld.unlinked.Lof.gatk.csv.gz'), quote = FALSE)
