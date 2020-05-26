# identify unlinked loci from ngsLD
# run after nsdLD_find_blocks.sh
# currently only set up for GATK loci


# functions
require(data.table)

findmid <- function(POS){ # function to return a value near the middle of a vector of positions
	mn <- mean(range(POS))
	return(POS[which.min(abs(POS - mn))])
}


# read in list of linkage blocks from ngsLD (out to 5k)
ld <- fread('analysis/ld.blocks.gatk.csv.gz', drop = 1) # linkage blocks from ngsLD_find_blocks.r

# read in inversion coords
inv <- fread('data/inversions.csv')

# mark the inversions
for(i in 1:nrow(inv)){
  ld[CHROM == inv$CHROM[i] & POS >= inv$POSstart[i] & POS <= inv$POSend[i], cluster_can := -i] # mark as cluster -i in CAN
  ld[CHROM == inv$CHROM[i] & POS >= inv$POSstart[i] & POS <= inv$POSend[i], cluster_lof := -i] # in LOF
}

# make a cluster id across both populations. brute force
ld[, cluster := NA_integer_]
lastclustcan <- NA # to detect when cluster id in can changes. no clusters were labeled 0
lastclustlof <- NA
clustid <- 0 # cluster id counter
for(i in 1:nrow(ld)){ # for each row
  if(i %% 1000 == 0) cat(i)  
  # if this locus is in a cluster (can or lof)
  if(!is.na(ld[i,cluster_can]) | !is.na(ld[i, cluster_lof])){ 
    
    # and this is a new cluster for both can and lof, increment cluster id and label it
    if((is.na(lastclustcan) | (!is.na(lastclustcan) & !is.na(ld[i, cluster_can]) & lastclustcan != ld[i, cluster_can])) & 
       (is.na(lastclustlof) | (!is.na(lastclustlof) & !is.na(ld[i, cluster_lof]) & lastclustlof != ld[i, cluster_lof]))){
      clustid <- clustid + 1
      ld[i, cluster := clustid]
    } else {
      # but use old label if not a new cluster in lof or can
      ld[i, cluster := clustid]
    }
  } 

  # update the cluster number for can and lof
  lastclustcan <- ld[i, cluster_can]
  lastclustlof <- ld[i, cluster_lof]
}

# find a locus near the middle of each block to keep
canclust <- ld[!is.na(cluster_can), .(CHROM = unique(CHROM), POS = findmid(POS), nloci = length(POS)), by = .(cluster = cluster_can)] # find the cluster midpoints
canclust <- rbind(canclust, ld[is.na(cluster_can), .(CHROM, POS, nloci = 1, cluster = NA_integer_)]) # add the loci not in blocks
setkey(canclust, CHROM, POS)
setcolorder(canclust, c('CHROM', 'POS', 'cluster', 'nloci'))

lofclust <- ld[!is.na(cluster_lof), .(CHROM = unique(CHROM), POS = findmid(POS), nloci = length(POS)), by = .(cluster = cluster_lof)] # find the cluster midpoints
lofclust <- rbind(lofclust, ld[is.na(cluster_lof), .(CHROM, POS, nloci = 1, cluster = NA_integer_)]) # add the loci not in blocks
setkey(lofclust, CHROM, POS)
setcolorder(lofclust, c('CHROM', 'POS', 'cluster', 'nloci'))

allclust <- ld[!is.na(cluster), .(CHROM = unique(CHROM), POS = findmid(POS), nloci = length(POS)), by = .(cluster = cluster)] # find the cluster midpoints
allclust <- rbind(allclust, ld[is.na(cluster), .(CHROM, POS, nloci = 1, cluster = NA_integer_)]) # add the loci not in blocks
setkey(allclust, CHROM, POS)
setcolorder(allclust, c('CHROM', 'POS', 'cluster', 'nloci'))


# compare
nrow(ld)
nrow(canclust)
nrow(lofclust)
nrow(allclust)

# write out loci trimmed to unlinked blocks. no quotes so can process easily on command line
write.csv(canclust, gzfile('analysis/ld.unlinked.Can.gatk.csv.gz'), quote = FALSE)
write.csv(lofclust, gzfile('analysis/ld.unlinked.Lof.gatk.csv.gz'), quote = FALSE)
write.csv(allclust, gzfile('analysis/ld.unlinked.gatk.csv.gz'), quote = FALSE)
