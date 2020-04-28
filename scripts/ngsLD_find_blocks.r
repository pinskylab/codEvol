# identify clusters of linked loci from ngsLD
# run after nsdLD_bypop.sh
# currently only set up for GATK loci

#parameters
minr2 <- 0.5 # consider pairs of loci with r2 > minr2 to be linked

# functions
require(data.table)

# read in list of loci
gatk <- fread('data_31_01_20/GATK_filtered_SNP_set.tab')


# read in LD data (can do this for gatk loci only at this time). files from ngsLD_bypop.sh
ldCan40 <- fread('analysis/ld.Can_40.gatk.gz')
ldCan14 <- fread('analysis/ld.Can_14.gatk.gz')
ld07 <- fread('analysis/ld.Lof_07.gatk.gz')
ld11 <- fread('analysis/ld.Lof_11.gatk.gz')
ld14 <- fread('analysis/ld.Lof_14.gatk.gz')

# add column names
nms <- c('pos1nm', 'pos2nm', 'dist', 'r2', 'D', 'Dprime', 'r2em')
setnames(ldCan14, nms)
setnames(ldCan40, nms)
setnames(ld14, nms)
setnames(ld11, nms)
setnames(ld07, nms)

# make a chromosome column
ldCan40[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
ldCan14[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
ld07[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
ld11[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
ld14[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]

# remove unplaced
ldCan40 <- ldCan40[!(chr == 'Unplaced'), ]
ldCan14 <- ldCan14[!(chr == 'Unplaced'), ]
ld07 <- ld07[!(chr == 'Unplaced'), ]
ld11 <- ld11[!(chr == 'Unplaced'), ]
ld14 <- ld14[!(chr == 'Unplaced'), ]

# make position columns
ldCan40[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
ldCan40[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
ldCan14[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
ldCan14[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
ld07[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
ld07[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
ld11[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
ld11[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
ld14[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
ld14[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
  

# cluster loci into linkage blocks
# average r2 within a population
ldCan <- merge(ldCan40[ ,.(chr, pos1, pos2, r21 = r2)], ldCan14[ ,.(chr, pos1, pos2, r22 = r2)], all = TRUE) 
ldCan[, r2 := rowMeans(cbind(r21, r22), na.rm = TRUE)]

ldLof <- merge(ld07[ ,.(chr, pos1, pos2, r21 = r2)], ld11[ ,.(chr, pos1, pos2, r22 = r2)], all = TRUE) 
ldLof <- merge(ldLof, ld14[ ,.(chr, pos1, pos2, r23 = r2)], all = TRUE) 
ldLof[, r2 := rowMeans(cbind(r21, r22, r23), na.rm = TRUE)]

rm(ldCan40, ldCan14, ld07, ld11, ld14)

# trim to locus pairs that are linked to cluster these
ldCan <- ldCan[r2 > minr2,]
ldLof <- ldLof[r2 > minr2,]

# order by r2
setorder(ldCan, chr, pos1, pos2)
setorder(ldLof, chr, pos1, pos2)

# find the most linked loci
ldcanpos1 <- ldCan[, .(numr21 = .N), by = .(chr, pos = pos1)] # number of linked loci for each locus (when listed in pos1)
ldcanpos2 <- ldCan[, .(numr22 = .N), by = .(chr, pos = pos2)] # for pos2
ldcanpos <- merge(ldcanpos1, ldcanpos2, all = TRUE) # merge
ldcanpos[, num := rowSums(cbind(numr21, numr22), na.rm = TRUE)] # total number of linked loci
setorder(ldcanpos, -num, chr, pos) # order decreasing by # linked loci

ldlofpos1 <- ldLof[, .(numr21 = .N), by = .(chr, pos = pos1)] # number of linked loci for each locus (when listed in pos1)
ldlofpos2 <- ldLof[, .(numr22 = .N), by = .(chr, pos = pos2)] # for pos2
ldlofpos <- merge(ldlofpos1, ldlofpos2, all = TRUE) # merge
ldlofpos[, num := rowSums(cbind(numr21, numr22), na.rm = TRUE)] # total number of linked loci
setorder(ldlofpos, -num, chr, pos) # decreasing by # linked loci

# label clusters in pairwise LD matrix
clustID <- 1 # for Can
ldCan[, cluster := NA_real_] # column for the cluster IDs in the pairwise LD dataset
nrow(ldcanpos)
for(i in 1:nrow(ldcanpos)){ # for each locus that is linked, working from most to least
  if(i %% 1000 == 0) cat(paste0(i, ' '))
  inds <- ldCan[, which(chr == ldcanpos$chr[i] & (pos1 == ldcanpos$pos[i] | pos2 == ldcanpos$pos[i]))] # get indices to pairwise ld that include this locus
  
  if(ldCan[inds, all(is.na(cluster))]){ # if no cluster IDs used yet for this locus
    ldCan[inds, cluster := clustID] # label with the next cluster ID
    clustID <- clustID + 1 # increment the cluster ID
  } else { # if a cluster ID (or more) has already been used
    thisclustIDs <- ldCan[inds, ][!is.na(cluster), unique(cluster)] # get the ID(s) already used
    if(length(thisclustIDs) == 1) ldCan[inds, cluster := thisclustIDs] # if only one, use it
    if(length(thisclustIDs) > 1){ # if the locus links two or more clusters
      thisclustID <- min(thisclustIDs) # pick the lowest cluster ID
      ldCan[cluster %in% thisclustIDs, cluster := thisclustID] # label all clusters with one ID 
      ldCan[inds, cluster := thisclustID] # label this locus
    }
  }
}

clustID <- 1 # for Lof
ldLof[, cluster := NA_real_] # column for the cluster IDs in the pairwise LD dataset
nrow(ldlofpos)
for(i in 1:nrow(ldlofpos)){ # for each locus that is linked, working from most to least
  if(i %% 1000 == 0) cat(paste0(i, ' '))
  inds <- ldLof[, which(chr == ldlofpos$chr[i] & (pos1 == ldlofpos$pos[i] | pos2 == ldlofpos$pos[i]))] # get indices to pairwise ld that include this locus
  
  if(ldLof[inds, all(is.na(cluster))]){ # if no cluster IDs used yet for this locus
    ldLof[inds, cluster := clustID] # label with the next cluster ID
    clustID <- clustID + 1 # increment the cluster ID
  } else { # if a cluster ID (or more) has already been used
    thisclustIDs <- ldLof[inds, ][!is.na(cluster), unique(cluster)] # get the ID(s) already used
    if(length(thisclustIDs) == 1) ldLof[inds, cluster := thisclustIDs] # if only one cluster ID, use it
    if(length(thisclustIDs) > 1){ # if the locus links two or more clusters
      thisclustID <- min(thisclustIDs) # pick the lowest cluster ID
      ldLof[cluster %in% thisclustIDs, cluster := thisclustID] # label all clusters with one ID 
      ldLof[inds, cluster := thisclustID] # label this locus
    }
  }
}

# label clusters in locus list
gatk[, cluster_can := NA_real_] # column for the CAN cluster IDs
clusts <- ldCan[, sort(unique(cluster), decreasing = TRUE)] # start from least linked locus cluster
for(i in 1:length(clusts)){
  inds <- gatk[, which(CHROM == ldCan[cluster == clusts[i], unique(chr)] & POS %in% ldCan[cluster == clusts[i], c(pos1, pos2)])]
  inds <- min(inds):max(inds) # make sure all loci from min to max position get labeled as part of this cluster
  gatk[inds, cluster_can := clusts[i]]
}

gatk[, cluster_lof := NA_real_] # column for the CAN cluster IDs
clusts <- ldLof[, sort(unique(cluster), decreasing = TRUE)] # start from least linked locus cluster
for(i in 1:length(clusts)){
  inds <- gatk[, which(CHROM == ldLof[cluster == clusts[i], unique(chr)] & POS %in% ldLof[cluster == clusts[i], c(pos1, pos2)])]
  inds <- min(inds):max(inds) # make sure all loci from min to max position get labeled as part of this cluster
  gatk[inds, cluster_lof := clusts[i]]
}


# write out
write.csv(gatk, gzfile('analysis/ld.blocks.gatk.csv.gz'))
