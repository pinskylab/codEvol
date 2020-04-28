# shuffle ANGSD FST A and B values across sites and calculate windowed FST to get a null distribution of max genome-wide FST
# run last part of angsd_fst.sh first to get the *.fst.AB.gz files
# if using GATK loci, groups into linkage blocks based on ngsLD output. Need to run ngsLD_find_blocks.sh/.r first

# to run on saga

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) != 1) stop("Have to specify whether all loci (0) or GATK loci (1)", call.=FALSE)

gatkflag <- as.numeric(args[1])

if(gatkflag == 1){
	print('Using only GATK loci. Will trim to unlinked blocks.')
} else if(gatkflag == 0){
	print('Using all loci. Note that this option does not trim to unlinked groups yet.')
} else {
	stop(paste(gatkflag, ' is not 0 or 1. Please only specify 0 (all loci) or 1 (gatk loci).'))
}


# parameters
winsz <- 50000 # window size
winstp <- 10000 # window step
nrep <- 1000 # number of reshuffles

outfilecan <- 'analysis/Can_40.Can_14.fst.siteshuffle.csv.gz' # used if all loci are used
outfilelof0711 <- 'analysis/Lof_07.Lof_11.fst.siteshuffle.csv.gz'
outfilelof0714 <- 'analysis/Lof_07.Lof_14.fst.siteshuffle.csv.gz'
outfilelof1114 <- 'analysis/Lof_11.Lof_14.fst.siteshuffle.csv.gz'

# load functions
require(data.table)

# load fst A/B data
can <- fread('analysis/Can_40.Can_14.fst.AB.gz')
setnames(can, c('CHROM', 'POS', 'A', 'B'))

lof0711 <- fread('analysis/Lof_07.Lof_11.fst.AB.gz')
setnames(lof0711, c('CHROM', 'POS', 'A', 'B'))

lof0714 <- fread('analysis/Lof_07.Lof_14.fst.AB.gz')
setnames(lof0714, c('CHROM', 'POS', 'A', 'B'))

lof1114 <- fread('analysis/Lof_11.Lof_14.fst.AB.gz')
setnames(lof1114, c('CHROM', 'POS', 'A', 'B'))

# trim to gatk loci if requested
if(gatkflag == 1){
	# list of loci to use
	gatk <- fread('data_31_01_20/GATK_filtered_SNP_set.tab')

	# trim
	can <- can[gatk, on = c('CHROM', 'POS')]
	lof0711 <- lof0711[gatk, on = c('CHROM', 'POS')]
	lof0714 <- lof0714[gatk, on = c('CHROM', 'POS')]
	lof1114 <- lof1114[gatk, on = c('CHROM', 'POS')]

	# set new outfile names
	outfilecan <- 'analysis/Can_40.Can_14.gatk.fst.siteshuffle.csv.gz'
	outfilelof0711 <- 'analysis/Lof_07.Lof_11.gatk.fst.siteshuffle.csv.gz'
	outfilelof0714 <- 'analysis/Lof_07.Lof_14.gatk.fst.siteshuffle.csv.gz'
	outfilelof1114 <- 'analysis/Lof_11.Lof_14.gatk.fst.siteshuffle.csv.gz'
}

# remove unplaced
can <- can[!(CHROM %in% 'Unplaced'), ]
lof0711 <- lof0711[!(CHROM %in% 'Unplaced'), ]
lof0714 <- lof0714[!(CHROM %in% 'Unplaced'), ]
lof1114 <- lof1114[!(CHROM %in% 'Unplaced'), ]


# trim to unlinked loci if GATK loci were requested
if(gatkflag == 1){
	ld <- fread('analysis/ld.blocks.gatk.csv.gz') # linkage blocks from ngsLD_find_blocks.r
	
	can <- merge(can, ld[, .(CHROM, POS, cluster = cluster_can)], all.x = TRUE)
	lof0711 <- merge(lof0711, ld[, .(CHROM, POS, cluster = cluster_lof)], all.x = TRUE)
	lof0714 <- merge(lof0714, ld[, .(CHROM, POS, cluster = cluster_lof)], all.x = TRUE)
	lof1114 <- merge(lof1114, ld[, .(CHROM, POS, cluster = cluster_lof)], all.x = TRUE)

	# average Fst in linkage blocks and find a locus near the middle of each block to keep
	findmid <- function(POS){ # function to return a value near the middle of a vector of positions
		mn <- mean(range(POS))
		return(POS[which.min(abs(POS - mn))])
	}

	can[, keep := 1] # column for marking which loci to keep
	can[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
	canclust <- can[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
	can <- merge(can, canclust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
	can[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
	can <- can[keep == 1, ] # drop all loci in clusters that aren't midpoints
	can[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns

	lof0711[, keep := 1] # column for marking which loci to keep
	lof0711[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
	lof0711clust <- lof0711[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
	lof0711 <- merge(lof0711, lof0711clust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
	lof0711[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
	lof0711 <- lof0711[keep == 1, ] # drop all loci in clusters that aren't midpoints
	lof0711[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns

	lof0714[, keep := 1] # column for marking which loci to keep
	lof0714[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
	lof0714clust <- lof0714[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
	lof0714 <- merge(lof0714, lof0714clust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
	lof0714[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
	lof0714 <- lof0714[keep == 1, ] # drop all loci in clusters that aren't midpoints
	lof0714[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns

	lof1114[, keep := 1] # column for marking which loci to keep
	lof1114[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
	lof1114clust <- lof1114[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
	lof1114 <- merge(lof1114, lof1114clust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
	lof1114[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
	lof1114 <- lof1114[keep == 1, ] # drop all loci in clusters that aren't midpoints
	lof1114[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns

	# write out fst trimmed to unlinked blocks
	write.csv(can, gzfile('analysis/Can_40.Can_14.fst.AB.ldtrim.gz'))
	write.csv(lof0711, gzfile('analysis/Lof_07.Lof_11.fst.AB.ldtrim.gz'))
	write.csv(lof0714, gzfile('analysis/Lof_07.Lof_14.fst.AB.ldtrim.gz'))
	write.csv(lof1114, gzfile('analysis/Lof_11.Lof_14.fst.AB.ldtrim.gz'))

}




# shuffle and recalc windowed FST
# CAN
print('Starting Can')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(can), nrow(can), replace = FALSE)
	temp <- cbind(can[, .(CHROM, POS)], can[inds, .(A, B)]) # shuffle FSTs across positions
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc fst for each window
	for(j in 1:(winsz/winstp - 1)){
		if(j ==1) tempfsts <- temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	# exclude windows with negative midpoints
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfilecan), row.names = FALSE)

rm(maxfst)



# Lof0711
print('Starting Lof0711')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(lof0711), nrow(lof0711), replace = FALSE)
	temp <- cbind(lof0711[, .(CHROM, POS)], lof0711[inds, .(A, B)]) # shuffle FSTs across positions
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc fst for each window
	for(j in 1:(winsz/winstp - 1)){
		if(j ==1) tempfsts <- temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfilelof0711), row.names = FALSE)

rm(maxfst)


# Lof0714
print('Starting Lof0714')
for(i in 1:nrep){
 	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(lof0714), nrow(lof0714), replace = FALSE)
	temp <- cbind(lof0714[, .(CHROM, POS)], lof0714[inds, .(A, B)]) # shuffle FSTs across positions
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc fst for each window
	for(j in 1:(winsz/winstp - 1)){
		if(j ==1) tempfsts <- temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfilelof0714), row.names = FALSE)

rm(maxfst)


# Lof1114
print('Starting Lof1114')
for(i in 1:nrep){
 	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(lof1114), nrow(lof1114), replace = FALSE)
	temp <- cbind(lof1114[, .(CHROM, POS)], lof1114[inds, .(A, B)]) # shuffle FSTs across positions
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc fst for each window
	for(j in 1:(winsz/winstp - 1)){
		if(j ==1) tempfsts <- temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfilelof1114), row.names = FALSE)