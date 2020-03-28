# shuffle ANGSD FST A and B values across sites and calculate windowed FST to get a null distribution of max genome-wide FST
# run last part of angsd_fst.sh first to get the *.fst.AB.gz files

# parameters
winsz <- 50000 # window size
winstp <- 10000 # window step

# load functions
require(data.table)

# load data
can <- fread('analysis/Can_40.Can_14.fst.AB.gz')
setnames(can, c('CHROM', 'POS', 'A', 'B'))

lof0711 <- fread('analysis/Lof_07.Lof_11.fst.AB.gz')
setnames(lof0711, c('CHROM', 'POS', 'A', 'B'))

lof0714 <- fread('analysis/Lof_07.Lof_14.fst.AB.gz')
setnames(lof0714, c('CHROM', 'POS', 'A', 'B'))

lof1114 <- fread('analysis/Lof_11.Lof_14.fst.AB.gz')
setnames(lof1114, c('CHROM', 'POS', 'A', 'B'))

# remove unplaced
can <- can[!(CHROM %in% 'Unplaced'), ]
lof0711 <- lof0711[!(CHROM %in% 'Unplaced'), ]
lof0714 <- lof0714[!(CHROM %in% 'Unplaced'), ]
lof1114 <- lof1114[!(CHROM %in% 'Unplaced'), ]

# shuffle and recalc windowed FST
# CAN
for(i in 1:1000){
	cat(i)
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
	if(i == 1) maxfst <- tempfsts[, max(fst)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[, max(fst)])
}

write.csv(maxfst, gzfile('analysis/Can_40.Can_14.fst.siteshuffle.csv.gz'), row.names = FALSE)

rm(maxfst)

# Lof0711
for(i in 1:1000){
	cat(i)
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
	if(i == 1) maxfst <- tempfsts[, max(fst)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[, max(fst)])
}

write.csv(maxfst, gzfile('analysis/Lof_07.Lof_11.fst.siteshuffle.csv.gz'), row.names = FALSE)

rm(maxfst)

# Lof0714
for(i in 1:1000){
	cat(i)
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
	if(i == 1) maxfst <- tempfsts[, max(fst)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[, max(fst)])
}

write.csv(maxfst, gzfile('analysis/Lof_07.Lof_14.fst.siteshuffle.csv.gz'), row.names = FALSE)

rm(maxfst)

# Lof1114
for(i in 1:1000){
	cat(i)
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
	if(i == 1) maxfst <- tempfsts[, max(fst)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[, max(fst)])
}

write.csv(maxfst, gzfile('analysis/Lof_11.Lof_14.fst.siteshuffle.csv.gz'), row.names = FALSE)