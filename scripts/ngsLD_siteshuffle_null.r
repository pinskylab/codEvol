# shuffle ngsLD values across sites and calculate windowed LD to get a null distribution of max genome-wide LD
# run ngsLD_bypop.sh first to get input files

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) != 1) stop("Have to specify whether all loci (0) or GATK loci (1)", call.=FALSE)

gatkflag <- as.numeric(args[1])

if(gatkflag == 1){
	print('Using only GATK loci')
} else if(gatkflag == 0){
	stop('Using all loci, but code not set up for that yet')
} else {
	stop(paste(gatkflag, ' is not 0 or 1. Please only specify 0 (all loci) or 1 (gatk loci).'))
}


# parameters
winsz <- 50000 # window size
winstp <- 10000 # window step
nrep <- 1000 # number of reshuffles

outfilecan <- 'analysis/ld.siteshuffle.Can_40.Can_14.csv.gz' # used if all loci are used
outfilelof0711 <- 'analysis/ld.siteshuffle.Lof_07.Lof_11.csv.gz'
outfilelof0714 <- 'analysis/ld.siteshuffle.Lof_07.Lof_14.csv.gz'
outfilelof1114 <- 'analysis/ld.siteshuffle.Lof_11.Lof_14.csv.gz'

# load functions
require(data.table)


# load gatk loci LD calcs if requested
if(gatkflag == 1){
	# read in data
	datCan40 <- fread('analysis/ld.Can_40.gatk.gz')
	datCan14 <- fread('analysis/ld.Can_14.gatk.gz')
	dat07 <- fread('analysis/ld.Lof_07.gatk.gz')
	dat11 <- fread('analysis/ld.Lof_11.gatk.gz')
	dat14 <- fread('analysis/ld.Lof_14.gatk.gz')

	# add column names
	nms <- c('pos1nm', 'pos2nm', 'dist', 'r2', 'D', 'Dprime', 'r2em')
	setnames(datCan40, nms)
	setnames(datCan14, nms)
	setnames(dat14, nms)
	setnames(dat11, nms)
	setnames(dat07, nms)

	# set new outfile names
	outfilecan <- 'analysis/ld.siteshuffle.Can_40.Can_14.gatk.csv.gz' # used if all loci are used
	outfilelof0711 <- 'analysis/ld.siteshuffle.Lof_07.Lof_11.gatk.csv.gz'
	outfilelof0714 <- 'analysis/ld.siteshuffle.Lof_07.Lof_14.gatk.csv.gz'
	outfilelof1114 <- 'analysis/ld.siteshuffle.Lof_11.Lof_14.gatk.csv.gz'
}

# remove unplaced
datCan40 <- datCan40[grep('Unplaced', pos1nm, invert = TRUE), ]
datCan14 <- datCan14[grep('Unplaced', pos1nm, invert = TRUE), ]
dat07 <- dat07[grep('Unplaced', pos1nm, invert = TRUE), ]
dat11 <- dat11[grep('Unplaced', pos1nm, invert = TRUE), ]
dat14 <- dat14[grep('Unplaced', pos1nm, invert = TRUE), ]

# merge LD values by locus pair
can <- merge(datCan40[, .(pos1nm, pos2nm, r2_1 = r2)], datCan14[, .(pos1nm, pos2nm, r2_2 = r2)])
lof0711 <- merge(dat07[, .(pos1nm, pos2nm, r2_1 = r2)], dat11[, .(pos1nm, pos2nm, r2_2 = r2)])
lof0714 <- merge(dat07[, .(pos1nm, pos2nm, r2_1 = r2)], dat14[, .(pos1nm, pos2nm, r2_2 = r2)])
lof1114 <- merge(dat11[, .(pos1nm, pos2nm, r2_1 = r2)], dat14[, .(pos1nm, pos2nm, r2_2 = r2)])

# calculate change by locus pair
can[, r2d := r2_2 - r2_1]
lof0711[, r2d := r2_2 - r2_1]
lof0714[, r2d := r2_2 - r2_1]
lof1114[, r2d := r2_2 - r2_1]

# make a chromosome column
can[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
lof0711[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
lof0714[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]
lof1114[, chr := vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 1)]

# make position columns
can[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
can[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
lof0711[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
lof0711[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
lof0714[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
lof0714[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]
lof1114[, pos1 := as.numeric(vapply(strsplit(pos1nm, ":", fixed = TRUE), "[", "", 2))]
lof1114[, pos2 := as.numeric(vapply(strsplit(pos2nm, ":", fixed = TRUE), "[", "", 2))]


# shuffle and recalc windowed LD
# CAN
print('Starting Can')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(can), nrow(can), replace = FALSE)
	temp <- cbind(can[, .(chr, pos1, pos2)], can[inds, .(r2d)]) # shuffle LDs across position pairs
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win1_', j)) := floor((pos1 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
		temp[, (paste0('win2_', j)) := floor((pos2 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc ave LD change for each window
	for(j in 1:(winsz/winstp - 1)){
		nm1 <- paste0('win1_', j)
		nm2 <- paste0('win2_', j)
		if(j ==1) templds <- temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)]
		if(j > 1) templds <- rbind(templds, temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)])
	}

	# save the max windowed ld change
	# exclude windows with negative midpoints
	if(i == 1) maxld <- templds[POS > 0, max(r2dave, na.rm = TRUE)]	
	if(i > 1) maxld <- c(maxld, templds[POS > 0, max(r2dave, na.rm = TRUE)])
}

print(paste('Max:', max(maxld, na.rm = TRUE), '; 95th:', quantile(maxld, prob = 0.95, na.rm = TRUE)))

write.csv(maxld, gzfile(outfilecan), row.names = FALSE)

rm(maxld)



# Lof0711
print('Starting Lof0711')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(lof0711), nrow(lof0711), replace = FALSE)
	temp <- cbind(lof0711[, .(chr, pos1, pos2)], lof0711[inds, .(r2d)]) # shuffle LDs across position pairs
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win1_', j)) := floor((pos1 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
		temp[, (paste0('win2_', j)) := floor((pos2 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc ave LD change for each window
	for(j in 1:(winsz/winstp - 1)){
		nm1 <- paste0('win1_', j)
		nm2 <- paste0('win2_', j)
		if(j ==1) templds <- temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)]
		if(j > 1) templds <- rbind(templds, temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)])
	}

	# save the max windowed ld change
	# exclude windows with negative midpoints
	if(i == 1) maxld <- templds[POS > 0, max(r2dave, na.rm = TRUE)]	
	if(i > 1) maxld <- c(maxld, templds[POS > 0, max(r2dave, na.rm = TRUE)])
}

print(paste('Max:', max(maxld, na.rm = TRUE), '; 95th:', quantile(maxld, prob = 0.95, na.rm = TRUE)))

write.csv(maxld, gzfile(outfilelof0711), row.names = FALSE)

rm(maxld)


# Lof0714
print('Starting Lof0714')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(lof0714), nrow(lof0714), replace = FALSE)
	temp <- cbind(lof0714[, .(chr, pos1, pos2)], lof0714[inds, .(r2d)]) # shuffle LDs across position pairs
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win1_', j)) := floor((pos1 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
		temp[, (paste0('win2_', j)) := floor((pos2 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc ave LD change for each window
	for(j in 1:(winsz/winstp - 1)){
		nm1 <- paste0('win1_', j)
		nm2 <- paste0('win2_', j)
		if(j ==1) templds <- temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)]
		if(j > 1) templds <- rbind(templds, temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)])
	}

	# save the max windowed ld change
	# exclude windows with negative midpoints
	if(i == 1) maxld <- templds[POS > 0, max(r2dave, na.rm = TRUE)]	
	if(i > 1) maxld <- c(maxld, templds[POS > 0, max(r2dave, na.rm = TRUE)])
}

print(paste('Max:', max(maxld, na.rm = TRUE), '; 95th:', quantile(maxld, prob = 0.95, na.rm = TRUE)))

write.csv(maxld, gzfile(outfilelof0714), row.names = FALSE)

rm(maxld)



# Lof1114
print('Starting Lof1114')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(lof1114), nrow(lof1114), replace = FALSE)
	temp <- cbind(lof1114[, .(chr, pos1, pos2)], lof1114[inds, .(r2d)]) # shuffle LDs across position pairs
	
	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp - 1)){
		temp[, (paste0('win1_', j)) := floor((pos1 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
		temp[, (paste0('win2_', j)) := floor((pos2 - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc ave LD change for each window
	for(j in 1:(winsz/winstp - 1)){
		nm1 <- paste0('win1_', j)
		nm2 <- paste0('win2_', j)
		if(j ==1) templds <- temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)]
		if(j > 1) templds <- rbind(templds, temp[!is.na(r2d) & get(nm1) == get(nm2), .(chr, r2d, POS = get(nm1))][, .(r2dave = mean(r2d)), by = .(chr, POS)])
	}

	# save the max windowed ld change
	# exclude windows with negative midpoints
	if(i == 1) maxld <- templds[POS > 0, max(r2dave, na.rm = TRUE)]	
	if(i > 1) maxld <- c(maxld, templds[POS > 0, max(r2dave, na.rm = TRUE)])
}

print(paste('Max:', max(maxld, na.rm = TRUE), '; 95th:', quantile(maxld, prob = 0.95, na.rm = TRUE)))

write.csv(maxld, gzfile(outfilelof1114), row.names = FALSE)

rm(maxld)