# shuffle theta values across sites and calculate windowed averages to get a null distribution of max genome-wide thetas

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) != 1) stop("Have to specify whether all loci (0) or GATK loci (1)", call.=FALSE)

gatkflag <- as.numeric(args[1])

if(gatkflag == 1){
  stop('Using GATK loci, but code not set up for that yet')
} else if(gatkflag == 0){
  print('Using all loci')
} else {
	stop(paste(gatkflag, ' is not 0 or 1. Please only specify 0 (all loci) or 1 (gatk loci).'))
}


# parameters
winsz <- 50000 # window size
winstp <- 10000 # window step
nrep <- 1000 # number of reshuffles

nchrCan40 <- 42 # sample size in # chromosomes
nchrCan14 <- 48
nchrLof07 <- 44
nchrLof11 <- 48
nchrLof14 <- 44

outfilecan <- 'analysis/theta.siteshuffle.Can_40.Can_14.csv.gz' # used if all loci are used
outfilelof0711 <- 'analysis/theta.siteshuffle.Lof_07.Lof_11.csv.gz'
outfilelof0714 <- 'analysis/theta.siteshuffle.Lof_07.Lof_14.csv.gz'
outfilelof1114 <- 'analysis/theta.siteshuffle.Lof_11.Lof_14.csv.gz'

####################
# load functions
require(data.table)

# after https://github.com/ANGSD/angsd/blob/master/misc/stats.cpp
a1f <- function(nsam) return(sum(1/seq(1, nsam-1)))
a2f <- function(nsam) return(sum(1/(seq(1, nsam-1)*seq(1, nsam-1))))
b1f <- function(nsam) return((nsam + 1)/(3*(nsam-1)))
b2f <- function(nsam) return((2*(nsam*nsam + nsam + 3))/(9*nsam*(nsam - 1)))
c1f <- function(a1, b1) return(b1 - (1/a1))
c2f <- function(nsam, a1, a2, b2) return(b2 - ((nsam + 2)/(a1*nsam)) + (a2/(a1 * a1)))
e1f <- function(a1, c1) return(c1/a1)
e2f <- function(a1, a2, c2) return(c2/((a1*a1) + a2))

# Tajima's D calculation
# nsam: sample size
# thetaW: Watterson's theta (# segregating sites / a1)
# sumk: theta pi (average number of SNPs in pairwise comparisons)
# after https://github.com/ANGSD/angsd/blob/master/misc/stats.cpp
tajd <- function(nsam, thetaW, sumk){
	a1 <- a1f(nsam)
	segsites <- thetaW * a1
	if(segsites == 0) return(0)
	a2 <- a2f(nsam)
	b1 <- b1f(nsam)
  	b2 <- b2f(nsam)
	c1 <- c1f(a1, b1)
	c2 <- c2f(nsam, a1, a2, b2)
	e1 <- e1f(a1, c1)
	e2 <- e2f(a1, a2, c2)
	res <- (sumk - (thetaW))/sqrt((e1*segsites) + ((e2*segsites)*(segsites-1)))
	return(res)
}

# calculate stats from reshuffled data
shufflestat <- function(dat1, dat2, nchr1, nchr2){

	inds1 <- sample(1:nrow(dat1), nrow(dat1), replace = FALSE)
	temp1 <- cbind(dat1[, .(Chromo, Pos)], dat1[inds1, .(Watterson, Pairwise)]) # shuffle theta across position pairs
	inds2 <- sample(1:nrow(dat2), nrow(dat2), replace = FALSE)
	temp2 <- cbind(dat2[, .(Chromo, Pos)], dat2[inds2, .(Watterson, Pairwise)]) # second time-point

	# create new columns as indices for windows
	# need multiple columns because overlapping windows
	for(j in 1:(winsz/winstp)){
		temp1[, (paste0('win_', j)) := floor((Pos - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
		temp2[, (paste0('win_', j)) := floor((Pos - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	}
	
	# calc ave theta for each window
	for(j in 1:(winsz/winstp)){
		nm1 <- paste0('win_', j)
		if(j ==1){
		  bins1 <- temp1[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW1 = sum(exp(Watterson), na.rm = TRUE),
		                                                                       tP1 = sum(exp(Pairwise), na.rm = TRUE)), 
		                                                                   by = .(Chromo, WinCenter)]
		  bins2 <- temp2[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW2 = sum(exp(Watterson), na.rm = TRUE),
		                                                                             tP2 = sum(exp(Pairwise), na.rm = TRUE)), 
		                                                                         by = .(Chromo, WinCenter)]
		}
		if(j > 1){
		  bins1 <- rbind(bins1, temp1[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW1 = sum(exp(Watterson), na.rm = TRUE),
		                                                                             tP1 = sum(exp(Pairwise), na.rm = TRUE)), 
		                                                                         by = .(Chromo, WinCenter)])
		  bins2 <- rbind(bins2, temp2[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW2 = sum(exp(Watterson), na.rm = TRUE),
		                                                                             tP2 = sum(exp(Pairwise), na.rm = TRUE)), 
		                                                                         by = .(Chromo, WinCenter)])
		}
	}

	# calculate Tajima's D
	bins1[, tD1 := tajd(nchr1, tW1, tP1), by = 1:nrow(bins1)]
	bins2[, tD2 := tajd(nchr2, tW2, tP2), by = 1:nrow(bins2)]

	# merge by window and calculate change
	bins <- merge(bins1, bins2, by = c("Chromo", "WinCenter"))
	bins[, ':='(tWd = tW2 - tW1, tPd = tP2 - tP1, tDd = tD2 - tD1)]

	# return the max and min windowed change
	# exclude windows with negative midpoints
	return(bins[WinCenter >= winsz/2, .(mintWd = min(tWd, na.rm = TRUE), 
								maxtWd = max(tWd, na.rm = TRUE),
								mintPd = min(tPd, na.rm = TRUE), 
								maxtPd = max(tPd, na.rm = TRUE),
								mintDd = min(tDd, na.rm = TRUE), 
								maxtDd = max(tDd, na.rm = TRUE))])
}



# load all loci theta calcs if requested
if(gatkflag == 0){
	# read in data
	datCan40 <- fread('analysis/thetas.Can_40.pestPG.gz')
	datCan14 <- fread('analysis/thetas.Can_14.pestPG.gz')
	dat07 <- fread('analysis/thetas.Lof_07.pestPG.gz')
	dat11 <- fread('analysis/thetas.Lof_11.pestPG.gz')
	dat14 <- fread('analysis/thetas.Lof_14.pestPG.gz')
}

if(gatkflag == 1){
  # read in files when available

  # set new outfile names
	outfilecan <- 'analysis/theta.siteshuffle.Can_40.Can_14.gatk.csv.gz' # used if all loci are used
	outfilelof0711 <- 'analysis/theta.siteshuffle.Lof_07.Lof_11.gatk.csv.gz'
	outfilelof0714 <- 'analysis/theta.siteshuffle.Lof_07.Lof_14.gatk.csv.gz'
	outfilelof1114 <- 'analysis/theta.siteshuffle.Lof_11.Lof_14.gatk.csv.gz'
}

# fix name
setnames(datCan40, '#Chromo', 'Chromo')
setnames(datCan14, '#Chromo', 'Chromo')
setnames(dat07, '#Chromo', 'Chromo')
setnames(dat11, '#Chromo', 'Chromo')
setnames(dat14, '#Chromo', 'Chromo')

# remove unplaced
datCan40 <- datCan40[grep('Unplaced', Chromo, invert = TRUE), ]
datCan14 <- datCan14[grep('Unplaced', Chromo, invert = TRUE), ]
dat07 <- dat07[grep('Unplaced', Chromo, invert = TRUE), ]
dat11 <- dat11[grep('Unplaced', Chromo, invert = TRUE), ]
dat14 <- dat14[grep('Unplaced', Chromo, invert = TRUE), ]



# shuffle and recalc windowed LD
# CAN
print('Starting Can')
for(i in 1:nrep){
	cat(i); cat(' ')
	tempminmax <- shufflestat(datCan40, datCan14, nchrCan40, nchrCan14)

	if(i == 1) minmax <- tempminmax
	if(i > 1) minmax <- rbind(minmax, tempminmax)
}

write.csv(minmax, gzfile(outfilecan), row.names = FALSE)
rm(minmax)



# Lof0711
print('Starting Lof0711')
for(i in 1:nrep){
	cat(i); cat(' ')
	tempminmax <- shufflestat(datLof07, datLof11, nchrLof07, nchrLof11)

	if(i == 1) minmax <- tempminmax
	if(i > 1) minmax <- rbind(minmax, tempminmax)
}

write.csv(minmax, gzfile(outfilelof0711), row.names = FALSE)
rm(minmax)


# Lof0714
print('Starting Lof0714')
for(i in 1:nrep){
	cat(i); cat(' ')
	tempminmax <- shufflestat(datLof07, datLof14, nchrLof07, nchrLof14)

	if(i == 1) minmax <- tempminmax
	if(i > 1) minmax <- rbind(minmax, tempminmax)
}

write.csv(minmax, gzfile(outfilelof0714), row.names = FALSE)
rm(minmax)



# Lof1114
print('Starting Lof1114')
for(i in 1:nrep){
	cat(i); cat(' ')
	tempminmax <- shufflestat(datLof07, datLof14, nchrLof07, nchrLof14)

	if(i == 1) minmax <- tempminmax
	if(i > 1) minmax <- rbind(minmax, tempminmax)
}

write.csv(minmax, gzfile(outfilelof1114), row.names = FALSE)
rm(minmax)