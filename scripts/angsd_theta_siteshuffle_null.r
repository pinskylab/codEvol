# shuffle theta values across sites and calculate windowed averages to get a null distribution of max genome-wide thetas

#################################
# read command line arguments
#################################
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) != 2) stop("Have to specify whether all loci (0) or GATK nodam2 unlinked loci (1), and which pop (1=Can, 2=Lof0711, 3=Lof0714, 4=Lof1114)", call.=FALSE)

gatkflag <- as.numeric(args[1])
popflag <- as.numeric(args[2])

if(gatkflag == 1){
  print('Using GATK nodam2 unlinked loci')
} else if(gatkflag == 0){
  print('Using all loci')
} else {
	stop(paste(gatkflag, ' is not 0 or 1. Please only specify 0 (all loci) or 1 (gatk loci).'))
}
if(!(popflag %in% 1:4)) stop('Have to specify pop in 1, 2, 3 or 4 (1=Can, 2=Lof0711, 3=Lof0714, 4=Lof1114)')

#################################
# parameters
#################################
winsz <- 50000 # window size
winstp <- 10000 # window step
nrep <- 1000 # number of reshuffles

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

# calculate thetas in windows from data
calcstat <- function(temp1, temp2, nchr1, nchr2){
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
  
  return(bins)
}

# calculate min and max change from reshuffled data
shufflestat <- function(dat1, dat2, nchr1, nchr2){
	inds1 <- sample(1:nrow(dat1), nrow(dat1), replace = FALSE)
	temp1 <- cbind(dat1[, .(Chromo, Pos)], dat1[inds1, .(Watterson, Pairwise)]) # shuffle theta across position pairs
	inds2 <- sample(1:nrow(dat2), nrow(dat2), replace = FALSE)
	temp2 <- cbind(dat2[, .(Chromo, Pos)], dat2[inds2, .(Watterson, Pairwise)]) # second time-point

	bins <- calcstat(temp1, temp2, nchr1, nchr2)

	# return the max and min windowed change
	# exclude windows with negative midpoints
	return(bins[WinCenter >= winsz/2, .(mintWd = min(tWd, na.rm = TRUE), 
	                                    maxtWd = max(tWd, na.rm = TRUE),
	                                    mintPd = min(tPd, na.rm = TRUE), 
	                                    maxtPd = max(tPd, na.rm = TRUE),
	                                    mintDd = min(tDd, na.rm = TRUE), 
	                                    maxtDd = max(tDd, na.rm = TRUE))])
	
}


######################
# Load data
######################

# load all loci theta calcs if requested
if(gatkflag == 0){
	if(popflag == 1){
		print('Starting Can')
		dat1 <- fread('analysis/thetas.Can_40.pestPG.gz')
		dat2 <- fread('analysis/thetas.Can_14.pestPG.gz')
		nchr1 <- 42 # sample size in # chromosomes
		nchr2 <- 48
		outfilesuff <- 'Can_40.Can_14'
	}
	if(popflag == 2){
		print('Starting Lof0711')
		dat1 <- fread('analysis/thetas.Lof_07.pestPG.gz')
		dat2 <- fread('analysis/thetas.Lof_11.pestPG.gz')
		nchr1 <- 44
		nchr2 <- 48
		outfilesuff <- 'Lof_07.Lof_11'
	}
	if(popflag == 3){
		print('Starting Lof0714')
		dat1 <- fread('analysis/thetas.Lof_07.pestPG.gz')
		dat2 <- fread('analysis/thetas.Lof_14.pestPG.gz')
		nchr1 <- 44
		nchr2 <- 44
		outfilesuff <- 'Lof_07.Lof_14'
	}
	if(popflag == 4){
		print('Starting Lof1114')
		dat1 <- fread('analysis/thetas.Lof_11.pestPG.gz')
		dat2 <- fread('analysis/thetas.Lof_14.pestPG.gz')
		nchr1 <- 48
		nchr2 <- 44
		outfilesuff <- 'Lof_11.Lof_14'
	}
}

if(gatkflag == 1){
	if(popflag == 1){
		print('Starting Can')
		dat1 <- fread('analysis/thetas.Can_40.gatk.pestPG.gz')
		dat2 <- fread('analysis/thetas.Can_14.gatk.pestPG.gz')
		nchr1 <- 42 # sample size in # chromosomes
		nchr2 <- 48
		outfilesuff <- 'Can_40.Can_14.gatk'
	}
	if(popflag == 2){
		print('Starting Lof0711')
		dat1 <- fread('analysis/thetas.Lof_07.gatk.pestPG.gz')
		dat2 <- fread('analysis/thetas.Lof_11.gatk.pestPG.gz')
		nchr1 <- 44
		nchr2 <- 48
		outfilesuff <- 'Lof_07.Lof_11.gatk'
	}
	if(popflag == 3){
		print('Starting Lof0714')
		dat1 <- fread('analysis/thetas.Lof_07.gatk.pestPG.gz')
		dat2 <- fread('analysis/thetas.Lof_14.gatk.pestPG.gz')
		nchr1 <- 44
		nchr2 <- 44
		outfilesuff <- 'Lof_07.Lof_14.gatk'
	}
	if(popflag == 4){
		print('Starting Lof1114')
		dat1 <- fread('analysis/thetas.Lof_11.gatk.pestPG.gz')
		dat2 <- fread('analysis/thetas.Lof_14.gatk.pestPG.gz')
		nchr1 <- 48
		nchr2 <- 44
		outfilesuff <- 'Lof_11.Lof_14.gatk'
	}
}

# fix name
setnames(dat1, '#Chromo', 'Chromo')
setnames(dat2, '#Chromo', 'Chromo')

# remove unplaced
dat1 <- dat1[grep('Unplaced', Chromo, invert = TRUE), ]
dat2 <- dat2[grep('Unplaced', Chromo, invert = TRUE), ]

# trim to nodam2 loci if gatk
if(gatkflag == 1){
  print('Trimming to nodam2 loci')
  nodam2 <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab') # list of loci that pass nodam2 filter
	setnames(nodam2, c('Chromo', 'Pos', 'REF', 'ALT'))
	dat1 <- merge(dat1, nodam2[, .(Chromo, Pos)])
	dat2 <- merge(dat2, nodam2[, .(Chromo, Pos)])
}

# mark unlinked loci if gatk
if(gatkflag == 1){
  print('Trimming to unlinked loci')
	if(popflag == 1) unl <- fread('analysis/ld.unlinked.Can.gatk.nodam.csv.gz') # list of unlinked loci for Canada
	if(popflag %in% 2:4) unl <- fread('analysis/ld.unlinked.Lof.gatk.nodam.csv.gz') # for Lofoten
	dat1 <- merge(dat1, unl[, .(Chromo = CHROM, Pos = POS, keep = 1)], all.x = TRUE)
	dat2 <- merge(dat2, unl[, .(Chromo = CHROM, Pos = POS, keep = 1)], all.x = TRUE)
	dat1[is.na(keep), keep := 0]
	dat2[is.na(keep), keep := 0]
}


#######################################
# Calculate change on unshuffled data
#######################################

# full dataset
out <- calcstat(dat1, dat2, nchr1, nchr2)
outfile1 <- paste0('analysis/theta_change_region_', winsz, '.', outfilesuff, '.csv.gz')
write.csv(out, file = gzfile(outfile1))
print(paste0('Wrote ', outfile1))

# unlinked (if GATK)
if(gatkflag == 1){
  out2 <- calcstat(dat1[keep == 1, ], dat2[keep == 1, ], nchr1, nchr2)
  outfile2 <- paste0('analysis/theta_change_region_', winsz, '.', outfilesuff, '.unlinked.csv.gz')
  write.csv(out2, file = gzfile(outfile2))
  print(paste0('Wrote ', outfile2))
}


#################################
# Trim to unlinked sites if GATK
#################################

if(gatkflag == 1){
  dat1 <- dat1[keep == 1, ]
  dat2 <- dat2[keep == 1, ]
}

################################
# Run reshuffle calculations
# shuffle and recalc windowed thetas
################################

for(i in 1:nrep){
	cat(i); cat(' ')
	tempminmax <- shufflestat(dat1, dat2, nchr1, nchr2)

	if(i == 1) minmax <- tempminmax # save the results from this iteration
	if(i > 1) minmax <- rbind(minmax, tempminmax)
}

write.csv(minmax, gzfile(paste0('analysis/theta.siteshuffle.', outfilesuff, '.csv.gz')), row.names = FALSE) # write out all iterations

rm(minmax)
