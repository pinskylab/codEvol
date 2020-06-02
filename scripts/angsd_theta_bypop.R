# average theta values across sites

#############
# Parameters
#############
# number of chromosomes in each sample
nchrCan40 <- 42 # sample size in # chromosomes
nchrCan14 <- 48
nchrLof07 <- 44
nchrLof11 <- 48
nchrLof14 <- 44

nboot <- 1000


####################
# load functions
require(data.table)
require(boot) # for bootstrap CIs


# For Tajima's D calcs. After https://github.com/ANGSD/angsd/blob/master/misc/stats.cpp
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

# calc thetas
calcthetas <- function(dat, nchr, nloci){
	# calc ave theta
	thetas <- as.numeric(dat[, .(tW = sum(exp(Watterson)/nloci, na.rm = TRUE), tP = sum(exp(Pairwise)/nloci, na.rm = TRUE))])

	# calculate Tajima's D
	thetas[3] <- tajd(nchr, thetas[1], thetas[2])
	
	#return
	names(thetas) <- c('tW', 'tP', 'tD')
	return(thetas)
}

# calculate stats from specified LGs for block bootstrapping across LGs
thetablock <- function(lgs, indices, alldata, nchr, regs){
	# make bootstrapped dataset
	mydata <- do.call("rbind", lapply(indices, function(n) subset(alldata, Chromo==lgs[n])))
	
	# calculate number of callable loci, given the LGs in this bootstrapped sample
	nloci <- regs[Chromo %in% lgs, .(len = Pos2 - Pos1 + 1), by = .(Chromo, Pos1)][, sum(len)] # sum of bp in the callable region
	
	# calc thetas
	thetas <- calcthetas(mydata, nchr, nloci)
	
	# return
	return(thetas)
}



######################
# Load data
######################

# load all loci theta calcs
datCan40 <- fread('analysis/thetas.Can_40.pestPG.gz')
datCan14 <- fread('analysis/thetas.Can_14.pestPG.gz')
datLof07 <- fread('analysis/thetas.Lof_07.pestPG.gz')
datLof11 <- fread('analysis/thetas.Lof_11.pestPG.gz')
datLof14 <- fread('analysis/thetas.Lof_14.pestPG.gz')

# gatk loci
datCan40gatk <- fread('analysis/thetas.Can_40.gatk.pestPG.gz')
datCan14gatk <- fread('analysis/thetas.Can_14.gatk.pestPG.gz')
datLof07gatk <- fread('analysis/thetas.Lof_07.gatk.pestPG.gz')
datLof11gatk <- fread('analysis/thetas.Lof_11.gatk.pestPG.gz')
datLof14gatk <- fread('analysis/thetas.Lof_14.gatk.pestPG.gz')

# fix name
setnames(datCan40, '#Chromo', 'Chromo')
setnames(datCan14, '#Chromo', 'Chromo')
setnames(datLof07, '#Chromo', 'Chromo')
setnames(datLof11, '#Chromo', 'Chromo')
setnames(datLof14, '#Chromo', 'Chromo')

setnames(datCan40gatk, '#Chromo', 'Chromo')
setnames(datCan14gatk, '#Chromo', 'Chromo')
setnames(datLof07gatk, '#Chromo', 'Chromo')
setnames(datLof11gatk, '#Chromo', 'Chromo')
setnames(datLof14gatk, '#Chromo', 'Chromo')

# list of callable regions
regs <- fread('data_2020.05.07/Callable_bases_gadmor2.bed')
setnames(regs, c('Chromo', 'Pos1', 'Pos2'))

# list of no damage sites
nodam <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab')
setnames(nodam, c('Chromo', 'Pos', 'REF', 'ALT'))


# remove unplaced
datCan40 <- datCan40[grep('Unplaced', Chromo, invert = TRUE), ]
datCan14 <- datCan14[grep('Unplaced', Chromo, invert = TRUE), ]
datLof07 <- datLof07[grep('Unplaced', Chromo, invert = TRUE), ]
datLof11 <- datLof11[grep('Unplaced', Chromo, invert = TRUE), ]
datLof14 <- datLof14[grep('Unplaced', Chromo, invert = TRUE), ]

datCan40gatk <- datCan40gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datCan14gatk <- datCan14gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datLof07gatk <- datLof07gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datLof11gatk <- datLof11gatk[grep('Unplaced', Chromo, invert = TRUE), ]
datLof14gatk <- datLof14gatk[grep('Unplaced', Chromo, invert = TRUE), ]

regs <- regs[Chromo != 'Unplaced', ]
nodam <- nodam[Chromo != 'Unplaced', ]

###################################################
# create table of loci trimmed to no damage sites
###################################################

datCan40gatknd <- merge(datCan40, nodam[, .(Chromo, Pos)], by = c('Chromo', 'Pos')) # trim
datCan14gatknd <- merge(datCan14, nodam[, .(Chromo, Pos)], by = c('Chromo', 'Pos')) # trim
datLof07gatknd <- merge(datLof07, nodam[, .(Chromo, Pos)], by = c('Chromo', 'Pos')) # trim
datLof11gatknd <- merge(datLof11, nodam[, .(Chromo, Pos)], by = c('Chromo', 'Pos')) # trim
datLof14gatknd <- merge(datLof14, nodam[, .(Chromo, Pos)], by = c('Chromo', 'Pos')) # trim


################################
# Run theta calculations
################################
nloci <- regs[, .(len = Pos2 - Pos1 + 1), by = .(Chromo, Pos1)][, sum(len)] # sum of bp in the callable region

# all loci
calcthetas(datCan40, nchrCan40, nloci)
calcthetas(datCan14, nchrCan14, nloci)
calcthetas(datLof07, nchrLof07, nloci)
calcthetas(datLof11, nchrLof11, nloci)
calcthetas(datLof14, nchrLof14, nloci)

# gatk loci
calcthetas(datCan40gatk, nchrCan40, nloci)
calcthetas(datCan14gatk, nchrCan14, nloci)
calcthetas(datLof07gatk, nchrLof07, nloci)
calcthetas(datLof11gatk, nchrLof11, nloci)
calcthetas(datLof14gatk, nchrLof14, nloci)

# gatk no damage loci
calcthetas(datCan40gatknd, nchrCan40, nloci)
calcthetas(datCan14gatknd, nchrCan14, nloci)
calcthetas(datLof07gatknd, nchrLof07, nloci)
calcthetas(datLof11gatknd, nchrLof11, nloci)
calcthetas(datLof14gatknd, nchrLof14, nloci)




# block bootstrapping across LGs
lgs <- datCan40[, sort(unique(Chromo))]
datlist <- list(datCan40, datCan14, datLof07, datLof11, datLof14, 
				datCan40gatk, datCan14gatk, datLof07gatk, datLof11gatk, datLof14gatk,
				datCan40gatknd, datCan14gatknd, datLof07gatknd, datLof11gatknd, datLof14gatknd)
names(datlist) <- c('Can40 all loci', 'Can14 all loci', 'Lof07 all loci', 'Lof11 all loci', 'Lof14 all loci', 
					'Can40 gatk loci', 'Can14 gatk loci', 'Lof07 gatk loci', 'Lof11 gatk loci', 'Lof14 gatk loci',
					'Can40 gatk no dam loci', 'Can14 gatk no dam loci', 'Lof07 gatk no dam loci', 'Lof11 gatk no dam loci', 
					'Lof14 gatk no dam loci')
nchrlist <- list(nchrCan40, nchrCan14, nchrLof07, nchrLof11, nchrLof14, nchrCan40, nchrCan14, nchrLof07, nchrLof11, nchrLof14, nchrCan40, nchrCan14, nchrLof07, nchrLof11, nchrLof14)

thetabootout <- data.frame(type = names(datlist), tW = NA, tWl95 = NA, tWu95 = NA, tP = NA, tPl95 = NA, tPu95 = NA, tD = NA, tDl95 = NA, tDu95 = NA)

for(i in 1:length(datlist)){
	print(names(datlist)[i])
	bootlg <- boot(lgs, thetablock, nboot,  alldata = datlist[[i]], nchr = nchrlist[[i]], regs = regs)
	
	print(bootlg)
	ciW <- boot.ci(bootlg, type = c('perc'), index = 1)
	ciP <- boot.ci(bootlg, type = c('perc'), index = 2)
	ciD <- boot.ci(bootlg, type = c('perc'), index = 3)
	
	thetabootout$tW[i] <- bootlg$t0[1] # the point estimates
	thetabootout$tP[i] <- bootlg$t0[2]	
	thetabootout$tD[i] <- bootlg$t0[3]

	thetabootout$tWl95[i] <- ciW$percent[4] # the confidence intervals
	thetabootout$tWu95[i] <- ciW$percent[5]

	thetabootout$tPl95[i] <- ciP$percent[4]
	thetabootout$tPu95[i] <- ciP$percent[5]

	thetabootout$tDl95[i] <- ciD$percent[4]
	thetabootout$tDu95[i] <- ciD$percent[5]
}

# save
write.csv(thetabootout, file = 'analysis/thetas.boot.cis.csv')
