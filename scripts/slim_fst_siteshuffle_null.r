# shuffle FST A and B values across sites and calculate windowed FST to get a null distribution of max genome-wide FST from SLiM simulations
# run slim_calcfst_reynolds.sh first to get the *.fst.csv.gz files
# also used by slim_fst_siteshuffle_null.sh
# called as
#	Rscript slim_fst_siteshuffle_null.r infile outfile

# to run on saga

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
infile <- args[1]
outfile <- args[2]


# parameters
winsz <- 50000 # window size
winstp <- 10000 # window step
nrep <- 1000 # number of reshuffles
minloci <- 2 # minimum number of loci per window to consider


# load functions
require(data.table)


#############
# Prep data
#############
# load fst A/B data
dat <- fread(infile)



# create new columns as indices for windows
for(j in 1:(winsz/winstp)){
	dat[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
}

# mark windows with < minloci for removal
rem <- 0 # number of windows removed
for(j in 1:(winsz/winstp)){
	datwin <- dat[, .(nsnps = length(POS)), by = .(win = get(paste0('win', j)))] # calc num snps per window
	rem[1] <- rem[1] + datwin[, sum(nsnps < minloci)] # record number to be removed
	datwin[, (paste0('win', j, 'keep')) := 1] # create col to mark which windows to keep
	datwin[nsnps < minloci, (paste0('win', j, 'keep')) := 0] # mark windows to remove
	datwin[, nsnps := NULL] # drop column
	setnames(datwin, "win", paste0('win', j)) # change col name
	dat <- merge(dat, datwin, by = paste0('win', j), all.x = TRUE) # merge keeper col back to full dataset

}

rem # number of windows removed for each comparison




####################################
# shuffle and recalc windowed FST
####################################
colnms <- c('CHROM', 'POS', paste0('win', 1:(winsz/winstp)), paste0('win', 1:(winsz/winstp), 'keep')) # list of column names we want out of the base data.table

for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(dat), nrow(dat), replace = FALSE)
	temp <- cbind(dat[, ..colnms], dat[inds, .(A, B)]) # shuffle FSTs across positions
		
	# calc fst for each window to keep
	for(j in 1:(winsz/winstp)){
		temp2 <- temp[get(paste0('win', j, 'keep')) == 1, ] # trim to windows to keep. can't combine with next line for some reason.
		if(j ==1) tempfsts <- temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst from the genome
	# exclude windows with negative midpoints
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfile), row.names = FALSE)
