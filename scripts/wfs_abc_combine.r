# run after wfs_abc.r or wfs_abc_function.r (latter called by wfs_abc_sbatch.sh or wfs_abc_sbatch_3times.sh, which in turn are called by script output by wfs_abc_makesbatchscript.r)
# combines the posterior samples into one file

# parameters
pop <- 'Lof'; myyr1 <- '07'; myyr2 <- '1114' # all 3 time points
#pop <- 'Can'; myyr1 <- '40'; myyr2 <- '13'

################################
# load functions and prep data
################################
if(!grepl('hpc.uio.no', Sys.info()["nodename"])){ # on something other than cod node
	require(RColorBrewer)
	require(data.table)
	require(boa)
	require(hexbin)
	require(plyr) # for summarize
}
if(grepl('hpc.uio.no', Sys.info()["nodename"])){ # on cod node
	require(RColorBrewer, lib.loc="/projects/cees/lib/R_packages/")
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(boa, lib.loc="/projects/cees/lib/R_packages/")
	require(plyr, lib.loc="/projects/cees/lib/R_packages/")
}


####################
# load data
####################
if(pop == 'Lof'){
	if(myyr2 != '1114'){
		targfile <- paste('data_2017.11.24/Frequency_table_Lof', myyr1, '_Lof', myyr2, '.txt', sep='')
	}
	if(myyr2 == '1114'){
		targfile <- paste('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', sep='')
		targfile2 <- paste('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', sep='')
	}
}
if(pop == 'Can'){
	targfile <- paste('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', sep='')
}
targ <- fread(targfile, header=TRUE)
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locus:=1:nrow(targ)] # add a locus number indicator

if(myyr2!='1114'){
	targ <- targ[,.(locus, CHROM, POS, f1samp, f2samp)]
}
if(myyr2=='1114'){
	targ2 <- fread(targfile2, header=TRUE)
	setnames(targ2, 3:7, c('alcnt1', 'f1samp', 'alcnt3', 'f3samp', 'ABS_DIFF13'))
	setkey(targ, CHROM, POS)
	setkey(targ2, CHROM, POS)
	targ <- targ[targ2,.(locus, CHROM, POS, f1samp, f2samp, f3samp)]
}


# trim to focal loci (outliers)
locs <- fread('analysis/outlier_annotation.csv') # the outliers
print(paste('Started with', nrow(targ), 'loci'))
if(pop == 'Lof'){
	locs <- locs[q3.Lof071114 !='' | q3.comb071114Can !='',.(CHROM, POS)]
	locs[, POS := as.numeric(POS)]
	setkey(locs, CHROM, POS)
	setkey(targ, CHROM, POS)
	targ <- merge(targ, locs)
}
if(pop == 'Can'){
	locs <- locs[q3.Can !='' | q3.comb071114Can !='',.(CHROM, POS)]
	locs[, POS := as.numeric(POS)]
	setkey(locs, CHROM, POS)
	setkey(targ, CHROM, POS)
	targ <- merge(targ, locs)
}
print(paste('Trimmed to', nrow(targ), 'loci by using only outliers in analysis/outlier_annotation.csv'))

# prep datatable to hold the posterior samples
if(myyr2!='1114'){
	posts <- data.table(CHROM=character(), POS=numeric(), f1samp=numeric(), f2samp=numeric(), ne=numeric(), f1=numeric(), f2=numeric(), s=numeric(), h=numeric())
}
if(myyr2=='1114'){
	posts <- data.table(CHROM=character(), POS=numeric(), f1samp=numeric(), f2samp=numeric(), f3samp=numeric(), ne=numeric(), f1=numeric(), f2=numeric(), f3=numeric(), s=numeric(), h=numeric())
}

#########################################################
# combine posterior samples
# loop over all the posterior files output by wfs_abc.r
#########################################################

# find files to read in
files <- list.files(path='analysis/temp', pattern='wfs_abc_sampsize*', full.names=TRUE)
print(paste(length(files), 'to process'))

# read in each file and extract posterior samples
for(i in 1:length(files)){
	print(paste(i, 'of', length(files), ':', files[i]))
	theseposts <- fread(paste('zcat', files[i]))
	theselocs <- sort(unique(theseposts$locus))
	inds <- which(targ$locus %in% theselocs)
	
	if(length(inds)>0){
		origlen <- nrow(theseposts)
		theseposts <- merge(theseposts, targ[,.(CHROM, POS, locus)], by='locus')
		
		if(nrow(theseposts) != origlen){
			warning(paste('# rows did not match for file', files[i]))
		} else {
			theseposts[,locus:=NULL] # remove column
			posts <- rbind(posts, theseposts)
		}
		
	} else {
		warning(paste('WARNING: Could not find locus(s)', paste(theselocs, collapse=', ')))
	}
}

# sample size per locus
posts[,.N, by=c('CHROM', 'POS')]


# write out
outfile <- paste('analysis/wfs_abc_posts_', pop, '_', myyr1, '-', myyr2, '.rds', sep='')
outfile
saveRDS(posts, file=outfile)



