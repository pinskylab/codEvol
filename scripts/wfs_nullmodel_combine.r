# Calculate probability of null model producing results as extreme as our observations
# run after wfs_nullmodel_function.r
# Script combines the output
# To be run on cod node


# load functions: assume this is run on a cod or abel node
require(data.table, lib.loc="/projects/cees/lib/R_packages/")


# load data
#targ <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); suffix='_07-14' # Lof 1907-2014
#targ <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); suffix='_07-11' # Lof 1907-2011
#targ <- fread('data_2017.11.24/Frequency_table_Lof11_Lof14.txt', header=TRUE); suffix='_11-14' # Lof 2011-2014
targ <- fread('data_2017.11.24/Frequency_table_Lof07_Lof11.txt', header=TRUE); targ2 <- fread('data_2017.11.24/Frequency_table_Lof07_Lof14.txt', header=TRUE); suffix='_07-11-14' # Lof 1907-2011-2014
#targ <- fread('data_2017.11.24/Frequency_table_CAN_40_TGA.txt', header=TRUE); suffix='_Can' # Can
#targ <- fread('analysis/Frequency_table_PowerSims_Lof_Ne46000_cnt46_44.txt', header=TRUE); suffix='_Power_Lof_Ne46000_cnt46_44' # Lof power analysis
#targ <- fread('analysis/Frequency_table_PowerSims_Can_Ne5900_cnt32_40.txt', header=TRUE); suffix='_Power_Can_Ne5900_cnt32_40' # Can power analysis
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locusnum:=1:nrow(targ)] # add a locus number indicator

if(suffix=='_07-11-14'){
	setnames(targ2, 3:7, c('alcnt1', 'f1samp', 'alcnt3', 'f3samp', 'ABS_DIFF13'))
	setkey(targ, CHROM, POS)
	setkey(targ2, CHROM, POS)
	targ <- targ[targ2,.(locusnum, CHROM, POS, alcnt1, alcnt2, alcnt3, f1samp, f2samp, f3samp)]
}

# how many loci at each sample size
if(suffix != '_07-11-14'){
	setkey(targ, alcnt1, alcnt2)
	nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2)]
}
if(suffix == '_07-11-14'){
	setkey(targ, alcnt1, alcnt2, alcnt3)
	nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2, alcnt3)]
}
	nrow(nloci) # 1907-2014: 30
				# 1907-2011: 
				# 2011-2014: 
				# 1907-2011-2014: 1145
				# Can: 90

# read in files
if(suffix != '_07-11-14') dat <- data.frame(cnt1=numeric(0), cnt2=numeric(0), locusnum=character(0), p=numeric(0), n=numeric(0))
if(suffix == '_07-11-14') dat <- data.frame(cnt1=numeric(0), cnt2=numeric(0), cnt3=numeric(0), locusnum=character(0), p=numeric(0), n=numeric(0))
for(i in 1:nrow(nloci)){
	if(i %% 25 == 0) print(i)
	# check if the simulation files for this sample size exists
	myalcnt1 <- nloci[i,alcnt1]
	myalcnt2 <- nloci[i,alcnt2]
	if(suffix == '_07-11-14') myalcnt3 <- nloci[i,alcnt3]
	
	if(suffix != '_07-11-14') thesefiles <- list.files(path='analysis/temp', pattern=paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''), full.names=TRUE)
	if(suffix == '_07-11-14') thesefiles <- list.files(path='analysis/temp', pattern=paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, myalcnt3, sep=','), '_locus*', sep=''), full.names=TRUE)

	if(length(thesefiles)>0){
		for(j in 1:length(thesefiles)){
			thisdat <- read.csv(thesefiles[j])
			if(ncol(thisdat) == 3){
				thisdat$cnt1 <- myalcnt1
				thisdat$cnt2 <- myalcnt2
				if(suffix == '_07-11-14') thisdat$cnt3 <- myalcnt3
				dat <- rbind(dat,thisdat)
			} else {
				if(nloci[i,nloci]==1){ # file organized as one column if only one locus
					thisdat <- data.frame(cnt1=myalcnt1, cnt2=myalcnt2, locusnum=thisdat[1,1], p=thisdat[2,1], n=thisdat[3,1])
					if(suffix == '_07-11-14') thisdat$cnt3 <- myalcnt3
					dat <- rbind(dat,thisdat)
				} else {
					print(paste('>1 locus and wrong number of columns for i=', i, sep=''))
				}
			}
		}
	} else {
		if(suffix != '_07-11-14') print(paste('Missing null model analysis for', paste(myalcnt1, myalcnt2, sep=',')))
		if(suffix == '_07-11-14') print(paste('Missing null model analysis for', paste(myalcnt1, myalcnt2, myalcnt3, sep=',')))
	}
}

# sort by locus number
dat <- dat[order(dat$locusnum),]

# all loci analyzed?
dim(targ)
nrow(targ)
nrow(targ) - nrow(dat) # number missing loci

# NA values?
sum(is.na(dat$p)) # 0

# adjusted p-values
dat$p.adj <- p.adjust(dat$p, method='fdr')


# write out
outfile <- paste('analysis/wfs_nullmodel_pvals', suffix, '.rdata', sep='')
outfile
save(dat, file=outfile)


