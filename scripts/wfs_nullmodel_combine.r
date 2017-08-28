# Calculate probability of null model producing results as extreme as our observations
# run after wfs_nullmodel_function.r
# Script combines the output
# To be run on cod node


# load functions: assume this is run on a cod or abel node
require(data.table, lib.loc="/projects/cees/lib/R_packages/")


# load data
#targ <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_25k.txt', header=TRUE); suffix='_07-14_25k' # Lof 1907-2014 25kmer
#targ <- fread('data_29.06.17/Frequency_table_Lof07_Lof14_150k.txt', header=TRUE); suffix='_07-14_150k' # Lof 1907-2014 150k
#targ <- fread('data_29.06.17/Frequency_table_Lof07_Lof11_25k.txt', header=TRUE); suffix='_07-11_25k' # Lof 1907-2011 25k
#targ <- fread('data_29.06.17/Frequency_table_Lof11_Lof14_25k.txt', header=TRUE); suffix='_11-14_25k' # Lof 2011-2014 25k
#targ <- fread('data_11.07.17/Frequency_table_Can_40_Can_25k.txt', header=TRUE); suffix='_Can_25k' # Can 25k
#targ <- fread('data_11.07.17/Frequency_table_Can_40_Can_150k.txt', header=TRUE); suffix='_Can_150k' # Can 150k
targ <- fread('analysis/Frequency_table_PowerSims_Lof_Ne46000_cnt46_44.txt', header=TRUE); suffix='_Power_Lof_Ne46000_cnt46_44' # Can 150k
setnames(targ, 3:7, c('alcnt1', 'f1samp', 'alcnt2', 'f2samp', 'ABS_DIFF'))
targ[,locusnum:=1:nrow(targ)] # add a locus number indicator

# how many loci at each sample size
setkey(targ, alcnt1, alcnt2)
nloci <- targ[,.(nloci=length(locusnum)), by=.(alcnt1, alcnt2)]
	nrow(nloci) # 1907-2014: 176 for 25kmer, 193 for 150kmer
				# 1907-2011: 249 for 25kmer
				# 2011-2014: 162 for 25kmer
				# Can: 		55 for 25kmer, 55 for 150kmer

# read in files
dat <- data.frame(cnt1=numeric(0), cnt2=numeric(0), locusnum=character(0), p=numeric(0), n=numeric(0))
for(i in 1:nrow(nloci)){
	if(i %% 25 == 0) print(i)
	# check if the simulation files for this sample size exists
	myalcnt1 <- nloci[i,alcnt1]
	myalcnt2 <- nloci[i,alcnt2]
	
	thesefiles <- list.files(path='analysis/temp', pattern=paste('wfs_nullmodel_sampsize', paste(myalcnt1, myalcnt2, sep=','), '_locus*', sep=''), full.names=TRUE)

	if(length(thesefiles)>0){
		for(j in 1:length(thesefiles)){
			thisdat <- read.csv(thesefiles[j])
			if(ncol(thisdat) == 3){
				thisdat$cnt1 <- myalcnt1
				thisdat$cnt2 <- myalcnt2
				dat <- rbind(dat,thisdat)
			} else {
				if(nloci[i,nloci]==1){ # file organized as one column if only one locus
					thisdat <- data.frame(cnt1=myalcnt1, cnt2=myalcnt2, locusnum=thisdat[1,1], p=thisdat[2,1], n=thisdat[3,1])
					dat <- rbind(dat,thisdat)
				} else {
					print(paste('>1 locus and wrong number of columns for i=', i, sep=''))
				}
			}
		}
	} else {
		print(paste('Missing null model analysis for', paste(myalcnt1, myalcnt2, sep=',')))
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

# create new p-value never ==0
i <- dat$p == 0
sum(i)
dat$pmax <- dat$p
dat$pmax[i] <- 1/dat$n[i] # set to highest possible based on sample size

# adjusted p-values
dat$p.adj <- p.adjust(dat$pmax, method='fdr')


# write out
outfile <- paste('analysis/wfs_nullmodel_pvals', suffix, '.rdata', sep='')
outfile
save(dat, file=outfile)


