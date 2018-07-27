# Trying popgenome R package to calculate diversity stats


	# cod node
	require(data.table, lib.loc="/projects/cees/lib/R_packages/")
	require(parallel, lib.loc="/projects/cees/lib/R_packages/")
	require(iterators, lib.loc="/projects/cees/lib/R_packages/") # used with foreach
	require(foreach, lib.loc="/projects/cees/lib/R_packages/") # for parallel procssing using foreach()
	require(doMC, lib.loc="/projects/cees/lib/R_packages/")
	require(doParallel, lib.loc="/projects/cees/lib/R_packages/") # needed for foreach
	require(doSNOW, lib.loc="/projects/cees/lib/R_packages/") # set up SOCK cluster


#########################################################################################
# read in data, one lg at a time, and concat together
# need tabix-ed vcf file e.g., run this on a cod node: module load tabix; tabix -p vcf All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz
#########################################################################################
# on my mac
require(PopGenome)
require(data.table)

# settings
lgs <- c('LG03', 'LG04', 'LG05', 'LG06', 'LG08', 'LG09', 'LG10', 'LG11', 'LG13', 'LG14', 'LG15', 'LG16', 'LG17', 'LG18', 'LG19', 'LG20', 'LG21', 'LG22', 'LG23') # skipping the inversions
frompos <- 1
topos <- 50000000; # more than largest LG

# read in kmer and depth statistics
# depth stat only calculated on loci that passed kmer25 filter, so just use that
dpLof <- fread('analysis/Outlier_sequencing_depth_statistic_All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz.txt')
dpCan <- fread('analysis/Outlier_sequencing_depth_statistic_All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz.txt')
	setnames(dpLof, c('CHROM', 'POS', 'dpstatLof'))
	setnames(dpCan, c('CHROM', 'POS', 'dpstatCan'))
depthstat <- merge(dpLof, dpCan, all=FALSE) # only care about the intersecting loci
	nrow(dpLof)
	nrow(dpCan)
	nrow(depthstat)
	
depthstat[,dpFlag := (dpstatLof < quantile(dpstatLof, na.rm=TRUE, probs=0.95)) & (dpstatCan < quantile(dpstatCan, na.rm=TRUE, probs=0.95))] # TRUE if not a depth outlier in either population (outlier are top 5%)
	depthstat[is.na(dpFlag), dpFlag := FALSE] # FALSE if locus not genotyped
	depthstat[,sum(dpFlag)] # number of loci genotyped in both that pass kmer and depth filter in both populations
depthstat <- depthstat[dpFlag==TRUE,] # trim just to loci that pass

# read in first lg
vnor <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid=lgs[1], approx=FALSE, gffpath='../genome_data/gadMor2_annotation_filtered_only_gene_models.gff', frompos=frompos, topos=topos, numcols=100000, include.unknown=TRUE)

vcan <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz', tid=lgs[1], approx=FALSE, gffpath='../genome_data/gadMor2_annotation_filtered_only_gene_models.gff', frompos=frompos, topos=topos, numcols=100000, include.unknown=TRUE)

# find the list of shared SNPs that pass kmer and depth filters
#g@region.data@biallelic.matrix[[1]][1:10,1:20] # look at the matrix of major and minor alleles. Two rows for each individual
locnor <- colnames(vnor@region.data@biallelic.matrix[[1]])
loccan <- colnames(vcan@region.data@biallelic.matrix[[1]])
locshared <- intersect(locnor, loccan)
locpass <- intersect(locshared, depthstat[CHROM==lgs[1], as.numeric(POS)])
	length(locpass)
locindsnor <- locnor %in% locpass
locindscan <- loccan %in% locpass

# simplify to genotype matrix (0,1,2) of the shared SNPs
dims <- dim(vnor@region.data@biallelic.matrix[[1]])
temp <- array(vnor@region.data@biallelic.matrix[[1]][,locindsnor], dim=c(2,dims[1]/2, sum(locindsnor))) # reshape so first dimension is within individual
gnor <- apply(X=temp[], MARGIN=c(2,3), FUN=sum) # not clear why temp[] is needed instead of temp. takes a few seconds
rownames(gnor) <- rownames(vnor@region.data@biallelic.matrix[[1]])[seq(1,dims[1],by=2)]
colnames(gnor) <- paste(lgs[1], '_', colnames(vnor@region.data@biallelic.matrix[[1]])[locindsnor], sep='')
	dim(gnor)
	
dims <- dim(vcan@region.data@biallelic.matrix[[1]])
temp <- array(vcan@region.data@biallelic.matrix[[1]][,locindscan], dim=c(2,dims[1]/2, sum(locindscan))) # reshape so first dimension is within individual
gcan <- apply(X=temp[], MARGIN=c(2,3), FUN=sum) # not clear why temp[] is needed instead of temp. takes a few seconds
rownames(gcan) <- rownames(vcan@region.data@biallelic.matrix[[1]])[seq(1,dims[1],by=2)]
colnames(gcan) <- paste(lgs[1], '_', colnames(vcan@region.data@biallelic.matrix[[1]])[locindscan], sep='')
	dim(gcan)

# loop through rest of LGs
for(i in 2:length(lgs)){
	print(lgs[i])
	vnor <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid=lgs[i], approx=FALSE, gffpath='../genome_data/gadMor2_annotation_filtered_only_gene_models.gff', frompos=frompos, topos=topos, numcols=100000, include.unknown=TRUE)
	vcan <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz', tid=lgs[i], approx=FALSE, gffpath='../genome_data/gadMor2_annotation_filtered_only_gene_models.gff', frompos=frompos, topos=topos, numcols=100000, include.unknown=TRUE)
	locnor <- colnames(vnor@region.data@biallelic.matrix[[1]])
	loccan <- colnames(vcan@region.data@biallelic.matrix[[1]])
	locshared <- intersect(locnor, loccan)
	locpass <- intersect(locshared, depthstat[CHROM==lgs[i], as.numeric(POS)])
		print(length(locpass))
	locindsnor <- locnor %in% locpass
	locindscan <- loccan %in% locpass

	dims <- dim(vnor@region.data@biallelic.matrix[[1]])
	temp <- array(vnor@region.data@biallelic.matrix[[1]][,locindsnor], dim=c(2,dims[1]/2, sum(locindsnor))) # reshape so first dimension is within individual
	temp2 <- apply(X=temp[], MARGIN=c(2,3), FUN=sum) # not clear why temp[] is needed instead of temp. takes a few seconds
	rownames(temp2) <- rownames(vnor@region.data@biallelic.matrix[[1]])[seq(1,dims[1],by=2)]
	colnames(temp2) <- paste(lgs[i], '_', colnames(vnor@region.data@biallelic.matrix[[1]])[locindsnor], sep='')

	gnor <- cbind(gnor, temp2)

	dims <- dim(vcan@region.data@biallelic.matrix[[1]])
	temp <- array(vcan@region.data@biallelic.matrix[[1]][,locindscan], dim=c(2,dims[1]/2, sum(locindscan))) # reshape so first dimension is within individual
	temp2 <- apply(X=temp[], MARGIN=c(2,3), FUN=sum) # not clear why temp[] is needed instead of temp. takes a few seconds
	rownames(temp2) <- rownames(vcan@region.data@biallelic.matrix[[1]])[seq(1,dims[1],by=2)]
	colnames(temp2) <- paste(lgs[i], '_', colnames(vcan@region.data@biallelic.matrix[[1]])[locindscan], sep='')

	gcan <- cbind(gcan, temp2)

	print(dim(gnor))
	print(dim(gcan))
}

# combine norway and canada
genos <- rbind(gnor, gcan)
	summary(as.numeric(genos))
	dim(genos)

# turn NaNs to NAs (required for gbm)
inds <- is.nan(genos)
genos[inds] <- NA # slow

# define populations
poplist <- list(LOF07=c('BM_209', 'BM_211', 'BM_213', 'BM_214', 'BM_216', 'BM_217', 'BM_218', 'BM_219', 'BM_220', 'BM_221', 'BM_222', 
'BM_223', 'BM_224', 'BM_225', 'BM_226', 'BM_227', 'BM_230', 'BM_231', 'BM_232', 'BM_234', 'BM_236', 'BM_237', 'BM_239', 'BM_240', 'BM_241', 'BM_242', 'BM_243', 'BM_244'), LOF11=c('LOF1103001', 'LOF1103002', 'LOF1103003', 'LOF1103004', 'LOF1103005', 'LOF1103006', 'LOF1103007', 'LOF1103008', 'LOF1103009', 'LOF1103010', 'LOF1103011', 'LOF1103012', 'LOF1103013', 'LOF1103014', 'LOF1103015', 'LOF1103016', 'LOF1103017', 'LOF1103018', 'LOF1103019', 'LOF1103020', 'LOF1103021', 'LOF1103022', 'LOF1103023', 'LOF1103024'), LOF14=c('LOF1403001', 'LOF1403002', 'LOF1403003', 'LOF1403004', 'LOF1403005', 'LOF1403006', 'LOF1403007', 'LOF1403008', 'LOF1403009', 'LOF1403010', 'LOF1403011', 'LOF1403012', 'LOF1403013', 'LOF1403014', 'LOF1403015', 'LOF1403016', 'LOF1403017', 'LOF1403018', 'LOF1403019', 'LOF1403020', 'LOF1403021', 'LOF1403022', 'LOF1403023', 'LOF1403024'), CAN40=c('BM_026', 'BM_027', 'BM_028', 'BM_029', 'BM_032', 'BM_105', 'BM_106', 'BM_107', 'BM_109', 'BM_111', 'BM_112', 'BM_115', 'BM_116', 'BM_117', 'BM_118', 'BM_119', 'BM_123', 'BM_124', 'BM_125', 'BM_126', 'BM_127', 'BM_128'), CANMod=c('TWI1307001', 'TWI1307002', 'TWI1307003', 'TWI1307004', 'TWI1307005', 'TWI1307006', 'TWI1307007', 'TWI1307008', 'TWI1307009', 'TWI1307010', 'TWI1307011', 'TWI1307012', 'TWI1307013', 'TWI1307014', 'TWI1307015', 'TWI1307016', 'TWI1307017', 'TWI1307018', 'TWI1307019', 'TWI1307020', 'TWI1307021', 'TWI1307022', 'TWI1307023', 'TWI1307024'))

pops <- rep(NA, nrow(genos))
pops[rownames(genos) %in% poplist[['LOF07']]] <- 'LOF'
pops[rownames(genos) %in% poplist[['LOF11']]] <- 'LOF'
pops[rownames(genos) %in% poplist[['LOF14']]] <- 'LOF'
pops[rownames(genos) %in% poplist[['CAN40']]] <- 'CAN'
pops[rownames(genos) %in% poplist[['CANMod']]] <- 'CAN'
pops <- as.factor(pops)

# define pre and post fishing
fish <- rep(NA, nrow(genos))
fish[rownames(genos) %in% poplist[['LOF07']]] <- 1
fish[rownames(genos) %in% poplist[['LOF11']]] <- 2
fish[rownames(genos) %in% poplist[['LOF14']]] <- 2
fish[rownames(genos) %in% poplist[['CAN40']]] <- 1
fish[rownames(genos) %in% poplist[['CANMod']]] <- 2

# create data frame with response and explanatory variables
dat <- cbind(pops, fish, genos)



# write out genotype matrix
saveRDS(dat, file='analysis/polygenic_test_genos.rds')


###################################
# Correct for population structure
###################################
# read in genotype matrix
dat <- data.frame(readRDS('analysis/polygenic_test_genos.rds'))

# Logging function for the parallel foreach loop
# to listen, first run "nc -l 4000" at the terminal
log.socket <- make.socket(port=4000)
Log <- function(text, ...) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  write.socket(log.socket, msg)
}

# Correct for population stratification before conducting RF to minimize the risk of false positive associations.
# Specifically, we will correct the genotypes and phenotypes using the approach of Zhao et al.(2012)
# Now correct the genotypes and phenotypes using the regression/residual method. We're using a standard linear regression because Zhao et al. 2012 found that the correction procedure is robust to selection of the link function

	# my mac
#	library(doMC)
#	registerDoMC(cores = 3)
	
	# cod node
	.libPaths('/projects/cees/lib/R_packages/')
	cl <- makeCluster(rep('localhost', 30), type='SOCK') # make the cluster on the localhost
	registerDoSNOW(cl) # register the cluster so accessible to foreach
	clusterCall(cl, function(x) .libPaths(x), .libPaths()) # add libpaths to each node	

n <- ncol(dat)
datcor <- foreach (d=iter(dat[,2:ncol(dat)], by='column'), i=icount(), .combine=cbind) %dopar% { # may need , .packages=c('doParallel')
	lmi <- lm(d ~ dat$pops, na.action=na.exclude) # na.exclue important for padding NAs
	return(residuals(lmi)) # apply linear model to all loci and the response. return residuals padded with residuals()
#	if(i%%50==0) print(i)
	Log("Processing %d of %d", i, n) # listen by first running "nc -l 4000" at the terminal
}
	dim(datcor)

stopCluster(cl)

datcor <- cbind(dat$pops, data.frame(datcor))
colnames(datcor) <-colnames(dat)[1:ncol(datcor)]

# Export a copy of the corrected genotypes and phenotypes for future reference (This corresponds to the data in Table S6)
saveRDS(datcor,file="analysis/polygenic_test_genos_corr.rds")



#########################
# Fit BRT
#########################

# see polygenic_brt_fit.r

#########################
# Examine BRT
#########################
require(caret)
require(gbm)
require(pdp) # for partial dependency plots


# load model
brtfit <- readRDS('analysis/polygenic_test_brtfit_nocorr_fact.rds') # for classification model on non-population-corrected data
#brtfit <- readRDS('analysis/polygenic_test_brtfit_corr_nofact.rds') # for regression model on population-corrected data

# load data
dat <- readRDS('analysis/polygenic_test_genos.rds'); dat$fish <- as.factor(dat$fish)
#dat <- readRDS('analysis/polygenic_test_genos_corr.rds')

# trim data to training/testing
inTrain <- rownames(dat) %in% rownames(brtfit$trainingData)
dat_train <- dat[inTrain,]
dat_test <- dat[!inTrain,]
	nrow(dat_train)
	nrow(dat_test)
	
# examine model
brtfit

### Examine how model fit (measured as Root Mean Square Error) varied with the tuning parameters.
ggplot(brtfit)


# See which explanatory variables were most important. Each level of the region factor is reported separately, but you can effectively consider them all the same variable.
varImp(brtfit) # scales most important to 100 (caret package)

varimp <- summary(brtfit, plotit=FALSE) # returns all the variables and their relative importances

	i = order(varimp$rel.inf)
	subs <- c(seq(1, 100000, by=1000), seq(100001, 170000, by=500), 170001:nrow(varimp)) # don't plot everything
	plot((1:nrow(varimp))[subs], varimp$rel.inf[i][subs])
	subs <- varimp$rel.inf[i] > 0
		sum(subs)
	plot((1:nrow(varimp))[subs], varimp$rel.inf[i][subs])
	plot((1:nrow(varimp))[subs], varimp$rel.inf[i][subs], xlim=c(170820, nrow(varimp))) # upper edge of elbow
	plot((1:nrow(varimp))[subs], varimp$rel.inf[i][subs], log='y')
	
	thresh <- quantile(varimp$rel.inf, probs=0.95)
	sum(varimp$rel.inf >= thresh)

	subs <- varimp$rel.inf[i] > 0
	plot((1:nrow(varimp))[subs], cumsum(varimp$rel.inf[i][subs]), cex=0.2)
		abline(h=20)

# Visualize the BRT fits for one variable. Slow.
partial(brtfit, pred.var='X17836249', plot=TRUE, rug=TRUE)


# partial(brtfit, pred.var=c('SBT.max', 'GRAINSIZE'), plot=TRUE, rug=FALSE)

# Test the model against the data we held aside
preds <- predict(brtfit, dat_test) # figure out the model predictions, given the explanatory variables in the testing dataset.
predsNor <- predict(brtfit, dat_test[dat_test$pops=='LOF',])
predsCan <- predict(brtfit, dat_test[dat_test$pops=='CAN',])
predstrain <- predict(brtfit, dat_train) # training data
	
	# regression
	postResample(pred=preds, obs=dat_test$fish) # compare the predictions and the observations with RMSE, R2, and Mean Absolute Error.
	postResample(pred=preds, obs=dat_test$fish[dat_test$pops=='LOF'])
	postResample(pred=preds, obs=dat_test$fish[dat_test$pops=='CAN'])
	postResample(pred=predstrain, obs=dat_train$fish)
