# Fit a brt model to pre- and post- fishing loci
# Run this from the command line on a cod node: nohup R CMD BATCH scripts/polygenic_brt_fit.r &

Sys.time() # when did we start?

# set arguments
usecorrected <- FALSE; suff<-'_nocorr' # whether or not to use population-corrected phenotypes and genotypes
#usecorrected <- TRUE; suff<-'_corr' # whether or not to use population-corrected phenotypes and genotypes

usefactor <- TRUE; suff<-paste(suff,'_fact',sep='') # only makes sense with usecorrected=FALSE
#usefactor <- FALSE; suff<-paste(suff,'_nofact',sep='')

interaction.depth <- c(2,6,12); n.trees <- seq(500,10000,by=500); shrinkage <- c(0.01, 0.001); n.minobsinnode <- 10; p.train <- 0.66; ncv <- 10; nrep <- 3; ncores <- 30 # big
#interaction.depth <- c(2); n.trees <- 100; shrinkage <- 0.01; n.minobsinnode <- 10; p.train <- 0.66; ncv <- 3; nrep <- 1; ncores <- 2 # test

# load required packages
.libPaths('/projects/cees/lib/R_packages/') # set package locations on cod node
require(caret, lib.loc="/projects/cees/lib/R_packages/")
require(gbm, lib.loc="/projects/cees/lib/R_packages/")
require(e1071, lib.loc="/projects/cees/lib/R_packages/") # for factor response variable
require(doSNOW, lib.loc="/projects/cees/lib/R_packages/") # set up SOCK cluster

# load data
if(usecorrected){ # load the corrected data
	dat <- readRDS('analysis/polygenic_test_genos_corr.rds')
}
if(!usecorrected){ # load the uncorrected data
	dat <- as.data.frame(readRDS('analysis/polygenic_test_genos.rds'))
}

# Use factor for response?
if(usefactor){
	dat$fish <- as.factor(dat$fish)
}

	dat[1:10,1:5] # examine to make sure we got the right one
	summary(dat[,1:5]) # examine to make sure formatted correctly

# Divide into training and testing datasets
# ALSO STRATIFY BY POPULATION?
inTrain <- createDataPartition(y = dat$fish, p = p.train, list = FALSE)
dat_train <- dat[inTrain,] # make the training dataset

nrow(dat_train) # size of the training dataset
nrow(dat) - nrow(dat_train) # size of the testing dataset

# Define the tuning parameter grid over which to search (number of interaction levels, number of trees, and learning rate)
fitGrid <- expand.grid(interaction.depth=interaction.depth, n.trees=n.trees, shrinkage=shrinkage, n.minobsinnode=n.minobsinnode)

nrow(fitGrid) # number of parameter sets to analyze

# Define the tuning and training methods
fitControl <- trainControl(method="repeatedcv", number=ncv, repeats=nrep) 

# Set up a cluster
cl <- makeCluster(rep('localhost', ncores), type='SOCK') # make the cluster on the localhost
registerDoSNOW(cl) # register the cluster so accessible to foreach
clusterCall(cl, function(x) .libPaths(x), .libPaths()) # add libpaths to each node	

# Train BRT models with each combination of tuning parameters
# Uses a default bag fraction of 0.5
brtfit <- train(x=dat_train[,3:ncol(dat)], y=dat_train$fish, method="gbm", tuneGrid=fitGrid, trControl=fitControl, verbose=TRUE)
	
# print some basic stats
brtfit

# save the model
saveRDS(brtfit, file=paste('analysis/polygenic_test_brtfit', suff, '.rds', sep=''))

#shut down cluster
stopCluster(cl)

Sys.time() # when did we end?
