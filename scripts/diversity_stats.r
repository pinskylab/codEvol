# Trying popgenome R package to calculate diversity stats

require(PopGenome)

# settings
lg <- 'LG05'; frompos <- 1; topos <- 24067848 # what to read in

# need tabix-ed vcf file
# cod node: module load tabix; tabix -p vcf All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz
g <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid=lg, approx=FALSE, gffpath='../genome_data/gadMor2_annotation_filtered_only_gene_models.gff', frompos=frompos, topos=topos, numcols=100000, include.unknown=TRUE)

#g <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz', tid='LG05', approx=FALSE, gffpath='../genome_data/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=1000000)

g

# look at data
head(g@region.data@biallelic.sites[[1]]) # location of SNPs
dim(g@region.data@biallelic.matrix[[1]]) # size of the matrix of major and minor alleles. Two rows for each individual
g@region.data@biallelic.matrix[[1]][1:10,1:20] # look at the matrix of major and minor alleles. Two rows for each individual

# define populations
g <- set.populations(g, list(LOF07=c('BM_209', 'BM_211', 'BM_213', 'BM_214', 'BM_216', 'BM_217', 'BM_218', 'BM_219', 'BM_220', 'BM_221', 'BM_222', 'BM_223', 'BM_224', 'BM_225', 'BM_226', 'BM_227', 'BM_230', 'BM_231', 'BM_232', 'BM_234', 'BM_236', 'BM_237', 'BM_239', 'BM_240', 'BM_241', 'BM_242', 'BM_243', 'BM_244'), LOF11=c('LOF1103001', 'LOF1103002', 'LOF1103003', 'LOF1103004', 'LOF1103005', 'LOF1103006', 'LOF1103007', 'LOF1103008', 'LOF1103009', 'LOF1103010', 'LOF1103011', 'LOF1103012', 'LOF1103013', 'LOF1103014', 'LOF1103015', 'LOF1103016', 'LOF1103017', 'LOF1103018', 'LOF110301', 'LOF1103020', 'LOF1103021', 'LOF1103022', 'LOF1103023', 'LOF1103024'), LOF14=c('LOF1403001', 'LOF1403002', 'LOF1403003', 'LOF1403004', 'LOF1403005', 'LOF1403006', 'LOF1403007', 'LOF1403008', 'LOF1403009', 'LOF1403010', 'LOF1403011', 'LOF1403012', 'LOF1403013', 'LOF1403014', 'LOF1403015', 'LOF1403016', 'LOF1403017', 'LOF1403018', 'LOF1403019', 'LOF1403020', 'LOF1403021', 'LOF1403022', 'LOF1403023', 'LOF1403024')), diploid=TRUE)

# set reference. not sure yet what this is used for
g <- set.synnonsyn(g, ref.chr='../genome_data/gadMor2.fasta', save.codons=TRUE)

# diversity stats by population
g <- diversity.stats(g)
get.diversity(g)
g@nuc.diversity.within

# try filtering by missing data
g2 <- set.filter(g, missing.freqs=TRUE, miss.lower.bound=0, miss.upper.bound=0.2)
length(g2@region.data@included[[1]])
sum(g2@region.data@included[[1]])

#############################
# sliding window analysis
#############################
# set up
winsz <- 10000
slide <- sliding.window.transform(g, width=winsz, 10000, type=2)

# filter out some data (could base on kmer25, GATK callability, etc.?)
slide <- set.filter(slide, missing.freqs=TRUE, miss.lower.bound=0, miss.upper.bound=0.2)

# number of windows
nwin <- length(slide@region.names)

# statistics
slide <- diversity.stats(slide)
nucdiv <- slide@nuc.diversity.within

# normalize by number of nucleotides in each window
# very roughly calculated as # SNPs that pass filter/total number of SNPS * length of window in BP
nsitesperwin <- sapply(slide@region.data@included, FUN=length) # number of SNPs
nincperwin <- sapply(slide@region.data@included, FUN=sum) # number of SNPs flagged as "included" in each window
nucdiv <- nucdiv/(nincperwin/nsitesperwin*winsz)
	head(nucdiv)
	summary(nucdiv)
	
# remove infinite values
nucdiv[is.infinite(nucdiv)] <- NA

# smoothing lines
ids <- 1:nwin
loess.nucdiv1 <- loess(nucdiv[,1] ~ ids, span=0.03)
loess.nucdiv2 <- loess(nucdiv[,2] ~ ids, span=0.03)
loess.nucdiv3 <- loess(nucdiv[,3] ~ ids, span=0.03)

plot(predict(loess.nucdiv1), type = "l", xaxt="n", xlab="position (Mb)", ylab="nucleotide diversity", main = "Chromosome (10kb windows)", ylim=c(0,0.015))

lines(predict(loess.nucdiv2), col="blue")
lines(predict(loess.nucdiv3), col="red")
axis(1,c(1,1000,2000,3000,4000,5000), c("0","10","20","30","40","50"))

# create the legend
legend("topright",c("LF07","LOF11","LOF14"),col=c("black","blue","red"), lty=c(1,1,1))