# Annotate outlier loci and output a table

require(PopGenome)
require(data.table)
require(rtracklayer) # source("https://bioconductor.org/biocLite.R"); biocLite("rtracklayer")

# read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
nullmod <- fread("gzcat analysis/wfs_nullmodel_outliers_07-11-14_Can_25k.tsv.gz") # on mac 108MB file
	outl <- nullmod[outlierLof071114_q3==1 | outlierCan_q3==1 | outlierLof071114_Can_q3==1, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, q3.Lof071114, q3.Can, q3.comb071114Can)] # the outliers
	dim(outl)

# NOT WORKING: read in genome data with readData(). folder has:
# 	All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf
# 	All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf
#	gLof <- readData('../genome_data/VCFLof', format='VCF', gffpath='../genome_data/gff')
#	gCan <- readData('../genome_data/VCFCan', format='VCF', gffpath='../genome_data/gff')
#	
#	get.sum.data(gLof) # seems to have read only 2058 snps
#
#	# set reference. not sure yet what this is used for
#	gLof <- set.synnonsyn(gLof, ref.chr='../genome_data/gadMor2.fasta', save.codons=TRUE) # error: no coding snps
#
#	# look at data
#	length(gLof@region.data@biallelic.sites[[1]]) # location of SNPs
#	head(gLof@region.data@biallelic.sites[[1]]) # location of SNPs
#	length(gLof@region.data@biallelic.matrix)
#	dim(gLof@region.data@biallelic.matrix[[1]]) # size of the matrix of major and minor alleles. Two rows for each individual
#	gLof@region.data@biallelic.matrix[[1]][1:10,1:20] # look at the matrix of major and minor alleles. Two rows for each individual
#
#	get.codons(gLof)
	
# read in genome data with readVCV() (one chr at a time... so annoying)
# need tabix-ed vcf file: on cod node: module load tabix; tabix -p vcf All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz
	g3 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG03', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=29451055, numcols=100000, include.unknown=TRUE) # LG03
	g4 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG04', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=34805322, numcols=100000, include.unknown=TRUE)
	g5 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG05', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=24067848, numcols=100000, include.unknown=TRUE)
	g6 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG06', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25464620, numcols=100000, include.unknown=TRUE)
	g8 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG08', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=26796886, numcols=100000, include.unknown=TRUE)
	g9 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG09', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25382314, numcols=100000, include.unknown=TRUE)
	g10 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG10', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25304306, numcols=100000, include.unknown=TRUE)
	g11 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG11', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=28942968, numcols=100000, include.unknown=TRUE)
	g13 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG13', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25676735, numcols=100000, include.unknown=TRUE)
	g16 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG16', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=31093243, numcols=100000, include.unknown=TRUE)
	g18 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG18', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=22554255, numcols=100000, include.unknown=TRUE)
	g20 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG20', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=24149133, numcols=100000, include.unknown=TRUE)
	g21 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG21', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=22510304, numcols=100000, include.unknown=TRUE)
	g22 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG22', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=21735703, numcols=100000, include.unknown=TRUE)
	g23 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG23', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=23264654, numcols=100000, include.unknown=TRUE)

	#g <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz', tid='LG05', approx=FALSE, gffpath='../genome_data/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=1000000)


	# look at data
	g3@region.names
	get.sum.data(g3) # 68716 biallelic sites
	get.sum.data(g4) # 73859 biallelic sites
	get.sum.data(g5) # 52908 biallelic sites
	get.sum.data(g6) # 60603 biallelic sites
	get.sum.data(g16) # 71846 biallelic sites
	head(g@region.data@biallelic.sites[[1]]) # location of SNPs
	dim(g@region.data@biallelic.matrix[[1]]) # size of the matrix of major and minor alleles. Two rows for each individual
	g@region.data@biallelic.matrix[[1]][1:10,1:20] # look at the matrix of major and minor alleles. Two rows for each individual. not working.

	# define populations
#	g5 <- set.populations(g, list(LOF07=c('BM_209', 'BM_211', 'BM_213', 'BM_214', 'BM_216', 'BM_217', 'BM_218', 'BM_219', 'BM_220', 'BM_221', 'BM_222', 'BM_223', 'BM_224', 'BM_225', 'BM_226', 'BM_227', 'BM_230', 'BM_231', 'BM_232', 'BM_234', 'BM_236', 'BM_237', 'BM_239', 'BM_240', 'BM_241', 'BM_242', 'BM_243', 'BM_244'), LOF11=c('LOF1103001', 'LOF1103002', 'LOF1103003', 'LOF1103004', 'LOF1103005', 'LOF1103006', 'LOF1103007', 'LOF1103008', 'LOF1103009', 'LOF1103010', 'LOF1103011', 'LOF1103012', 'LOF1103013', 'LOF1103014', 'LOF1103015', 'LOF1103016', 'LOF1103017', 'LOF1103018', 'LOF110301', 'LOF1103020', 'LOF1103021', 'LOF1103022', 'LOF1103023', 'LOF1103024'), LOF14=c('LOF1403001', 'LOF1403002', 'LOF1403003', 'LOF1403004', 'LOF1403005', 'LOF1403006', 'LOF1403007', 'LOF1403008', 'LOF1403009', 'LOF1403010', 'LOF1403011', 'LOF1403012', 'LOF1403013', 'LOF1403014', 'LOF1403015', 'LOF1403016', 'LOF1403017', 'LOF1403018', 'LOF1403019', 'LOF1403020', 'LOF1403021', 'LOF1403022', 'LOF1403023', 'LOF1403024')), diploid=TRUE)

	# set reference. used to define codons. 
	# very annoying that each chr must be own file. used scripts/split_fasta.sh to make them.
	g3 <- set.synnonsyn(g3, ref.chr='../genome_data/LG03.fasta', save.codons=TRUE)
	g4 <- set.synnonsyn(g4, ref.chr='../genome_data/LG04.fasta', save.codons=TRUE)
	g5 <- set.synnonsyn(g5, ref.chr='../genome_data/LG05.fasta', save.codons=TRUE)
	g6 <- set.synnonsyn(g6, ref.chr='../genome_data/LG06.fasta', save.codons=TRUE)
	g8 <- set.synnonsyn(g8, ref.chr='../genome_data/LG08.fasta', save.codons=TRUE)
	g9 <- set.synnonsyn(g9, ref.chr='../genome_data/LG09.fasta', save.codons=TRUE)
	g10 <- set.synnonsyn(g10, ref.chr='../genome_data/LG10.fasta', save.codons=TRUE)
	g11 <- set.synnonsyn(g11, ref.chr='../genome_data/LG11.fasta', save.codons=TRUE)
	g13 <- set.synnonsyn(g13, ref.chr='../genome_data/LG13.fasta', save.codons=TRUE)
	g16 <- set.synnonsyn(g16, ref.chr='../genome_data/LG16.fasta', save.codons=TRUE)
	g18 <- set.synnonsyn(g18, ref.chr='../genome_data/LG18.fasta', save.codons=TRUE)
	g20 <- set.synnonsyn(g20, ref.chr='../genome_data/LG20.fasta', save.codons=TRUE)
	g21 <- set.synnonsyn(g21, ref.chr='../genome_data/LG21.fasta', save.codons=TRUE)
	g22 <- set.synnonsyn(g22, ref.chr='../genome_data/LG22.fasta', save.codons=TRUE)
	g23 <- set.synnonsyn(g23, ref.chr='../genome_data/LG23.fasta', save.codons=TRUE)

	g <- list(g3, g4, g5, g6, g8, g9, g10, g11, g13, g16, g18, g20, g21, g22, g23)
	names(g) <- c('LG03', 'LG04', 'LG05', 'LG06', 'LG08', 'LG09', 'LG10', 'LG11', 'LG13', 'LG16', 'LG18', 'LG20', 'LG21', 'LG22', 'LG23')

# get codons and SNP effects
codons <- lapply(g, get.codons, 1) # fails for end of LG16 for some reason

	dim(codons[['LG05']]) # only 3685 on LG05? because only codons
	head(codons[['LG05']])
	
	# print codons for the outliers
	outlcodons <- NULL
	for(chr in sort(unique(outl[,CHROM]))){
		this <- merge(codons[[chr]], outl[CHROM==chr,], by.x='Position', by.y='POS')
		if(nrow(this)>0){
			if(is.null(outlcodons)) outlcodons <- cbind(data.frame(CHROM=chr), this)
			else outlcodons <- rbind(outlcodons, cbind(data.frame(CHROM=chr), this))
		}
	}	
	outlcodons
	

# DOESN'T WORK: get GFF information by SNPs using PopGenome package
# can position be a range?? I don't think so
# doesn't seem to return much useful information
#	anno <- NULL
#	for(i in 1:nrow(outl)){
#		cat(paste(i, ''))
#		this <- get_gff_info(object=FALSE, gff.file='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', chr=outl[i,CHROM], position=outl[i,POS])
#		if(length(this)>0){
#			if(is.null(anno)) anno <- data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS],gff=this)
#			else anno <- rbind(anno, data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS],gff=this))
#		}
#	}
	
# Read in GFF file
gff <- fread("../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff", header=F)
setnames(gff, c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'))

	anno <- NULL
	for(i in 1:nrow(outl)){
		j <- gff[,which(seqname==outl[i,CHROM] & start < outl[i,POS] & end > outl[i,POS] & feature!='contig')]
		if(length(j)>0){
			if(is.null(anno)) anno <- cbind(data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS]), gff[j,])
			else anno <- rbind(anno, cbind(data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS]), gff[j,]))
		} else {
			if(is.null(anno)) anno <- cbind(data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS], seqname='', source='', feature='', start='', end='', score='', strand='', frame='', attribute=''))
			else anno <- rbind(anno, cbind(data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS], seqname='', source='', feature='', start='', end='', score='', strand='', frame='', attribute='')))
		
		}
	}
	anno