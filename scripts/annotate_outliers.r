# Annotate outlier loci and output a table

require(PopGenome)
require(data.table)
#require(rtracklayer) # source("https://bioconductor.org/biocLite.R"); biocLite("rtracklayer")

##############################
# Read in WFS drift null model outlier values for all populations (combined by wfs_nullmodel_analysis_Canada,1907-2011&2014.r)
##############################
nullmod <- fread("gzcat analysis/wfs_nullmodel_outliers_07-11-14_Can.csv.gz") # on mac 65MB file
outl <- nullmod[q3.Lof0711<0.2 | q3.Lof0714<0.2 | q3.Can<0.2 | q3.comb0711Can<0.2 | q3.comb0714Can<0.2, .(CHROM, POS, Freq_07, Freq_11, Freq_14, Freq_Can40, Freq_CanMod, q3.Lof0711, q3.Lof0714, q3.Can, q3.comb0711Can, q3.comb0714Can)] # trim to just the outliers
	dim(outl) # 35


###################
# Get codon information
###################
outl[,table(CHROM)] # which LGs to load?

	
# read in genome data with readVCV() (one chr at a time... so annoying)
# need tabix-ed vcf file: on cod node: module load tabix; tabix -p vcf All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz
# 	g3 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG03', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=29451055, numcols=100000, include.unknown=TRUE) # LG03
	g4 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG04', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=34805322, numcols=100000, include.unknown=TRUE)
	g5 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG05', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=24067848, numcols=100000, include.unknown=TRUE)
# 	g6 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG06', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25464620, numcols=100000, include.unknown=TRUE)
 	g8 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG08', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=26796886, numcols=100000, include.unknown=TRUE)
	g9 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG09', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25382314, numcols=100000, include.unknown=TRUE)
	g10 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG10', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25304306, numcols=100000, include.unknown=TRUE)
	g11 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG11', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=28942968, numcols=100000, include.unknown=TRUE)
	g13 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG13', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25676735, numcols=100000, include.unknown=TRUE)
 	g14 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG14', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=29296932, numcols=100000, include.unknown=TRUE)
	g15 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG15', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=26597959, numcols=100000, include.unknown=TRUE)
	g16 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG16', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=31093243, numcols=100000, include.unknown=TRUE)
 	g17 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG17', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=19149207, numcols=100000, include.unknown=TRUE)
	g18 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG18', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=22554255, numcols=100000, include.unknown=TRUE)
 	g19 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG19', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=21176260, numcols=100000, include.unknown=TRUE)
# 	g20 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG20', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=24149133, numcols=100000, include.unknown=TRUE)
#	g21 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG21', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=22510304, numcols=100000, include.unknown=TRUE)
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
	head(g3@region.data@biallelic.sites[[1]]) # location of SNPs
	dim(g3@region.data@biallelic.matrix[[1]]) # size of the matrix of major and minor alleles. Two rows for each individual
	g3@region.data@biallelic.matrix[[1]][1:10,1:20] # look at the matrix of major and minor alleles. Two rows for each individual. not working.

	# set reference. used to define codons. 
	# very annoying that each chr must be own file. used scripts/split_fasta.sh to make them.
	if(exists('g3')) g3 <- set.synnonsyn(g3, ref.chr='../genome_data/LG03.fasta', save.codons=TRUE)
	if(exists('g4')) g4 <- set.synnonsyn(g4, ref.chr='../genome_data/LG04.fasta', save.codons=TRUE)
	if(exists('g5')) g5 <- set.synnonsyn(g5, ref.chr='../genome_data/LG05.fasta', save.codons=TRUE)
	if(exists('g6')) g6 <- set.synnonsyn(g6, ref.chr='../genome_data/LG06.fasta', save.codons=TRUE)
	if(exists('g8')) g8 <- set.synnonsyn(g8, ref.chr='../genome_data/LG08.fasta', save.codons=TRUE)
	if(exists('g9')) g9 <- set.synnonsyn(g9, ref.chr='../genome_data/LG09.fasta', save.codons=TRUE)
	if(exists('g10')) g10 <- set.synnonsyn(g10, ref.chr='../genome_data/LG10.fasta', save.codons=TRUE)
	if(exists('g11')) g11 <- set.synnonsyn(g11, ref.chr='../genome_data/LG11.fasta', save.codons=TRUE)
	if(exists('g13')) g13 <- set.synnonsyn(g13, ref.chr='../genome_data/LG13.fasta', save.codons=TRUE)
	if(exists('g14')) g14 <- set.synnonsyn(g14, ref.chr='../genome_data/LG14.fasta', save.codons=TRUE)
	if(exists('g15')) g15 <- set.synnonsyn(g15, ref.chr='../genome_data/LG15.fasta', save.codons=TRUE)
	if(exists('g16')) g16 <- set.synnonsyn(g16, ref.chr='../genome_data/LG16.fasta', save.codons=TRUE)
	if(exists('g17')) g17 <- set.synnonsyn(g17, ref.chr='../genome_data/LG17.fasta', save.codons=TRUE)
	if(exists('g18')) g18 <- set.synnonsyn(g18, ref.chr='../genome_data/LG18.fasta', save.codons=TRUE)
	if(exists('g19')) g19 <- set.synnonsyn(g19, ref.chr='../genome_data/LG19.fasta', save.codons=TRUE)
	if(exists('g20')) g20 <- set.synnonsyn(g20, ref.chr='../genome_data/LG20.fasta', save.codons=TRUE)
	if(exists('g21')) g21 <- set.synnonsyn(g21, ref.chr='../genome_data/LG21.fasta', save.codons=TRUE)
	if(exists('g22')) g22 <- set.synnonsyn(g22, ref.chr='../genome_data/LG22.fasta', save.codons=TRUE)
	if(exists('g23')) g23 <- set.synnonsyn(g23, ref.chr='../genome_data/LG23.fasta', save.codons=TRUE)

	gs <- c('g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'g15', 'g16', 'g17', 'g18', 'g19', 'g20', 'g21', 'g22', 'g23')
	lgs <- c('LG03', 'LG04', 'LG05', 'LG06', 'LG07', 'LG08', 'LG09', 'LG10', 'LG11', 'LG12', 'LG13', 'LG14', 'LG15', 'LG16', 'LG17', 'LG18', 'LG19', 'LG20', 'LG21', 'LG22', 'LG23')
	g <- vector('list', 0)
	ii <- 1
	for(i in 1:length(gs)){
		if(sum(outl$CHROM==lgs[i])>0){
			g[[ii]] <- get(gs[i])
			names(g)[ii] <- lgs[i]
			ii <- ii+1
		}
	}
	
	# g <- list(g4, g5, g8, g9, g10, g11, g13, g14, g15, g16, g17, g18, g19, g22, g23)

# get codons and SNP effects
codons <- lapply(g, get.codons, 1)

	dim(codons[['LG05']]) # only 3685 on LG05? because only codons
	head(codons[['LG05']])
	
	# print codons for the outliers
	outlcodons <- NULL
	for(chr in sort(unique(outl[,CHROM]))){
		this <- merge(codons[[chr]], outl[CHROM==chr,c('CHROM', 'POS')], by.x='Position', by.y='POS')
		if(nrow(this)>0){
			if(is.null(outlcodons)){
				outlcodons <- this
			} else {
				outlcodons <- rbind(outlcodons, this)
			}
		}
	}	
	outlcodons
	

##########################
# Get annotation from GFF
##########################
# Read in GFF file to get gene names and annotations
gff <- fread("../genome_data/gadMor2_annotation_filtered_only_gene_models.gff", header=F)
setnames(gff, c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'))

	# function to search for a string, remove it, and return the rest. Return NA if not found.
	grepstrip <- function(x, strip, na=NA){
		val <- grep(strip, x, value=TRUE)
		if(length(val)>0){
			return(gsub(strip, '', val))
		} else {
			return(na)
		}
	}

# Annotate with gene names and nearby genes
	anno <- NULL
	for(i in 1:nrow(outl)){
		j <- gff[,which(seqname==outl[i,CHROM] & start < outl[i,POS] & end > outl[i,POS] & feature!='contig')] # find all entries that overlap this position
		k <- gff[,which(seqname==outl[i,CHROM] & start < (outl[i,POS]+25000) & end > (outl[i,POS]-25000) & feature!='contig')] # find all entries within 25kb of this position
		k <- setdiff(k,j)
		if(length(j)>0){
			pieces <- strsplit(gff[j,attribute], split=';')
			fts <- paste(gff[j,feature], collapse=';')
			sts <- paste(gff[j,start], collapse=';')
			ends <- paste(gff[j,end], collapse=';')
			strnds <- paste(gff[j,strand], collapse=';')
			frms <- paste(gff[j,frame], collapse=';')
			ids <- paste(sapply(pieces, FUN=grepstrip, 'ID='), collapse=';')
			nms <- paste(sapply(pieces, FUN=grepstrip, 'Name='), collapse=';')
			prnt <- paste(sapply(pieces, FUN=grepstrip, 'Parent='), collapse=';')
			sim <- paste(setdiff(unique(sapply(pieces, FUN=grepstrip, strip='Note=Similar to ', na='')), ''), collapse=';')
			ingo <- paste(setdiff(unique(sapply(pieces, FUN=grepstrip, strip='Ontology_term=', na='')), ''), collapse=',')
		} else {
			fts <- sts <- ends <- strnds <- frms <- ids <- nms <- prnt <- sim <- ingo <- ''
		}
		if(length(k)>0 & length(j)==0){
			# find the closest one
			minstartindex <- which.min(abs(gff[k,start]-outl[i,POS]))
			minstart <- abs(gff[k,start]-outl[i,POS])[minstartindex]
			minendindex <- which.min(abs(gff[k,end]-outl[i,POS]))
			minend <- abs(gff[k,end]-outl[i,POS])[minendindex]
			if(minend<=minstart){
				pieces2 <- strsplit(gff[k,attribute][minendindex], split=';')
				neargene <- paste(setdiff(unique(sapply(pieces2, FUN=grepstrip, strip='ID=', na='')), ''), collapse=';')
				nearanno <- paste(setdiff(unique(sapply(pieces2, FUN=grepstrip, strip='Note=Similar to ', na='')), ''), collapse=';')
				neargo <- paste(setdiff(unique(sapply(pieces2, FUN=grepstrip, strip='Ontology_term=', na='')), ''), collapse=',')
			}
			if(minend>minstart){
				pieces2 <- strsplit(gff[k,attribute][minstartindex], split=';')
				neargene <- paste(setdiff(unique(sapply(pieces2, FUN=grepstrip, strip='ID=', na='')), ''), collapse=';')
				nearanno <- paste(setdiff(unique(sapply(pieces2, FUN=grepstrip, strip='Note=Similar to ', na='')), ''), collapse=';')
				neargo <- paste(setdiff(unique(sapply(pieces2, FUN=grepstrip, strip='Ontology_term=', na='')), ''), collapse=',')
			}
		} else {
			neargene <- nearanno <- neargo <- ''
		}			
					
		if(is.null(anno)){
			anno <- data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS], feature=fts, start=sts, end=ends, strand=strnds, frame=frms, ID=ids, Names=nms, Parent=prnt, WithinAnno=sim, WithinGO=ingo, NearGene=neargene, NearAnno=nearanno, NearGO=neargo)
		} else {
			anno <- rbind(anno, data.frame(CHROM=outl[i,CHROM], POS=outl[i,POS], feature=fts, start=sts, end=ends, strand=strnds, frame=frms, ID=ids, Names=nms, Parent=prnt, WithinAnno=sim, WithinGO=ingo, NearGene=neargene, NearAnno=nearanno, NearGO=neargo))
		}
	}
	
	
	anno
	

##########################################	
# Combine SNP effects and annotations
##########################################
# optional: read out outlcodons from analysis/outlier_annotation.csv (if anno only was updated)
# outlcodons <- fread('analysis/outlier_annotation.csv', select=c('CHROM', 'POS', "Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)"))
# setnames(outlcodons, 'POS', 'Position')

if(nrow(outlcodons)>0){
	anno2 <- merge(anno, outlcodons, by.y=c('CHROM', 'Position'), by.x=c('CHROM', 'POS'), all.x=TRUE)
} else {
	anno2 <- anno
	anno2[,c("Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)")] <- NA
}
anno2 <- merge(anno2, outl)
anno2 <- as.data.frame(apply(anno2, MARGIN=2, function(x){x[is.na(x)] <- ''; return(x)})) # turn NA to ''
anno2 <- anno2[, c("CHROM", "POS", "q3.Lof0711", "q3.Lof0714", "q3.Can", "q3.comb0711Can", "q3.comb0714Can", "Freq_07", "Freq_11", "Freq_14", "Freq_Can40", "Freq_CanMod", "Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)", "feature", "start", "end", "strand", "frame", "ID", "Names", "Parent", "WithinAnno", "WithinGO", "NearGene", "NearAnno", "NearGO")] # reorder columns
anno2$POS <- as.numeric(as.character(anno2$POS)) # convert POS to numeric
anno2 <- anno2[order(anno2$CHROM, anno2$POS),]

anno2


###################
# Label clusters
###################
# find distance among outlier loci
anno2 <- anno2[order(anno2$CHROM, anno2$POS),] # sort by position
anno2$ndist <- NA # nearest neighbor at an earlier position (measure only to left)
for(i in 1:nrow(anno2)){
	j <- which(anno2$CHROM == anno2$CHROM[i]) # other loci on same chromosome
	j <- j[j<i] # remove focal locus and any later loci
	if(length(j)>0) anno2$ndist[i] <- min(abs(anno2$POS[i] - anno2$POS[j]))
} 

# label
indx <- 1
anno2$cluster <- as.numeric(NA)
for(i in 1:nrow(anno2)){
	anno2$cluster[i] <- indx
	if(i < nrow(anno2)){
		if(anno2$ndist[i+1]>=25000 | is.na(anno2$ndist[i+1])) indx <- indx+1
	}
}

anno2[,c('CHROM', 'POS', 'ndist', 'cluster')]

length(unique(anno2$cluster)) # how many clusters?
table(anno2$cluster)
max(table(anno2$cluster)) # largest cluster
aggregate(list(maxdist=anno2$POS), by=list(cluster=anno2$cluster), FUN=function(x){ return(max(x) - min(x))}) # size in bp of each cluster

# reorder columns, drop ndist
anno2 <- anno2[, c("CHROM", "POS", "cluster", "q3.Lof0711", "q3.Lof0714", "q3.Can", "q3.comb0711Can", "q3.comb0714Can", "Freq_07", "Freq_11", "Freq_14", "Freq_Can40", "Freq_CanMod", "Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)", "feature", "start", "end", "strand", "frame", "ID", "Names", "Parent", "WithinAnno", "WithinGO", "NearGene", "NearAnno", "NearGO")]


########################
# Summary statistics
########################
	nrow(anno2)
	length(unique(anno2$CHROM)) # number of LGs
	length(grep('gene', anno2$feature)); length(grep('gene', anno2$feature))/nrow(anno2) # fraction in genes
	sum(anno2$WithinAnno != ''); sum(anno2$WithinAnno != '')/nrow(anno2) # fraction in named genes
	length(grep('CDS', anno2$feature)); length(grep('CDS', anno2$feature))/nrow(anno2) # fraction in coding sequence
	length(grep('three_prime_UTR', anno2$feature)); length(grep('three_prime_UTR', anno2$feature))/nrow(anno2) # fraction in 3'UTR
	length(grep('five_prime_UTR', anno2$feature)); length(grep('five_prime_UTR', anno2$feature))/nrow(anno2) # fraction in 5'UTR
	sum(anno2$synonymous==' TRUE', na.rm=TRUE); sum(anno2$synonymous==' TRUE', na.rm=TRUE)/nrow(anno2) # synonymous snps
	sum(anno2$synonymous=='FALSE', na.rm=TRUE); sum(anno2$synonymous=='FALSE', na.rm=TRUE)/nrow(anno2) # non-synonymous snps
		anno2[anno2$synonymous=='FALSE',]
	
	sum(anno2$NearGene != ''); sum(anno2$NearGene != '')/nrow(anno2) # fraction within 25kb of gene (but not in a gene)
	sum(anno2$NearAnno != ''); sum(anno2$NearAnno != '')/nrow(anno2) # fraction within 25kb of annotated gene (but not in a gene)
	length(unique(c(grep('gene', anno2$feature), which(anno2$NearGene != '')))); length(unique(c(grep('gene', anno2$feature), which(anno2$NearGene != ''))))/nrow(anno2) # fraction in OR within 25kb of genes
	
	length(unique(anno2$cluster)) # how many clusters?
	min(table(anno2$cluster)) # smallest cluster
	max(table(anno2$cluster)) # largest cluster
	aggregate(list(maxdist=anno2$POS), by=list(cluster=anno2$cluster), FUN=function(x){ return(max(x) - min(x))}) # size in bp of each cluster

	sum(as.numeric(as.character(anno2$q3.Lof0711)) < 0.2, na.rm=TRUE) # # outliers in Lof 11
	sum(as.numeric(as.character(anno2$q3.Lof0714)) < 0.2, na.rm=TRUE) # # outliers in Lof 14
	sum(as.numeric(as.character(anno2$q3.Can)) < 0.2, na.rm=TRUE) # # outliers in Can
	sum(as.numeric(as.character(anno2$q3.comb0711Can)) < 0.2, na.rm=TRUE) # # outliers combined
	sum(as.numeric(as.character(anno2$q3.comb0714Can)) < 0.2, na.rm=TRUE) # # outliers combined
	
	# table of GO terms for genes containing outlier SNPs
	uniq <- !duplicated(anno2$Names) & anno2$Names != '' # list of unique genes
	gos <- unlist(strsplit(paste(anno2$WithinGO[uniq], collapse=','), split=','))
	gos2 <- gos[gos !='']
	t(t(sort(table(gos2), decreasing=TRUE)))
	
############
# Write out
############

write.csv(anno2, file='tables/outlier_annotation.csv', row.names=FALSE)