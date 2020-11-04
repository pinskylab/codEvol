# Annotate loci in outlier regions and output a table

require(PopGenome)
require(data.table)
#require(rtracklayer) # source("https://bioconductor.org/biocLite.R"); biocLite("rtracklayer")

##############################
# Read in data
##############################
# outlier regions
reg <- fread('analysis/outlier_50kregions_shared_07-11-14_Can.csv.gz'); window <- 50000
	reg

reg[,region:=1:nrow(reg)]

# read in SNP positions (gatk nodam2 loci)
gatk <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab', col.names = c('CHROM', 'POS', 'REF', 'ALT'))
	
#############
# get SNPs in each outlier region
#############
outl <- gatk[CHROM==reg$chr[1] & POS>=(reg$WinCenter[1] - window/2) & POS<=(reg$WinCenter[1] + window/2), .(CHROM, POS, region=reg$region[1])]
if(nrow(reg)>1){
	for(i in 2:nrow(reg)){
		temp <- gatk[CHROM==reg$chr[i] & POS>=(reg$WinCenter[i] - window/2) & POS<=(reg$WinCenter[i] + window/2), .(CHROM, POS, region=reg$region[i])]
		outl <- rbind(outl, temp)
	}
}

	nsnps <- outl[, .(nsnps = .N), by=.(CHROM, region)]
	setkey(nsnps, nsnps); nsnps
	
###################
# Get codon information
###################
lgstoload <- outl[,sort(unique(CHROM))] # which LGs to load?


# read in genome data with readVCV() (one chr at a time... so annoying)
# need tabix-ed vcf file: on cod node: module load tabix; tabix -p vcf All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Canada.vcf.gz
if('LG01' %in% lgstoload) error('need to write lg01!')
if('LG02' %in% lgstoload) g2 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG02', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=24054406, numcols=100000, include.unknown=TRUE)
if('LG03' %in% lgstoload) g3 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG03', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=29451055, numcols=100000, include.unknown=TRUE)
if('LG04' %in% lgstoload) g4 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG04', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=34805322, numcols=100000, include.unknown=TRUE)
if('LG05' %in% lgstoload) g5 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG05', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=24067848, numcols=100000, include.unknown=TRUE)
if('LG06' %in% lgstoload) g6 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG06', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25464620, numcols=100000, include.unknown=TRUE)
if('LG07' %in% lgstoload) error('need to write lg07!')
if('LG08' %in% lgstoload) g8 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG08', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=26796886, numcols=100000, include.unknown=TRUE)
if('LG09' %in% lgstoload) g9 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG09', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25382314, numcols=100000, include.unknown=TRUE)
if('LG10' %in% lgstoload) g10 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG10', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25304306, numcols=100000, include.unknown=TRUE)
if('LG11' %in% lgstoload) g11 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG11', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=28942968, numcols=100000, include.unknown=TRUE)
if('LG12' %in% lgstoload) error('need to write lg12!')
if('LG13' %in% lgstoload) g13 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG13', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=25676735, numcols=100000, include.unknown=TRUE)
if('LG14' %in% lgstoload) g14 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG14', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=29296932, numcols=100000, include.unknown=TRUE)
if('LG15' %in% lgstoload) g15 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG15', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=26597959, numcols=100000, include.unknown=TRUE)
if('LG16' %in% lgstoload) g16 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG16', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=31093243, numcols=100000, include.unknown=TRUE)
if('LG17' %in% lgstoload) g17 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG17', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=19149207, numcols=100000, include.unknown=TRUE)
if('LG18' %in% lgstoload) g18 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG18', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=22554255, numcols=100000, include.unknown=TRUE)
if('LG19' %in% lgstoload) g19 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG19', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=21176260, numcols=100000, include.unknown=TRUE)
if('LG20' %in% lgstoload) g20 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG20', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=24149133, numcols=100000, include.unknown=TRUE)
if('LG21' %in% lgstoload) g21 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG21', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=22510304, numcols=100000, include.unknown=TRUE)
if('LG22' %in% lgstoload) g22 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG22', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=21735703, numcols=100000, include.unknown=TRUE)
if('LG23' %in% lgstoload) g23 <- readVCF('../genome_data/All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz', tid='LG23', approx=FALSE, gffpath='../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff', frompos=1, topos=23264654, numcols=100000, include.unknown=TRUE)

# set reference. used to define codons. 
# very annoying that each chr must be own file. used scripts/split_fasta.sh to make them.
g <- vector('list', 0); ii <- 1
if(exists('g1')){ g[[ii]] <- set.synnonsyn(g1, ref.chr='../genome_data/LG01.fasta', save.codons=TRUE); names(g)[ii] <- 'LG01'; ii <- ii+1}
if(exists('g2')){ g[[ii]] <- set.synnonsyn(g2, ref.chr='../genome_data/LG02.fasta', save.codons=TRUE); names(g)[ii] <- 'LG02'; ii <- ii+1}
if(exists('g3')){ g[[ii]] <- set.synnonsyn(g3, ref.chr='../genome_data/LG03.fasta', save.codons=TRUE); names(g)[ii] <- 'LG03'; ii <- ii+1}
if(exists('g4')){ g[[ii]] <- set.synnonsyn(g4, ref.chr='../genome_data/LG04.fasta', save.codons=TRUE); names(g)[ii] <- 'LG04'; ii <- ii+1}
if(exists('g5')){ g[[ii]] <- set.synnonsyn(g5, ref.chr='../genome_data/LG05.fasta', save.codons=TRUE); names(g)[ii] <- 'LG05'; ii <- ii+1}
if(exists('g6')){ g[[ii]] <- set.synnonsyn(g6, ref.chr='../genome_data/LG06.fasta', save.codons=TRUE); names(g)[ii] <- 'LG06'; ii <- ii+1}
if(exists('g7')){ g[[ii]] <- set.synnonsyn(g7, ref.chr='../genome_data/LG07.fasta', save.codons=TRUE); names(g)[ii] <- 'LG07'; ii <- ii+1}
if(exists('g8')){ g[[ii]] <- set.synnonsyn(g8, ref.chr='../genome_data/LG08.fasta', save.codons=TRUE); names(g)[ii] <- 'LG08'; ii <- ii+1}
if(exists('g9')){ g[[ii]] <- set.synnonsyn(g9, ref.chr='../genome_data/LG09.fasta', save.codons=TRUE); names(g)[ii] <- 'LG09'; ii <- ii+1}
if(exists('g10')){ g[[ii]] <- set.synnonsyn(g10, ref.chr='../genome_data/LG10.fasta', save.codons=TRUE); names(g)[ii] <- 'LG10'; ii <- ii+1}
if(exists('g11')){ g[[ii]] <- set.synnonsyn(g11, ref.chr='../genome_data/LG11.fasta', save.codons=TRUE); names(g)[ii] <- 'LG11'; ii <- ii+1}
if(exists('g12')){ g[[ii]] <- set.synnonsyn(g12, ref.chr='../genome_data/LG12.fasta', save.codons=TRUE); names(g)[ii] <- 'LG12'; ii <- ii+1}
if(exists('g13')){ g[[ii]] <- set.synnonsyn(g13, ref.chr='../genome_data/LG13.fasta', save.codons=TRUE); names(g)[ii] <- 'LG13'; ii <- ii+1}
if(exists('g14')){ g[[ii]] <- set.synnonsyn(g14, ref.chr='../genome_data/LG14.fasta', save.codons=TRUE); names(g)[ii] <- 'LG14'; ii <- ii+1}
if(exists('g15')){ g[[ii]] <- set.synnonsyn(g15, ref.chr='../genome_data/LG15.fasta', save.codons=TRUE); names(g)[ii] <- 'LG15'; ii <- ii+1}
if(exists('g16')){ g[[ii]] <- set.synnonsyn(g16, ref.chr='../genome_data/LG16.fasta', save.codons=TRUE); names(g)[ii] <- 'LG16'; ii <- ii+1}
if(exists('g17')){ g[[ii]] <- set.synnonsyn(g17, ref.chr='../genome_data/LG17.fasta', save.codons=TRUE); names(g)[ii] <- 'LG17'; ii <- ii+1}
if(exists('g18')){ g[[ii]] <- set.synnonsyn(g18, ref.chr='../genome_data/LG18.fasta', save.codons=TRUE); names(g)[ii] <- 'LG18'; ii <- ii+1}
if(exists('g19')){ g[[ii]] <- set.synnonsyn(g19, ref.chr='../genome_data/LG19.fasta', save.codons=TRUE); names(g)[ii] <- 'LG19'; ii <- ii+1}
if(exists('g20')){ g[[ii]] <- set.synnonsyn(g20, ref.chr='../genome_data/LG20.fasta', save.codons=TRUE); names(g)[ii] <- 'LG20'; ii <- ii+1}
if(exists('g21')){ g[[ii]] <- set.synnonsyn(g21, ref.chr='../genome_data/LG21.fasta', save.codons=TRUE); names(g)[ii] <- 'LG21'; ii <- ii+1}
if(exists('g22')){ g[[ii]] <- set.synnonsyn(g22, ref.chr='../genome_data/LG22.fasta', save.codons=TRUE); names(g)[ii] <- 'LG22'; ii <- ii+1}
if(exists('g23')){ g[[ii]] <- set.synnonsyn(g23, ref.chr='../genome_data/LG23.fasta', save.codons=TRUE); names(g)[ii] <- 'LG23'; ii <- ii+1}


# get codons and SNP effects
codons <- lapply(g, get.codons, 1)

head(codons[['LG10']])

# print codons for the outliers
## STOPPED UPDATING BELOW HERE ON 11/4/2020
## WAS CONVERTING TO USE SHARED FST OUTLIER REGIONS
## WILL INSTEAD INTEGRATE INTO annotate_outliers.r
outlcodons <- NULL
for(chr in outl[,sort(unique(CHROM))]){
  this <- merge(codons[[chr]], outl[CHROM==chr, .(CHROM, midPos)], by.x='Position', by.y='midPos') # trim to outlier SNPs with codons
  
  reglns <- outl[CHROM == chr & region == 1, ] # any outlier regions on this chromsomes
  if(nrow(reglns)>0){
    for(j in 1:nrow(reglns)){
      onereg <- codons[[chr]][codons[[chr]]$Position >= (reglns$midPos[j] - windsz/2) & codons[[chr]]$Position <= (reglns$midPos[j] + windsz/2), ]
      if(j == 1) thisreg <- onereg
      if(j > 1) thisreg <- rbind(thisreg, onereg)
    }
  }
  
  if(nrow(this) > 0 & nrow(reglns) > 0) this <- rbind(this, thisreg)
  if(nrow(this) == 0 & nrow(reglns) > 0) this <- thisreg
  
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

if(!is.null(outlcodons)){
	if(nrow(outlcodons)>0){
		anno2 <- merge(anno, outlcodons, by.y=c('CHROM', 'Position'), by.x=c('CHROM', 'POS'), all.x=TRUE)
	} else {
		anno2 <- anno
		anno2[,c("Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)")] <- NA	
	}
} else {
	anno2 <- anno
	anno2[,c("Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)")] <- NA
}
anno2 <- merge(anno2, outl)
anno2 <- as.data.frame(apply(anno2, MARGIN=2, function(x){x[is.na(x)] <- ''; return(x)})) # turn NA to ''
anno2 <- anno2[, c("CHROM", "POS", "region", "Freq_07", "Freq_11", "Freq_14", "Freq_Can40", "Freq_CanMod", "Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)", "feature", "start", "end", "strand", "frame", "ID", "Names", "Parent", "WithinAnno", "WithinGO", "NearGene", "NearAnno", "NearGO")] # reorder columns (no p-values)
# anno2 <- anno2[, c("CHROM", "POS", "region", "Freq_07", "Freq_11", "Freq_14", "Freq_Can40", "Freq_CanMod", "pLof0711", "pLof0714", "pCan", "p.comb0711Can", "p.comb0714Can", "Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)", "feature", "start", "end", "strand", "frame", "ID", "Names", "Parent", "WithinAnno", "WithinGO", "NearGene", "NearAnno", "NearGO")] # reorder columns (with p-values)
anno2$POS <- as.numeric(as.character(anno2$POS)) # convert POS to numeric
anno2 <- anno2[order(anno2$CHROM, anno2$POS),]

anno2


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
		
	# table of GO terms for genes containing outlier SNPs
	uniq <- !duplicated(anno2$Names) & anno2$Names != '' # list of unique genes
	gos <- unlist(strsplit(paste(anno2$WithinGO[uniq], collapse=','), split=','))
	gos2 <- gos[gos !='']
	t(t(sort(table(gos2), decreasing=TRUE)))

	# table of GO terms for genes close to outlier SNPs
	uniq <- !duplicated(anno2$NearAnno) & anno2$NearAnno != '' # list of unique genes
	gos <- unlist(strsplit(paste(anno2$NearGO[uniq], collapse=','), split=','))
	gos2 <- gos[gos !='']
	t(t(sort(table(gos2), decreasing=TRUE)))
	
############
# Write out
############

write.csv(anno2, file='tables/outlierregion_freqshared_annotation.csv', row.names=FALSE)