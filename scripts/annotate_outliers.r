# Annotate outlier loci and output a table

require(PopGenome)
require(data.table)
library(RcppCNPy) # for reading pcangsd output

###################################################
# Read in per-locus or per-region results from the various tests
###################################################

nodam2 <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab') # list of loci that pass nodam2 filter
setnames(nodam2, c('CHROM', 'POS', 'REF', 'ALT'))

# read in pcangsd results
pcangsdcan <- as.data.table(npyLoad('analysis/pcangsd_can_gatk_no_dam2_unlink.selection.npy')) # from angsd_pcangsdoutlier.sh
pcangsdlof0711 <- as.data.table(npyLoad('analysis/pcangsd_lof0711_gatk_no_dam2_unlink.selection.npy'))
pcangsdlof0714 <- as.data.table(npyLoad('analysis/pcangsd_lof0714_gatk_no_dam2_unlink.selection.npy'))
pcangsdcan_sites <- fread('analysis/pcangsd_can_gatk_no_dam2_unlink.sites', header = FALSE)
pcangsdlof0711_sites <- fread('analysis/pcangsd_lof0711_gatk_no_dam2_unlink.sites', header = FALSE)
pcangsdlof0714_sites <- fread('analysis/pcangsd_lof0714_gatk_no_dam2_unlink.sites', header = FALSE)

pcangsdcan_sites[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]
pcangsdcan_sites[, POS := as.integer(POS)]
pcangsdlof0711_sites[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]
pcangsdlof0711_sites[, POS := as.integer(POS)]
pcangsdlof0714_sites[, c('CHROM', 'POS') := tstrsplit(V1, '_', fixed = TRUE)]
pcangsdlof0714_sites[, POS := as.integer(POS)]

pcangsdcan <- cbind(pcangsdcan_sites[, .(CHROM, POS)], pcangsdcan)
pcangsdlof0711 <- cbind(pcangsdlof0711_sites[, .(CHROM, POS)], pcangsdlof0711)
pcangsdlof0714 <- cbind(pcangsdlof0714_sites[, .(CHROM, POS)], pcangsdlof0714)

pcangsdcan[, comp := 'Canada']
pcangsdlof0711[, comp := 'Norway 1907-2011']
pcangsdlof0714[, comp := 'Norway 1907-2014']

pcangsd <- rbind(pcangsdcan, pcangsdlof0711, pcangsdlof0714)

pcangsd[, p := pchisq(q = V1, df = 1, lower.tail = FALSE), by = comp]
pcangsd[, testval := p.adjust(p, method = 'fdr'), by = comp]
pcangsd[, testvaltype := 'FDR-corrected p-value']
pcangsd[, region := 0] # a SNP, not a region
pcangsd[, test := 'pcangsd']
pcangsd[, midPos := POS]
pcangsd[, POS := as.character(POS)]

pcangsd[, c('V1', 'p') := NULL]


# wfs null model
wfsLof <- readRDS(file=paste('analysis/wfs_nullmodel_pos&pvals_07-11-14.rds', sep='')) # from wfs_nullmodel_combine.r/wfs_nullmodel_lowp_rerun.r
wfsCan <- readRDS(file=paste('analysis/wfs_nullmodel_pos&pvals_Can.rds', sep=''))

wfsLof[, comp := 'Norway 1907-2011-2014']
wfsCan[, comp := 'Canada']
wfs <- rbind(wfsLof[, .(CHROM, POS, p, comp)], wfsCan[, .(CHROM, POS, p, comp)])
wfs <- merge(wfs, nodam2[, .(CHROM, POS)], by = c('CHROM', 'POS')) # trim to focal loci
wfs[, testval := p.adjust(p), by = comp]
wfs[, p := NULL]
wfs[, testvaltype := 'FDR-corrected p-value']
wfs[, region := 0] # a SNP, not a region
wfs[, test := 'Wright-Fisher null model']
wfs[, midPos := POS]
wfs[, POS := as.character(POS)]


#site-shuffle fst (region)
windsz <- 50000
fstshuff <- fread('output/fst_siteshuffle.angsd.gatk.csv.gz') # from angsd_fst_siteshuffle_null_stats.r
fstshuff <- fstshuff[nloci >1, ]
fstshuff[, POS := paste0(midPos - windsz/2, '-', midPos + windsz/2), by = 1:nrow(fstshuff)]
fstshuff[pop == 'can', comp := 'Canada']
fstshuff[pop == 'lof0711', comp := 'Norway 1907-2011']
fstshuff[pop == 'lof0714', comp := 'Norway 1907-2014']
setnames(fstshuff, 'p', 'testval')
fstshuff[, testvaltype := 'Genome-wide p-value for the region']
fstshuff[, region := 1] # a region, not a SNP
fstshuff[, test := 'FST outlier region']

fstshuff[, c('fst', 'nloci', 'pop') := NULL]


#site-shuffle pi and D (region)
# assumes same windsz as for fst
thetashuff <- fread('analysis/theta_siteshuffle.angsd.gatk.csv.gz') # from angsd_theta_siteshuffle_null_stats.r
thetashuff[, POS := paste0(WinCenter - windsz/2, '-', WinCenter + windsz/2), by = 1:nrow(thetashuff)]
thetashuff[pop == 'can', comp := 'Canada']
thetashuff[pop == 'lof0711', comp := 'Norway 1907-2011']
thetashuff[pop == 'lof0714', comp := 'Norway 1907-2014']
thetashuff[, testvaltype := 'Genome-wide p-value for the region']
thetashuff[, region := 1] # a region, not a SNP
pishuff <- thetashuff[, .(CHROM = Chromo, POS, comp, midPos = WinCenter, testval = tPd.p, testvaltype, region, test = 'Pi outlier region')]
Dshuff <- thetashuff[, .(CHROM = Chromo, POS, comp, midPos = WinCenter, testval = tDd.p, testvaltype, region, test = 'D outlier region')]




##############################################
# Combine the results and choose the outliers
##############################################
outl <- rbind(pcangsd[testval < 0.05, ], wfs[testval < 0.65, ], fstshuff[testval < 0.1, ], pishuff[testval < 0.05,], Dshuff[testval < 0.05,])


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

	head(codons[['LG05']])
	
# print codons for the outliers
outlcodons <- NULL
for(chr in outl[,sort(unique(CHROM))]){
	this <- merge(codons[[chr]], outl[CHROM==chr & region == 0,. (CHROM, midPos)], by.x='Position', by.y='midPos') # trim to outlier
	
	reglns <- outl[CHROM == chr & region == 1, ] # any outlier regions on this chromsomes
	if(nrow(reglns)>0){
	  for(j in 1:nrow(reglns)){
	    onereg <- codons[[chr]][codons[[chr]]$Position >= (reglns$midPos[j] - windsz/2) & codons[[chr]]$Position <= (reglns$midPos[j] + windsz/2), ]
	    if(j == 1) thisreg <- onereg
	    if(j > 1) thisreg <- rbind(thisreg, onereg)
	  }
	}
	
	this <- rbind(this, thisreg)
	
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
  if(outl$region[i] == 0){ # if not a region outlier
    j <- gff[,which(seqname==outl[i,CHROM] & start < outl[i,midPos] & end > outl[i,midPos] & feature!='contig')] # find all entries that overlap this position
    k <- gff[,which(seqname==outl[i,CHROM] & start < (outl[i,midPos]+25000) & end > (outl[i,midPos]-25000) & feature!='contig')] # find all entries within 25kb of this position
    k <- setdiff(k,j)
  } else {
    j <- gff[,which(seqname==outl[i,CHROM] & start < outl[i,midPos + windsz/2] & end > outl[i,midPos - windsz/2] & feature!='contig')] # find all entries that overlap this region
    k <- gff[,which(seqname==outl[i,CHROM] & start < (outl[i,midPos + windsz/2]+25000) & end > (outl[i,midPos - windsz/2]-25000) & feature!='contig')] # find all entries within 25kb of this region
    k <- setdiff(k,j)
  }
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
		minstartindex <- which.min(abs(gff[k,start]-outl[i,midPos]))
		minstart <- abs(gff[k,start]-outl[i,midPos])[minstartindex]
		minendindex <- which.min(abs(gff[k,end]-outl[i,midPos]))
		minend <- abs(gff[k,end]-outl[i,midPos])[minendindex]
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

# remove any duplicate rows
anno <- anno[!duplicated(anno), ]

head(anno)

##########################################	
# Combine SNP effects and annotations
##########################################
# optional: read out outlcodons from analysis/outlier_annotation.csv (if anno only was updated)
# outlcodons <- fread('analysis/outlier_annotation.csv', select=c('CHROM', 'POS', "Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)"))
# setnames(outlcodons, 'POS', 'Position')

if(!is.null(outlcodons)){
	anno2 <- merge(anno, outlcodons, by.y=c('CHROM', 'Position'), by.x=c('CHROM', 'POS'), all.x=TRUE)
} else {
	anno2 <- anno
	anno2[,c("Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", "Polarity (major)", "Polarity (minor)")] <- NA
}
anno2 <- merge(anno2, outl, by = c('CHROM', 'POS'))
anno2 <- as.data.frame(apply(anno2, MARGIN=2, function(x){x[is.na(x)] <- ''; return(x)})) # turn NA to ''
anno2 <- anno2[, c("CHROM", "POS", "midPos", "comp", "test", "testval", "testvaltype", "region", 
                   "Codons (minor)", "Codons (major)", "Protein (minor)", "Protein (major)", "synonymous", 
                   "Polarity (major)", "Polarity (minor)", "feature", "start", "end", "strand", "frame", "ID", 
                   "Names", "Parent", "WithinAnno", "WithinGO", "NearGene", "NearAnno", "NearGO")] # reorder columns
anno2 <- anno2[order(anno2$CHROM, anno2$midPos),]

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

############
# Write out
############

write.csv(anno2, file='tables/outlier_annotation.csv', row.names=FALSE)
	