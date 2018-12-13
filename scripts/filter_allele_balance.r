# to run on cod node

require(data.table, lib.loc="/projects/cees/lib/R_packages/")

gt <- fread('Allele_balance/Allele_balance.GT.FORMAT') # takes 45 sec or so
ad <- fread('Allele_balance/Allele_balance.AD.FORMAT') # takes 45 sec or so

dim(gt)
dim(ad)

setnames(gt, 3:ncol(gt), paste(names(gt)[3:ncol(gt)], '_gt', sep=''))

setkey(gt, CHROM, POS)
setkey(ad, CHROM, POS)
ad2 <- merge(ad, gt)
	
# mask uncalled or homozygote individuals (set to NA)
indivs <- setdiff(names(ad), c('CHROM', 'POS'))
for(i in indivs){
	ad2[get(paste(i, '_gt', sep='')) %in% c('0/0', './.', '1/1'), eval(i):=NA]
}

	# check that it worked
#	ad2[1:10, .(BM_026_gt, BM_026)]
#	ad2[BM_026_gt=='0/1', .(BM_026_gt, BM_026)]
#
#	ad2[1:10, .(BM_027_gt, BM_027)]
#	ad2[BM_027_gt=='0/1', .(BM_027_gt, BM_027)]

# remove gt columns (to clean up)
for(i in indivs){
	ad2[, eval(paste(i, '_gt', sep='')):=NULL]
}

# extract ref and alt alleles for het individuals (slow: takes 30? min or so)
for(i in indivs){
	cat(paste(i, ' ', sep=''))
	ad2[, eval(paste(i, '_ref', sep='')) := as.numeric(sapply(strsplit(get(i), split=','), '[', 1))]
	ad2[, eval(paste(i, '_alt', sep='')) := as.numeric(sapply(strsplit(get(i), split=','), '[', 2))]
}

# sum ref and alt alleles across het individuals
ad2[, sumRef := rowSums(.SD, na.rm=TRUE), .SDcol = grep('_ref', names(ad2))]
ad2[, sumAlt := rowSums(.SD, na.rm=TRUE), .SDcol = grep('_alt', names(ad2))]


# binomial test for allele balance == 0.5
ad2[,r:=1:.N] # row number
ad2[sumAlt>0 & sumRef>0 & !is.na(sumAlt) & !is.na(sumRef),binomp := binom.test(sumRef, sumRef+sumAlt, p=0.5, alternative='two.sided')$p.value, by=r]

	# ad2[,.(sumAlt, sumRef, binomp)]

# calc FDR
ad2[,binompFDR := p.adjust(binomp, method='fdr')]

# write out
#write.table(ad2[,.(CHROM, POS, sumRef, sumAlt)], file='Allele_balance/Allele_balance.binomp.tsv', row.names=FALSE, sep='/t')

write.table(ad2[,.(CHROM, POS, sumRef, sumAlt, binomp, binompFDR)], file='Allele_balance/Allele_balance.binomp.tsv', row.names=FALSE, sep='\t', quote=FALSE)