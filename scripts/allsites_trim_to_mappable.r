# Trim the list of all callable sites (from GATK) to just those that are also uniquely mappable at kmer 25

require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node


# read in all callable sites information
allsites <- fread('zcat all_sites_data_29_06_18/AllSites.kept.sites.gz') # slow 1.2G gzipped file
setkey(allsites, CHROM, POS)

# read in mappable regions
# has regions that pass the mappability test. Doesnt' include the last bp (POS2)
mappable <- fread('Mappability_tracks/Gadmor_unique_25k-mer_150bp_from_edge.bed')
setnames(mappable, c('CHROM', 'POS1', 'POS2', 'MAP')) 

# trim allsites to only sites in mappable regions
setkey(allsites, CHROM, POS)
allsites[,rem:=0]
	# mark them
for(i in 2:nrow(mappable)){ # this loop is slow. 20 min?
	if(i %% 1000 == 0) print(paste(i, 'of', nrow(mappable)))
	
	if(mappable$CHROM[i] == mappable$CHROM[i-1]){ # if this line and previous are on same CHROM
		remPOS <- mappable$POS2[i-1]:(mappable$POS1[i]-1) # POS to remove. includes end of previous interval, but not start of next
		allsites[.(mappable$CHROM[i], remPOS), rem:=1] # POS to remove. includes end of previous interval, but not start of next
	} else { # if previous line is on a different CHROM
		remPOS1 <- mappable$POS2[i-1]:allsites[mappable$CHROM[i-1], max(POS)] # POS to remove from end of last CHROM
		remPOS2 <- 1:(mappable$POS1[i]-1) # POS to remove from start of this CHROM
		allsites[.(mappable$CHROM[i-1], remPOS1), rem:=1]
		allsites[.(mappable$CHROM[i], remPOS2), rem:=1]
	}
}

allsites[,sum(rem)/]
allsites[,sum(rem)/.N] # will remove 78% of sites! wow

	# remove the marked rows
	dim(allsites)
allsites <- allsites[rem==0,]
	dim(allsites)
allsites[,rem:=NULL]


# write out
outfile <- 'all_sites_data_29_06_18/AllSites.kept.sites.kmer25.gz'
write.table(allsites, file=gzfile(outfile), sep='\t', row.names=FALSE, quote=FALSE)
