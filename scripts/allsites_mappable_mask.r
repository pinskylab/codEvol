# create mask file for vcftools
# for running on cod node

require(data.table, lib.loc="/projects/cees/lib/R_packages/") # on cod node

allsites <- fread('zcat all_sites_data_29_06_18/AllSites.kept.sites.kmer25.gz') # all callable and mappable sites

lens <- data.table(CHROM=c('LG01', 'LG02', 'LG03', 'LG04', 'LG05', 'LG06', 'LG07', 'LG08', 'LG09', 'LG10', 'LG11', 'LG12', 'LG13', 'LG14', 'LG15', 'LG16', 'LG17', 'LG18', 'LG19', 'LG20', 'LG21', 'LG22', 'LG23'), LEN=c(28303952, 24054406, 29451055, 34805322, 24074055, 25464620, 31232877, 26796886, 25382314, 25304306, 28942968, 27297974, 25676735, 29296932, 26597959, 31093243, 19149207, 22554255, 21176260, 24149133, 22510304, 21735703, 23264654))

# write out
outfile <- 'all_sites_data_29_06_18/AllSites.kept.sites.kmer25.mask.gz'
gz1 <- gzfile(outfile, 'w') # open file for writing

for(i in 1:nrow(lens)){ # for each chrom
	print(i)
	# write CHROM name
	writeLines(paste('>', lens$CHROM[i], sep=''), con=gz1)

	# create mask vector
	out <- rep(1, lens$LEN[i]) # set to mask everything
	good <- allsites[CHROM==lens$CHROM[i], POS] # the sites that are ok
	out[good] <- 0 # turn the good sites to 0 so not masked out
	
	# write mask vector
	writeLines(paste(out, collapse=''), con=gz1)
}

close(gz1)