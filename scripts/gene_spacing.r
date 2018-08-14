# histogram of distances among genes

require(data.table)
gff <- fread("../genome_data/gff/gadMor2_annotation_filtered_only_gene_models.gff", header=F)
setnames(gff, c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'))
chroms <- c("LG01", "LG02", "LG03", "LG04", "LG05", "LG06", "LG07", "LG08", "LG09", "LG10", "LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", "LG21", "LG22", "LG23")

# calculate and plot distances between genes
starts <- gff[feature=='gene' & seqname %in% chroms, .(seqname, start)]
ends <- gff[feature=='gene' & seqname %in% chroms, .(seqname, end)]
dists <- starts[2:length(starts)] - ends[1:(length(ends)-1)] # produces negatives... why? on different chromosomes?
	summary(dists)

dists2 <- rep(NA, nrow(ends) - 1)
for(i in 1:nrow(ends)){
	if(i %% 1000 == 0) print(i)
	thesestarts <- starts[seqname == ends[i,seqname] & start > ends[i,end], start]
	if(length(thesestarts)>0){
		dists2[i] <- min(thesestarts - ends[i,end])
	} else {
		dists2[i] <- NA
	}
}


summary(dists2)
hist(dists2, col='grey', breaks=20)


# calculate and plot gene lengths
lens <- gff[feature=='gene' & seqname %in% chroms, end - start]

summary(lens)
hist(lens)