#!/bin/bash
# arguments:
# 1: input file
# 2: output file prefix
# 3: population file (list of individuals)
# 4: maximum r2 value for plink trimming

module load vcftools
module load plink

# convert vcf to plink
# trim out inversions, unplaced
# trim to only one population
# trim to relatively high quality loci
vcftools --gzvcf $1 --plink --out $2 --keep $3 --minGQ 15 --minDP 3 --not-chr LG01 --not-chr LG02 --not-chr LG07 --not-chr LG12 --not-chr Unplaced --chrom-map data/chrom_map

# trim to loci with r2 < $4
plink --noweb --file $2 --indep-pairwise 100 10 $4

# recode into a new plink file (.map and .ped)
plink --noweb --file $2 --extract plink.prune.in --recode --out $2

# use PGDspider to convert to genepop
java -Xmx1024m -Xms512m -jar /projects/cees/in_progress/historic_cod_malin/bin/PGDSpider_2.1.1.1/PGDSpider2-cli.jar -inputfile $2.ped -inputformat 'PED' -outputfile $2_r$4.txt -outputformat 'GENEPOP' -spid /projects/cees/in_progress/historic_cod_malin/bin/PGDSpider_2.1.1.1/ped_to_genepop.spid

# how many loci kept?
echo 'Number of loci kept:'
cat plink.prune.in | wc -l

# clean up intermediate plink files
rm plink.nosex
rm plink.prune.in
rm plink.prune.out
