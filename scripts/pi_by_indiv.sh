#!/bin/bash

# use export to make variables available to sub-processes
export indvfile="/projects/cees/in_progress/historic_cod_malin/data/all_indivs.txt"
export outfile="/projects/cees/in_progress/historic_cod_malin/analysis/angsd/pi_by_indiv.tsv"
export lg="LG03"
export vcfloc="/projects/cees/in_progress/historic_cod_malin/data/VCF_to_Malin_22.08.16.recode.vcf.gz"
export temploc="/projects/cees/in_progress/historic_cod_malin/analysis/temp/"
export failoc="/projects/cees/databases/gadus_morhua/gadMor2/gadMor2.fasta.fai"
export faloc="/projects/cees/databases/gadus_morhua/gadMor2/gadMor2.fasta"

module load vcftools/0.1.14.zlib.1.2.8
module load angsd/0.912
module load parallel/1.1

function calcpi {
	local indv=$1
	echo $indv

	# output individual to BEAGLE format
	echo "Converting to BEAGLE"
	vcftools --gzvcf $vcfloc --chr $lg --indv $indv --BEAGLE-PL -c | tr ':' '_' | gzip -c > ${temploc}${indv}.BEAGLE.PL.gz

	# calculate Site Allele Frequency from genotype probabilities (folded)
	angsd -beagle ${temploc}${indv}.BEAGLE.PL.gz -fai $failoc -anc $faloc -doSaf 4 -P 4 -fold 1 -intName 1 \
		-out ${temploc}SafFold_$indv

	# maximum likelihood estimate of the SFS
	realSFS ${temploc}SafFold_$indv.saf.idx -P 4 -r $lg > ${temploc}SafFold_$indv.sfs

	# calculate Waterson's and pairwise theta by site
	angsd -beagle ${temploc}${indv}.BEAGLE.PL.gz -fai $failoc -anc $faloc -intName 1 -out ${temploc}SafFold_$indv \
		-doThetas 1 -doSaf 4 -pest ${temploc}SafFold_$indv.sfs -fold 1

	# add indiv name to file contents, including a column header on line 1 (M,N applies sed to only certain lines)
	# double quotes around variable name so that it is evaluated
	zcat ${temploc}SafFold_$indv.thetas.gz | sed -e '1,1s/$/\tIndiv/' -e '2,$s/$/\t'"$indv"'/' | \
		gzip -c > ${temploc}SafFold_$indv.thetas2.gz

}
export -f calcpi # make calcpi available to sub-processes


for indv in $(cat $indvfile)
#for indv in $"Lofoten_BM_219"
do
parallel --no-notice --semaphore -j10 calcpi $indv # use 10 processes at once
done

parallel --no-notice --semaphore --wait # wait until finished

# create outfile header. -e to allow backslash escapes in echo
# /d in sed to remove all header lines (have #Chromo)
echo -e '#Chromo\tPos\tWatterson\tPairwise\tthetaSingleton\tthetaH\tthetaL\tIndiv' > $outfile
zcat ${temploc}SafFold_*.thetas2.gz | sed '/#Chromo/d' >> $outfile
gzip $outfile
