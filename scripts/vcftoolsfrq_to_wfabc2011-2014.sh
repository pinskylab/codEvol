#!/bin/bash

# Format data for WFABC

vcf='/projects/cees/in_progress/historic_cod_malin/data/VCF_to_Malin_22.08.16.recode.vcf.gz' # VCF file
pop1='/projects/cees/in_progress/historic_cod_malin/data/pop_LOF_S_11.txt' # file to list indivs in pop1
pop2='/projects/cees/in_progress/historic_cod_malin/data/pop_LOF_S_14.txt'
pre1='LOF_S_11_LG03'
pre2='LOF_S_14_LG03'
tempdir='/projects/cees/in_progress/historic_cod_malin/analysis/temp/' # temp directory
outdir='/projects/cees/in_progress/historic_cod_malin/analysis/' # output directory
gen=20 # number of generations
ntimes=2 # number of time points

f1=$tempdir$pre1'.frq' # first file of frequencies
f2=$tempdir$pre2'.frq'
#outfile=$outdir$pre1'_to_'$pre2'.wfabc' # output file
outfile=$outdir$pre1'_to_'$pre2'_notrim.wfabc' # output file

# load programs
module load vcftools/0.1.14.zlib.1.2.8
module load bcftools/1.3
module load samtools/1.3.1 # for bgzip

# extract focal populations from VCF
# only include loci with minor allele count > 0
# could add --minGQ 20 to trim on genotype quality
# vcftools --gzvcf $vcf --chr LG03 --keep $pop1 --mac 2 --recode --stdout | bgzip > $tempdir$pre1'.recode.vcf.gz'
# bcftools index $tempdir$pre1'.recode.vcf.gz' # trim to minor allele count 2
# vcftools --gzvcf $vcf --chr LG03 --keep $pop2 --mac 2 --recode --stdout | bgzip > $tempdir$pre2'.recode.vcf.gz'
# bcftools index $tempdir$pre2'.recode.vcf.gz'
vcftools --gzvcf $vcf --chr LG03 --keep $pop1 --recode --stdout | bgzip > $tempdir$pre1'.recode.vcf.gz'
bcftools index $tempdir$pre1'.recode.vcf.gz' # no trimming
vcftools --gzvcf $vcf --chr LG03 --keep $pop2 --recode --stdout | bgzip > $tempdir$pre2'.recode.vcf.gz'
bcftools index $tempdir$pre2'.recode.vcf.gz'

# intersect the two files and trim each to the same set of loci
bcftools isec $tempdir$pre1'.recode.vcf.gz' $tempdir$pre2'.recode.vcf.gz' -p $tempdir # output files 0002 and 0003 are only the matching loci
bgzip $tempdir'0002.vcf'
bgzip $tempdir'0003.vcf'
bcftools index $tempdir'0002.vcf.gz'
bcftools index $tempdir'0003.vcf.gz'
bcftools merge $tempdir'0002.vcf.gz' $tempdir'0003.vcf.gz' > $tempdir$pre1'_'$pre2'.recode.vcf'

# trim to minor allele freq > 0.025
# commented out: do this in WFABC instead
#vcftools --vcf $tempdir$pre1'_'$pre2'.recode.vcf' --maf 0.025 --recode --out $tempdir$pre1'_'$pre2'_maf0.025'

# extract allele frequencies from VCF
# vcftools --vcf $tempdir$pre1'_'$pre2'_maf0.025.recode.vcf' --keep $pop1 --freq --out $tempdir$pre1
# vcftools --vcf $tempdir$pre1'_'$pre2'_maf0.025.recode.vcf' --keep $pop2 --freq --out $tempdir$pre2
vcftools --vcf $tempdir$pre1'_'$pre2'.recode.vcf' --keep $pop1 --freq --out $tempdir$pre1
vcftools --vcf $tempdir$pre1'_'$pre2'.recode.vcf' --keep $pop2 --freq --out $tempdir$pre2

# OPTIONAL for testing: trim to just the head
#head $f1 > $f1'head'
#head $f2 > $f2'head'
#f1=$tempdir$pre1'.frqhead' # first file of frequencies
#f2=$tempdir$pre2'.frqhead'



# figure out number of sites
nsites=`cat $f1 | tail -n+2 | wc -l` # skip the first line

# extract sample sizes
# could probably use awk to avoid the temporary files
cat $f1 | tail -n+2 | cut -f 4 -d$'\t' > $tempdir'sampsize1'
cat $f2 | tail -n+2 | cut -f 4 -d$'\t' > $tempdir'sampsize2'

# extract frequency of the first allele
cat $f1 | tail -n+2 | cut -f 5 -d$'\t' | cut -f 2 -d":" > $tempdir'alfreq1'
cat $f2 | tail -n+2 | cut -f 5 -d$'\t' | cut -f 2 -d":" > $tempdir'alfreq2'

# multiply sample size by allele freq and round to 0 precision
paste $tempdir'sampsize1' $tempdir'alfreq1' | awk '{ printf "%.0f\n", $1 * $2}' > $tempdir'alcnt1'
paste $tempdir'sampsize2' $tempdir'alfreq2' | awk '{ printf "%.0f\n", $1 * $2}' > $tempdir'alcnt2'

# interleave the sample sizes and the allele freqs
paste -d "," $tempdir'sampsize1' $tempdir'sampsize2' > $tempdir'sampsizes'
paste -d "," $tempdir'alcnt1' $tempdir'alcnt2' > $tempdir'alcnts'

# write header
echo "$nsites $ntimes" > $outfile
echo "0 $gen" >> $outfile

# write samples sizes and allele freqs
paste -d '\n' $tempdir'sampsizes' $tempdir'alcnts' >> $outfile

# a start at awk code, but doesn't work
#awk 'FNR>1 {print "1: "$0; if(getline < "$f2") print "2: "$0}' $f1


# clean up temp files
rm $tempdir$pre1'_'$pre2'.recode.vcf' # the 1st intermediate VCF file
#rm $tempdir$pre1'_'$pre2'_maf0.025.recode.vcf' # the 2nd intermediate VCF file
#rm $tempdir$pre1'_'$pre2'.log' # the 1st log
#rm $tempdir$pre1'_'$pre2'_maf0.025.log'
rm $tempdir$pre1'.recode.vcf.gz'
rm $tempdir$pre2'.recode.vcf.gz'
rm $tempdir$pre1'.recode.vcf.gz.csi'
rm $tempdir$pre2'.recode.vcf.gz.csi'
rm $f1
rm $f2
rm $tempdir$pre1'.log' # log file from vcftools output
rm $tempdir$pre2'.log' 
rm $tempdir'sampsize1'
rm $tempdir'sampsize2'
rm $tempdir'alfreq1'
rm $tempdir'alfreq2'
rm $tempdir'alcnt1'
rm $tempdir'alcnt2'
rm $tempdir'sampsizes'
rm $tempdir'alcnts'
rm $tempdir'README.txt' # output by bcftools isec
rm $tempdir'sites.txt' # output by bcftools isec
rm $tempdir'0000.vcf' # output by bcftools isec
rm $tempdir'0001.vcf' # output by bcftools isec
rm $tempdir'0002.vcf.gz' # output by bcftools isec
rm $tempdir'0003.vcf.gz' # output by bcftools isec
rm $tempdir'0002.vcf.gz.csi' # output by bcftools isec
rm $tempdir'0003.vcf.gz.csi' # output by bcftools isec
