#!/bin/bash

# Job name:
#SBATCH --job-name=Variance
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=3G
#
#if you need long time
#   --partition=long

#SBATCH --ntasks-per-node=1



#Do the first steps once!
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

#Do the first steps once!
# Complex calculations that needs some preparation!
# Calculations run per SITE
# First generate an average coverage for all sites to be checked, called depth_all.idepth
#vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz --depth --out depth_all

# Second create a position file with the coordinates of all the SNPs
#vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz --site-depth --out depth_all
#awk '{print $1"\t"$2}' depth_all.ldepth > all.postions.txt



#for f in $(cat chrom_list);do grep $f all.postions.txt > $f;done


mkdir -p output2
mkdir -p error2

rm output2/var_${1}

#get file with postion
#for f in "$(head positions_25k_mer)";
#do echo "$f" > temp1;


####THEN run multiple times; for f in $(cat chrom_list);do sbatch Calculate_variance2.sh $f;done

###output depth DATA
/cluster/work/users/bastiaas/temp/bin/vcftools --gzvcf Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer.vcf.gz --extract-FORMAT-info DP --positions ${1}\
 --out output2/${1} 2> error2/stderr_${1}.txt


cat ${1} | while read line ; do
echo "$line" > temp2_${1};
## Get specific SNP position
x=$(awk '{print $2}' temp2_${1})
grep $x output2/${1}.DP.FORMAT > temp3_${1}

### Transpose
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
	print str
    }
}' temp3_${1} |awk 'NR > 1 {print $0}' > temp4_${1}

##Exponential
#paste depth_all.idepth  ${1}.idepth |awk 'NR >1 {print 10**($6/$3)}'\
##Non_exponential
#paste depth_all.idepth  ${1}.idepth |awk 'NR >1 {print ($6/$3)}'\

### Calculate depth variance
paste depth_all.idepth  temp4_${1} |awk 'NR >1 {print 10**($4/$3)}'\
 |awk  '{ for(i=1;i<=NF;i++) {total[i]+=$i ; sq[i]+=$i*$i ; } } END\
 { for(i=1;i<=NF;i++) printf "%f ",sq[i]/NR-(total[i]/NR)**2 ; printf "\n" ;}' >> output2/var_${1}


done


rm temp2_${1}
rm temp3_${1}
rm temp4_${1}
#rm ${1}.idepth


paste ${1} output2/var_${1} > output2/${1}_var


Concatenate all output.
Then 
source Get_top_5_percent_var.sh All_variance.tab
This is the list that gets excluded

