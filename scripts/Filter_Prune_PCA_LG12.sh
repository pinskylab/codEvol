#Note give one argument_to_name_resulting_file

#you have to adapt the filtering to what you want!
#you may have to change the plotting characterisitics

module purge
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load BCFtools/1.9-foss-2018b
module load PLINK/1.9b_6.13-x86_64


vcftools --gzvcf ../../GVCF/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz --plink --out $1 \
--chr LG12 --from-bp 1300000 --to-bp 13600000 \
--minGQ 20 \
--minDP 9 \
--chrom-map ../../GVCF/chrom_map

#--exclude-bed ../../data/SNP_regions_with_top5_perc_SNPdensity_150bp_flanking.tab \

#plink --file $1 --indep-pairwise 100 10 0.8


# 100 SNPs 10 step R 0.5

plink --file $1 --make-bed --out $1.pruned

#plink --file $1 --extract plink.prune.in --make-bed --out $1.pruned

#MDS

module purge
module load EIGENSOFT/7.2.1-intel-2019a


cp $1.pruned.bim $1.pruned.pedsnp
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $1.pruned.fam > $1.test
paste $1.test pop_assign.txt > $1.pruned.pedind


echo "genotypename: ${1}.pruned.bed" > $1.par
echo "snpname: ${1}.pruned.pedsnp" >> $1.par
echo "indivname: ${1}.pruned.pedind" >> $1.par
echo "evecoutname: ${1}.pca.evec" >> $1.par
echo "evaloutname: ${1}.eval" >> $1.par
echo "altnormstyle: NO" >> $1.par
echo "lsqproject: YES" >> $1.par
echo "poplistname: pop_project.txt" >> $1.par
echo "snpweightoutname: ${1}.SNP.loadings" >> $1.par
echo "numoutevec: 5" >> $1.par
echo "numoutlieriter: 0" >> $1.par
echo "numoutlierevec: 2" >> $1.par
echo "outliersigmathresh: 6" >> $1.par
echo "qtmode: 0" >> $1.par
echo "fastmode:  NO" >> $1.par

##PCA EIG
smartpca -p $1.par > $1.log





awk '{print $1"\t"$2"\t"$3}' $1.pca.evec > $1.pca.txt

paste pop_assign2.txt $1.pca.txt |awk 'NR>1 {print $0}' > temp


echo  -e "Pop\tInd\tC1\tC2" > tem

cat tem temp > $1.csv
rm tem temp


#plink --bfile $1.pruned --cluster --noweb --mds-plot 2

#plink --file $1 --cluster --noweb --mds-plot 2


#awk '{print $1"\t"$4"\t"$5}' plink.mds > $1.pruned.mds




#paste /projects/researchers/researchers01/bastiaas/aDNA_james_cod/Data_management/pop_plot_assigment_course.mds $1.pruned.mds |awk 'NR<2 {print $0;next}{print $0| "sort -k1"}' > $1.csv

module purge
module load R/3.6.2-foss-2019b

#PLOT



Rscript Plot_mds_less_than_12_groups.r $1.csv


#rm *prune*
#rm plink*
