# codEvol

Scripts for analysis of Atlantic cod genome sequencies through time (Northeast Arctic and Canada)

Paleomix_settings_historic_cod.yaml
Settings for handling raw historic reads in Paleomix (Adapter removal, Collapse of Paired End reads, BWA, and Duplicate Removal)

Paleomix_settings_modern_cod.yaml
Settings for handling raw modern reads in Paleomix (Adapter removal, Collapse of Paired End reads, BWA, and Duplicate Removal)

allele_balance_calc.R
Script to calculate bias in the number of allelic variants in heterozygote genotypes. I.e. is there consistent bias for specific variants (expection is 50/50)


allele_balance_plot.R
Script to plot histograms of read-depth of alternate alleles for heterozygote genotypes

filter_allele_balance.r
Script that performs a binomial test on heterozygote genotype calls

filter_allele_balance.saga.sh
Wrapper script to submit the filter_allele_balance.r to the SAGA compute cluster

Resample_ANGSD_slurm
Random resampling of individuals. We resample to get bootstrap intervals of windowed Fst by randomly resampling the temporal comparisons.

ANGSD_snp_calling
SNP calling using ANSGS. We following N. Therkildsen approaches: N. O. Therkildsen, et al., Contrasting genomic shifts underlie parallel phenotypic evolution in response to fishing. Science 365, 487–490 (2019).

Make_eigen_R_PCA_ANGSD.rscript
Script to rescale eigen values coming from PCAANGSD. Ignore R plotting code

PCAANGSD.sh
Running PCA ANGSD from beagle files

ANGSD_pop_freq.sh
Calculating population specific frequencies 

ANGSD_admixture
Running ANSGD admixture for different values of K

Plot_admixture.R
Basic R script to plot admixture output

Filter_prune_admixture
Script to run ADMIXTURE for a number of K on SAGA using GATK generated VCF

Filter_Prune_PCA_LG12.sh
Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status.

Filter_Prune_PCA_LG07.sh
Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status.

Filter_Prune_PCA_LG02.sh
Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status.

Filter_Prune_PCA_LG01.sh
Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status.

Filter_Prune_PCA.sh
Script to run PCA using a GATK generated VCF. Various filters are applied (E.g. a filter for LD, and removal of inversion chromosomes)

Run_bionomical_test.saga.sh
Wrapper script to run the bionomical test on SAGA

Calculate_variance2.sh
Script that calculates a measure of read-depth variation between individuals on a specfic site.

Filter_vcf_files.sh
Filter settings for GATK  VCF file as described. Please note that the script cannot be run in a single iteration: various subscripts (such as the binomial test or variance calculations) have to be run at intermediated stages. These generated lists of SNPs to be excluded. 

Genotype_Caller.sh
Running GATKs genotype caller. Note we run this per chromosome in parallel 

Haplotype_caller_GATK.sh
Running GATKs haplotype caller 

Remove_clipped_reads.sh
Script to remove clipped reads from BAM file. Any clipped reads (including entire pairs for which one reads was clipped) will has been removed

Create_mappability_track.sh.
We use the GEM mapper https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377 to create mappability tracks that can be uniquely mapped. K-mer size determines the size of the unique regions. 
