# codEvol

Scripts, supporting data, and figures from analysis of Atlantic cod genome sequences through time (Northeast Arctic and Canada), as published in

Pinsky, ML, AM Eikeset, C Helmerson, IR Bradbury, P Bentzen, C Morris, A Gondek, HT Baalsrud, MSO Brieuc, OS Kjesbu, JA Godiksen, JMI Barth, M Matschiner, NC Stenseth, KS Jakobsen, SJentoft, and B Star. Genomic stability through time despite decades of exploitation in cod on both sides of the Atlantic. *PNAS*.

These scripts contain hard links and specific environment settings. They would require adaptation to run elsewhere and on other data. The scripts were run on a Linux scientific workstation (cod) and two Linux clusters with SLURM job management (abel and saga) at the University of Oslo, plus on a MacBook Pro running R 3.*.

# data/
1. `kjesbu_fecundity/AFWG_Table3.11proportionmatureatage.csv`: Proportion of Northeast Arctic cod mature at each age. From Table 3.11 in [ICES AFWG Report 2016. ICES CM 2016/ACOM:06. Report of the Arctic Fisheries Working Group (AFWG).](http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2016/AFWG/01%20AFWG%20Report%202016.pdf)
1. `kjesbu_fecundity/AFWG_Table3.16stocknumbersatage.csv`: Numbers at age from the stock assessment. Table 3.16 in ICES AFWG Report 2016.
1. `kjesbu_fecundity/AFWG_Table3.6catchnumbersatage.csv`: Numbers at age from catch data. Table 3.16 in ICES AFWG Report 2016.
1 `kjesbu_fecundity/Fekunditet Andenes 1986_2006 Excel_2.csv`: Fecundity data from a local fishing port in northern Norway covering a period of 20 years. Use Main_series == Y (Yes), Oto_Age, whole body weight (in g) in Fish_Weight, liver weight (in  g) in Fish_Liver, Oto_type == 3, 4 or 5 (meaning that we exclude Coastal cod type 1 and 2), fecundity (in mill.) in Fec_Mill. Mean oocyte diameter is in Ooc_Mean_Dia in micrometers; females with small oocytes tend to show a decrease in fecundity up to spawning due to atresia (reabsorption of developing oocytes).  
1. `kjesbu_fecundity/Leverdata.csv`: Relative liver size, called the hepatosomatic index, from [Kjesbu et al. 2014 ICES Journal of Marine Science](https://doi.org/10.1093/icesjms/fsu030)
1. `kjesbu_fecundity/Rollesfsen1953_Fig1.csv`: Female age distribution from Fig. 1 in [Rollefsen 1953 Observations on the cod and cod fisheries of Lofoten. Rapports et Proces-Verbaux des Reunions, Conseil lnternationale pour l'Exploration de la Mer 136: 40-47](http://hdl.handle.net/11250/101129)
1. `phenotypes/AFWG_2019_3_Northeast Arctic_Cod_Table3.18.csv`: Average fishing mortality (F) for ages 5-10 on the Northeast Arctic stock, from the stock assessment. From Table 3.18 in [ICES AFWG Report 2019](http://ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/Fisheries%20Resources%20Steering%20Group/2019/AFWG/AFWG%202019_3_Northeast%20Arctic%20Cod.pdf)
1. `phenotypes/AFWG_Table 3.11.csv`: Proportion mature at each age, Northeast Arctic cod. From Table. 3.11 in [ ICES. 2020. Arctic Fisheries Working Group (AFWG). ICES Scientific Reports 2(52), 577 pp](http://doi.org/10.17895/ices.pub.6050)
1. `phenotypes/Brattey et al_2018_Table18_proportion_mature.csv`: Proportion mature at each age, Canadian population. From Table 18 in [Brattey, J., Cadigan, N., Dwyer, K. S., Healey, B. P., Ings, D. W., Lee, E. M., Maddock Parsons, D., Morgan, M. J., Regular, P., Rideout, R. M. 2018. Assessment of the Northern Cod (Gadus morhua) stock in NAFO Divisions 2J3KL in 2016. DFO Can. Sci. Advis. Sec. Res. Doc. 2018/018. v + 107p](https://waves-vagues.dfo-mpo.gc.ca/Library/40702017.pdf)
1. `phenotypes/Brattey_etal_2018_CSAS_TableA2-5_NorthernCod_fishingmortality.csv`: Fishing mortality (F) by age, from Brattey et al. 2018.
1. `phenotypes/Eikeset_Age_and_length_at_maturation.csv`: Age and length at maturity in Northeast Arctic cod, from [Eikeset et al. 2016 PNAS](https://doi.org/10.1073/pnas.1525749113)
1. `inversions.csv`: Genome locations of the cod inversions. From [Star et al. 2017 PNAS](https://doi.org/10.1073/pnas.1710186114)
1. `lg_length.csv`: Length in bp of each linkage group. From [Tørresen et al 2019 BMC Genomics](https://doi.org/10.1186/s12864-016-3448-x)


# scripts/
## Phenotypes
1. `phenotypes.r`: Calculates the age at 50% mature in Canadian and Northeast Arctic populations. Needs data in `data/phenotypes/`. Writes `output/age_50percmature.csv`.

## Bioinformatics
### Aligning and handling sequence reads:

1. `Paleomix_settings_historic_cod.yaml`: Settings for handling raw historic reads in Paleomix (Adapter removal, Collapse of Paired End reads, BWA, and Duplicate Removal). Requires (paired) fasta.gz and reference genome. Outcome: BAM files and summary data
1. `Paleomix_settings_modern_cod.yaml`: Settings for handling raw modern reads in Paleomix (Adapter removal, Collapse of Paired End reads, BWA, and Duplicate Removal) Requires (paired) fasta.gz and reference genome Outcome: BAM files and summary data
1. `Remove_clipped_reads.sh`: Script to remove clipped reads from BAM file. Any clipped reads (including entire pairs for which one reads was clipped) will has been removed. Required: BAM files. Outcome: BAM file with clipped reads removed.

### SNP calling
1. `Haplotype_caller_GATK.sh`: Running GATK's haplotype caller. Requires BAM files, reference genome (as used to create BAM) Outputs a g.vcf.gz file per individual for further analyses
1. `Genotype_Caller.sh`: Running GATK's genotype caller. Note we ran this per chromosome in parallel jointly on all .g.vcf files. Requires g.vcf.gz, list of chromosomes, reference genome. Outputs Unfiltered vcf.gz 
1. `ANGSD_snp_calling`: SNP calling using ANGSD. We followed N. O. Therkildsen, et al., Contrasting genomic shifts underlie parallel phenotypic evolution in response to fishing. *Science* 365, 487–490 (2019), for which N.O.T. kindly shared a template script. Specifically, we identified SNP locations that can be used for faster handling later on (and we used a subset of these based on GATK filtered data). Requires BAM files and reference genome. Outputs a range of standard ANGSD output files for further analyses, including a location file. This file is edited slightly following N.O.T.'s examples.

### Filtering of VCF files
1. `Create_mappability_track.sh`: We used the [GEM mapper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377) to create mappability       tracks that can be uniquely mapped. K-mer size determines the size of the unique regions. Requires reference genome and choice of K-mer size. Outputs BED track with unique sections of the genome. 
1. `Filter_vcf_files.sh`: Filter settings for GATK VCF file as described. Please note that the script cannot be run in a single iteration: various subscripts (such as the binomial test or variance calculations) have to be run at intermediate stages. These scripts generate lists of SNPs to be ex- or included. Requires Unfiltered VCF. Outputs filtered VCF.
1. `allele_balance_calc.R`: Script to calculate bias in the number of allelic variants in heterozygote genotypes. I.e. is there consistent bias for specific variants (expection is 50/50). Requires VCF GT and DP files. Wrapper script. Outputs list of SNPs to be filtered.
1. `allele_balance_plot.R`: Script to plot histograms of read-depth of alternate alleles for heterozygote genotypes. Requires VCF GT and DP files. Outputs a histogram. 
1. `filter_allele_balance.r`: Script that performs a binomial test on heterozygote genotype calls. Reequires VCF GT and DP files. Outputs list of SNPs to be filtered.
1. `filter_allele_balance.saga.sh`: Wrapper script to submit the `filter_allele_balance.r` to the SAGA compute cluster. Requires VCFfile, `allele_balance_calc.R`. Outputs list of SNPs to be filtered at certain thresholds.
1. `Calculate_variance2.sh`: Script that calculates a measure of read-depth variation between individuals on a specfic site. Requires VCF file. Outputs list of SNPs with a measure of variation in read depth between individuals. Top 5% were filtered.

## Basic ANGSD analyses
1. `Resample_ANGSD_slurm`: Not used. Random resampling of individuals. We resample to get bootstrap intervals of windowed Fst by randomly resampling the temporal comparisons. Requires BAM files, reference genome, list of individuals. Outputs population frequencies based on randomly resampled individuals. 
1. `ANGSD_pop_freq.sh`: Calculating population specific frequencies. Requires BAM files, reference genome, list of individuals. Outputs population frequencies 
1. `ngsLD_bypop.sh`: Calculate linkage disequilibrium using [ngsLD](https://github.com/fgvieira/ngsLD). Needs list of high qualty SNPs (`data_2020.05.07/GATK_filtered_SNP_no_dam2.tab`) and beagle files. Outputs `ld.*.gatk.gz`
1. `ngsLD_decay.r`: Make a plot of the genome-wide average LD decay. Run after `ngsLD_bypop.sh`. Writes a figure.
1. `ngsLD_find_blocks.r`: Finds blocks of linked SNPs. Run after `ngsLD_bypop.sh`. Writes `analysis/ld.blocks.gatk.nodam.csv.gz`.
1. `ngsLD_find_blocks.sh`: Submits a job with `ngsLD_find_blocks.r`.
1. `ngsLD_find_unlinked.r`: Identifies unlinked loci for use in downstream analyses. Run after `nsdLD_find_blocks.sh`. Writes `analysis/ld.unlinked.*.gatk.nodam.csv.gz`.

## Basic GATK analyses
1. `ld_decay.r`: Not used in the end. Plot linkage disequilibrium decay using output from vcftools, e.g., 
    vcftools --gzvcf data_2020.05.07/Historic_dataset_no_clip.vcf.gz_HF_GQ_HWE_MISS_IND_Kmer_VAR_Binom_No_Dam2.vcf.gz --geno-r2 --ld-window-bp 5000 --keep data_2020.05.07/popLof07.txt --minGQ 15 --minDP 3 --not-chr LG01 --not-chr LG02 --not-chr LG07 --not-chr LG12 --not-chr Unplaced --chrom-map data_2019_03_18/chrom_map --out analysis/LOF_07
 
## PCA and Admixture
### from ANGSD
1. `PCAANGSD.sh`: Runs PCA from ANGSD. Requires beagle files, reference, positions file. Outputs PCA.
1. `Make_eigen_R_PCA_ANGSD.r`: Script to rescale eigen values coming from PCAANGSD. Requires PCA ANGSD output. Outputs rescaled eigenvalues.
1. `ANGSD_admixture`: Runs ANSGD admixture for different values of K. Requires BAM files, reference genome, list of individuals. Outputs ADMIXTURE results. 
1. `Plot_admixture.R`: Basic R script to plot admixture output. Requires ANGSD Admixture output. Outputs Rplot with admixture results.
1. `angsd_pcangsd_pca.sh`: Also runs PCAngsd and admixutre a few ways, including after dropping individual BM_115 and removing linked loci. Requires beagle files, output from `ngsLD_find_unlinked.r`. Outputs PCA files.
1. `angsd_pcangsd_plot_pca.R`: Plot output from `angsd_pcangsd_pca.sh`

### from GATK
1. `Filter_prune_admixture`: Script to run ADMIXTURE for a number of K on SAGA using GATK generated VCF. Requires VCF file, chromosome names for plink conversion Outputs Admixture results.
1. `Filter_Prune_PCA_LG12.sh`: Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status. Requires VCF file, chromosome names for plink conversion, inversion boundaries Outputs PCA showing inversion haplotypes.
1. `Filter_Prune_PCA_LG07.sh`: Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status. Requires VCF file, chromosome names for plink conversion, inversion boundaries Outputs PCA showing inversion haplotypes.
1. `Filter_Prune_PCA_LG02.sh`: Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status.Requires VCF file, chromosome names for plink conversion, inversion boundaries Outputs PCA showing inversion haplotypes.
1. `Filter_Prune_PCA_LG01.sh`: Script to run PCA using a GATK generated VCF. This is specifically run on an inversion location to determine haplotype status. Requires VCF file, chromosome names for plink conversion, inversion boundaries. Outputs PCA showing inversion haplotypes.
1. `Filter_Prune_PCA.sh`: Script to run PCA using a GATK generated VCF. Various filters are applied (E.g. a filter for LD, and removal of inversion chromosomes). Requires VCF file, chromosome names for plink conversion, inversion boundaries Outputs PCA showing inversion haplotypes.

## Diversity and Ne from ANGSD
1. `angsd_theta.sh`: Prints the per-site pi and calculates windowed pi and Tajima's D by population. Needs *.thetas.idx and *.thetas.gz files from realSFS and thetaStat, e.g.,
    realSFS theta/$1.saf.idx -fold 1 -P 14 > theta/$1.sfs
    realSFS saf2theta theta/$1.saf.idx -sfs theta/$1.sfs -outname theta_out/$1
    thetaStat print theta_out/$1.thetas.idx 2>/dev/null | head
    thetaStat do_stat out.thetas.idx
    cat out.thetas.idx.pestPG
    thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz
1. `angsd_theta_bypop.R`: Calculate theta and Tajima's D for the whole genome and bootstrap across linkage groups. Run by angsd_theta_bypop.sh. Needs *.pestPG.gz files from angsd_theta.sh. Writes analysis/thetas.boot.cis.csv.
1. `angsd_theta_bypop.sh`: Submits `angsd_theta_bypop.R` to the SLURM system.
1. `kjesbu_calcs.R`: Analyses fecundity, abundance, maturity, and liver condition data to calculate generation time (average age of parents) for the Northeast Arctic population. Needs files in `data\kjesbu_fecundity`. 
1. `calcNe_ANGSD.r`: Calculates effective population size (Ne) from ANGSD output using the Jorde-Ryman temporal method. Bootstraps across loci. Based on equations in the manual from [NeEstimator 2.1](http://www.molecularfisherieslaboratory.com.au/neestimator-software/).

## Selection scans
### Windowed FST
1. `angsd_fst.sh`: Calculate Fst from ANGSD output (global, by site, and windowed). Outputs *.fst.AB.gz, *.slide, *.fst.idx, *.ml, *fst.gz
1. `angsd_fst_siteshuffle_null.r`: Shuffles FST A and B values across SNP positions to calculate a null distribution of maximum FST values per genome. Run after `angsd_fst.sh`. Outputs `analysis/*.*.gatk.nodam.fst.siteshuffle.csv.gz` and `analysis/*.*.gatk.nodam.ldtrim.fst.AB.csv.gz`
1. `angsd_fst_siteshuffle_null.sh`: Submits angsd_fst_siteshuffle_null.r to SLURM.
1. `angsd_fst_siteshuffle_null_stats.r`: Examine the results from shuffling FST values across SNP positions and calculate a p-value per window. Run after `angsd_fst_siteshuffle_null.r`. Writes `output/fst_siteshuffle.angsd.gatk.csv.gz`.

### Windowed Pi and Tajima's D
1. `angsd_theta_siteshuffle_null.r`: Shuffles per-locus theta value values across SNP positions to calculate a null distribution of maximum theta values per genome. Needs *.pestPG.gz files from `angsd_theta.sh`. Outputs `analysis/theta.siteshuffle.*.csv.gz`
1. `angsd_theta_siteshuffle_null.sh`: Submits `angsd_theta_siteshuffle_null.r` to SLURM.
1. `angsd_theta_siteshuffle_null_stats.R`: Examine the results from shuffling per-locus theta values across SNP positions and calculate a p-value per window. Run after `angsd_theta_siteshuffle_null.sh`/`.r`. Writes `analysis/theta_siteshuffle.angsd.gatk.csv.gz`.

### Windowed linkage disequilibrium (LD)
1. `ngsLD_region_change.r`: Calculate changes in LD in windows. Needs output from `ngsLD_bypop.sh`. Writes `analysis/ld_change_region_5e4_ngsLD.gatk.csv.gz`.

### Shared changes across populations
1. `region_change_compare.r`: Compare Fst, LD change, Tajima's D change, pi change across populations to find regions that change substantially in more than one. Needs output from ANGSD Fst, ngsLD, and ANGSD theta change calculations. Writes `analysis/outlier_50kregions_shared_07-11-14_Can.csv.gz`.

### Wright-Fisher drift and sampling null model
#### Make the simulations
1. `wfs_byf1samp.r`: Function to do Wright-Fisher simulation with selection and a finite sample size using priors on Ne and s and specifying initial sample allele frequency. Samples at two time points. (i.e., Canada) 
1. `wfs_byf1samp_3samps.r`: Function to do Wright-Fisher simulation with selection and a finite sample size using priors on Ne and s and specifying initial sample allele frequency. Samples at three time points (i.e., Northeast Arctic).
1. `wfs_make_sims_null_function.r`: Make simulations with the null model for a particular set of sample sizes. For two samples in time (like Canada). Uses starting sample allele frequencies that match the observed data. Uses the function in `wfs_byf1samp.r`. Used by `wfs_make_sims_null_sbatch.sh`.
1. `wfs_make_sims_null.r`: Superceded by `wfs_make_sims_null_function.r`. Make simulations with the null model for a particular set of sample sizes. For two samples in time (like Canada). Uses starting sample allele frequencies that match the observed data. Uses the function in `wfs_byf1samp.r`. 
1. `wfs_make_sims_null_3times.r`: Make simulations with the null model for a particular set of sample sizes. For three samples in time (like Northeast Arctic). Uses starting sample allele frequencies that match the observed data. Uses the function in `wfs_byf1samp_3samps.r`. 
1. `wfs_make_sims_null_sbatch.sh`: Submits a SLURM job to make simulations with `wfs_make_sims_null_function.r`.
1. `wfs_make_sims_null_makebatchscript.r`: Creates batch scripts to submit all the `wfs_make_sims_null_sbatch.sh` jobs needed to make the null model simulations for a particular temporal comparison.
1. `wfs_test.r`: Code to check that the simulations match theoretical expectations.
1. `wfs_ts.r`: Function to do Wright-Fisher forward simulations with drift and selection and save all the allele frequencies through time. Used by `wfs_ts_examples.r`.
1. `wfs_ts_examples.r`: Plot a null model and a selection simulation example. Uses `wfs_ts.r`.
1. `wfs_nullmodel_compare_sims_and_obs.r`: Not used. Code to compare the simulated and observed allele frequency changes and check that the simulations are appropriate.
1. `wfs.r`: Not used. Function to do Wright-Fisher forward simulations with drift and selection and save initial and final allele frequencies. Uses priors on Ne, s, and initial true allele frequency (as opposed to initial sample allele frequency in `wfs_byf1samp.r`). 

#### Calculate p-values
1. `wfs_nullmodel_function.r`: Helper script that compares null model simulations to the observations and calculates a per-locus p-value. Runs in parallel across many loci at once. 
1. `wfs_nullmodel_sbatch.sh`: Submits `wfs_nullmodel_function.r` to run for a set of loci in a population comparison with two time points.
1. `wfs_nullmodel_sbatch_3times.sh`: Submits `wfs_nullmodel_function.r` to run for a set of loci in a population comparison with three time points.
1. `wfs_nullmodel_makesbatchscript.r`: Creates batch scripts to submit all the `wfs_nullmodel_sbatch.sh` jobs for a particular temporal comparison. All of those scripts should then be run to calculate the per-locus p-values.
1. `wfs_nullmodel_combine.r`: 
1. `wfs_nullmodel_lowp_rerun.r`:
1. `wfs_nullmodel_analysis.r`: Examines the frequencies at which the null model produces results as extreme as our observations. Run after `wfs_nullmodel_function.r` and `wfs_nullmodel_combine.r`/`wfs_nullmodel_lowp_rerun.r`. Writes `analysis/wfs_nullmodel_padj.csv.gz`.
1. `wfs_nullmodel_analysis_Canada,1907-2011&2014.r`: Not used. Messy code to examine the p-values.

### PCAngsd
1. `angsd_cut_beagle_bypop.sh`: Divides the beagle genotype file by population comparison, to get ready for running PCAngsd. Needs beagle file from ANGSD trimmed to the GATK high quality SNPs (`All_ind_beagle.GATK_no_dam.gz`). Outputs beagle.gz files.
1. `angsd_pcangsdoutlier.sh`: Run the PCAngsd selection scan. Needs beagle files, output from `ngsLD_find_unlinked.r`. Outputs *.cov, *.selection.npy, *.sites.
1. `angsd_pcangsd_plot_selection.R`: Run after `angsd_pcangsdoutlier.sh` to plot the output.

### Examine potential outliers
1. `split_fasta.sh`: Split the reference genome fasta file apart by linkage group. Output is used by `annotate_outliers.r`.
1. `annotate_outliers.r`: Compiles annotation and other information about potential outlier SNPs. Needs the VCF file (`All_rerun_hist.vcf.gz_HF_GQ_HWE_MISS_0.6_IND_Norway.vcf.gz`) and the reference genome FASTA split apart by linkage group (see `split_fasta.sh`), plus `analysis/pcangsd_outlier.gatk.nodam.unlinked.csv.gz` from running PCAngsd, `analysis/wfs_nullmodel_padj.csv.gz` from the WFS null model, `output/fst_siteshuffle.angsd.gatk.csv.gz` from the windowed FST outlier scan, `analysis/theta_siteshuffle.angsd.gatk.csv.gz` from the windowed pi and Tajima's D outlier scan, and `analysis/outlier_50kregions_shared_07-11-14_Can.csv.gz` from the shared outlier regions across populations. Writes `tables/outlier_annotation.csv` and `output/outlierSNPs_in_outlierregions.csv`

### Power analyses with SLiM
1. `slim_makeburnin.py`: Helper script to make coalescent burn-ins for the SLiM simulations. Writes a vcf file.
1. `slim_runburnins.sh`: Set up and run the SLiM burn-ins by calling `slim_makeburnin.py`. Writes a set of vcf files.
1. `slim_softsweep.slim`: SLiM code to run a forward simulation and output vcf file samples before and after soft sweep.
1. `slim_run.R`: Set up and run the SLiM simulations by calling `slim_softsweep.slim`. Checks which simulation output already exists so that those are not run again. Uses burnin files from `slim_makeburnin.py`. Outputs vcf files.
1. `fst_reynolds_fromvcf.R`: Helper script. Calculates Reynolds FST from .vcf files and writes a table of per-locus A and B values. Useful for processing output from SLiM. ANGSD also uses Reynolds FST.
1. `slim_calcfst_reynolds.sh`: Calculate Reynolds FST components for slim output. Processes vcf files from SLiM and calls `fst_reynolds_fromvcf.R`.
1. `slim_fst_reynolds_plot.R`: Plot per-locus Reynolds Fst as output by `slim_calcfst_reynolds.sh`.
1. `slim_calcfst.sh`: Not used. Calculate sliding window FST using vcftools on SLiM vcf output. Writes slim_sim_n*_s*_f*_i*.windowed.weir.fst.
1. `slim_fst_plot.R`: Not used. Makes FST Manhattan plots from vcftools output to visualize the effect of the sweeps and make sure the sims worked well. Run after `slim_calcfst.sh`.
1. `vcftobeagle.R`: Helper script to convert a vcf file to beagle format so that PCAngsd can run. Used by `slim_calcpcangsdoutlier.sh`.
1. `slim_calcpcangsdoutlier.sh`: Run the PCAngsd outlier test on the vcf output from SLiM. Writes analysis/slim_sim/*.cov, *.selection.npy, and *.sites.
1 `slim_pcangsdoutlier_plot.R`: Plots the output from `slim_calcpcangsdoutlier.sh`.
1. `slim_fst_siteshuffle_null.r`: Conduct the windowed FST outlier test on the SLiM simulations. Needs output from `slim_calcfst_reynolds.sh`. Writes `analysis/slim_sim/slim_sim_n30000_*.fst.siteshuffle.csv.gz`.
1. `slim_fst_siteshuffle_null.sh`: Calls `slim_fst_siteshuffle_null.r` across all SLiM simulations.
1. `slim_fst_siteshuffle_plot.R`: Visualize some of the results from the windowed FST outlier test. Needs output from `slim_fst_siteshuffle_null.sh`.

## Make figures
1. `figures_for_paper.r`: Makes the figures for the paper. Needs output from many of the previous analyses. Writes figures/*.png and .pdf


## Unused helper files
1. `vcf_to_trimmed_genepop.sh`: Convert vcf file to a genepop format and trim out linked loci with plink. Uses PGDSpider.
1. `vcf_to_trimmedplink.sh`: Convert vcf file to a plink format and trim out linked loci with plink. Uses PGDSpider.

# output/
1. `age_50percmature.csv`: Age at 50% mature in Canada and NE Arctic cod. From `scripts/phenotypes.r`
1. `fst_siteshuffle.angsd.gatk.csv.gz`: p-values from the windowed FST outlier test in `angsd_fst_siteshuffle_null_stats.r`. 
1. `outlierSNPs_in_outlierregions.csv`: A list of SNPs with FST>0.2 in outlier regions. Produced by `annotate_outliers.r`.

# figures/
The script-produced figures from the paper

# tables/
1. `outlier_annotation.csv`: Annotated list of potential outlier SNPs and regions. Produced by `annotate_outliers.r`. Used by `figures_for_paper.r`
1. `tableS5.csv`: Nicely formated information on potential outlier SNPs and regions. Produced by `figures_for_paper.r`.