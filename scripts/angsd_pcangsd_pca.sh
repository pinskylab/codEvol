#!/bin/bash
# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Run pcangsd for all individuals
#SBATCH --job-name=pca
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-02:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=16G
#
# Number of cores:
#SBATCH --cpus-per-task=2

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error



#####################################
# Create input files
#####################################
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-intel-2018b
module load Java/1.8.0_212 

# create indices
zcat data_31_01_20/All_ind_beagle.GATK.gz | bgzip -i -I tmp/GATK.beagle.gz.gzi -c >tmp/GATK.beagle.gz # GATK SNPs

zcat data_2020.05.07/GATK_no_dam2.beagle.gz | bgzip -i -I tmp/GATK_no_dam2.beagle.gz.gzi -c >tmp/GATK_no_dam2.beagle.gz # GATK SNPs filtered for no aDNA damage

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | bgzip -i -I tmp/GATK_ex_inv.beagle.gz.gzi -c >tmp/GATK_ex_inv.beagle.gz # GATK SNPs w/out inversions


# Not using Canada BM_115 (#11) since a close relative of BM_111 (#9)
# indiv X has cols 4+3*(X-1) to 6+3*(X-1)
# cols 1-3 are basic information (marker, allele1, allele2)
zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,4-33,37-342 | bgzip -i -I tmp/GATK.norel.beagle.gz.gzi -c >tmp/GATK.norel.beagle.gz # GATK SNPs

zcat data_2020.05.07/GATK_no_dam2.beagle.gz | cut -f 1-3,4-33,37-342 | bgzip -i -I tmp/GATK_no_dam2.norel.beagle.gz.gzi -c >tmp/GATK_no_dam2.norel.beagle.gz # GATK SNPs filtered for no aDNA damage

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,4-33,37-342 | bgzip -i -I tmp/GATK_ex_inv.norel.beagle.gz.gzi -c >tmp/GATK_ex_inv.norel.beagle.gz # GATK SNPs w/out inversions




###########################################
# Trim to unlinked loci (only for GATK loci)
###########################################

# make list of unlinked loci in correct format
# uses output from ngsLD_find_unlinked.r
zcat analysis/ld.unlinked.gatk.csv.gz | tail -n +2 | awk --field-separator ',' '{print $2"_"$3}' > tmp/gatk.unlink.beagle.pos

# trim loci (all indivs). output has no header
zcat tmp/GATK.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.unlink.beagle.pos | bgzip -i -I tmp/GATK.unlink.beagle.gz.gzi >tmp/GATK.unlink.beagle.gz
zcat tmp/GATK_no_dam2.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.unlink.beagle.pos | bgzip -i -I tmp/GATK_no_dam2.unlink.beagle.gz.gzi >tmp/GATK_no_dam2.unlink.beagle.gz
zcat tmp/GATK_ex_inv.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.unlink.beagle.pos | bgzip -i -I tmp/GATK_ex_inv.unlink.beagle.gz.gzi >tmp/GATK_ex_inv.unlink.beagle.gz

# trim loci (w/out BM_115). output has no header
zcat tmp/GATK.norel.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.unlink.beagle.pos | bgzip -i -I tmp/GATK.unlink.norel.beagle.gz.gzi >tmp/GATK.unlink.norel.beagle.gz
zcat tmp/GATK_no_dam2.norel.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.unlink.beagle.pos | bgzip -i -I tmp/GATK_no_dam2.unlink.norel.beagle.gz.gzi >tmp/GATK_no_dam2.unlink.norel.beagle.gz
zcat tmp/GATK_ex_inv.norel.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/gatk.unlink.beagle.pos | bgzip -i -I tmp/GATK_ex_inv.unlink.norel.beagle.gz.gzi >tmp/GATK_ex_inv.unlink.norel.beagle.gz


##############
# Run PCAngsd: pca and admixture
##############
module --quiet purge  # Reset the modules to the system default
module load PCAngsd/200115-foss-2019a-Python-2.7.15 # 0.982

# All individuals
#pcangsd.py -beagle tmp/All_ind.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_all # all loci
pcangsd.py -beagle tmp/GATK.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk # gatk loci
pcangsd.py -beagle tmp/GATK_no_dam2.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_no_dam # gatk loci w/out aDNA damaged sites
pcangsd.py -beagle tmp/GATK_ex_inv.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_ex_inv # gatk loci w/out inversions
pcangsd.py -beagle tmp/GATK.unlink.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_unlink # unlinked gatk loci
pcangsd.py -beagle tmp/GATK_no_dam2.unlink.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_no_dam_unlink # unlinked gatk no dam loci
pcangsd.py -beagle tmp/GATK_ex_inv.unlink.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_ex_inv_unlink # unlinked gatk loci w/out inversions

# W/out BM_115
#pcangsd.py -beagle tmp/All_ind.norel.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_all.norel
pcangsd.py -beagle tmp/GATK.norel.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk.norel
pcangsd.py -beagle tmp/GATK_no_dam2.norel.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_no_dam.norel
pcangsd.py -beagle tmp/GATK_ex_inv.norel.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_ex_inv.norel
pcangsd.py -beagle tmp/GATK.unlink.norel.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_unlink.norel
pcangsd.py -beagle tmp/GATK_no_dam2.unlink.norel.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_no_dam_unlink.norel
pcangsd.py -beagle tmp/GATK_ex_inv.unlink.norel.beagle.gz -minMaf 0.05 -threads 2 -admix -o analysis/pcangsd_pca_gatk_ex_inv_unlink.norel


####################################
# clean up the intermediate files
####################################
rm tmp/*.beagle.gz tmp/*.beagle.gz.gzi tmp/gatk.unlink.beagle.pos
