#!/bin/bash
# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Run pcangsd selection scan by population
#SBATCH --job-name=pca_sel
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-06:00:00
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
# Create input files by population
#####################################
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-intel-2018b
module load Java/1.8.0_212 


# Canada (n = 21 and 24, individuals #1-10, 12-21 and #90-113 (see data_31_01_20/List_to_malin.tab)
# Not using BM_115 (#11) since a close relative of BM_111 (#9)
# indiv X has cols 4+3*(X-1) to 6+3*(X-1)
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,4-33,37-66,271-342 | bgzip -i -I tmp/Can_ind.beagle.gz.gzi -c >data_31_01_20/tmp.beagle.gz # all SNPs
zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,4-33,37-66,271-342 | bgzip -i -I tmp/Can_ind.GATK.beagle.gz.gzi -c >tmp/Can_ind.GATK.beagle.gz # GATK SNPs

zcat data_31_01_20/All_ind_beagle.GATK_no_dam2.gz | cut -f 1-3,4-33,37-66,271-342 | bgzip -i -I tmp/Can_ind.GATK_no_dam2.beagle.gz.gzi -c >tmp/Can_ind.GATK_no_dam2.beagle.gz # GATK SNPs filtered for no aDNA damage

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,4-33,37-66,271-342 | bgzip -i -I tmp/Can_ind.GATK_ex_inv.beagle.gz.gzi -c >tmp/Can_ind.GATK_ex_inv.beagle.gz # GATK SNPs w/out inversions

# Lof 07-11 (n = 22 and 24, individuals #22-43 and #44-67)
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I tmp/Lof0711_ind.beagle.gz.gzi -c >tmp/Lof0711_ind.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I tmp/Lof0711_ind.GATK.beagle.gz.gzi -c >tmp/Lof0711_ind.GATK.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_no_dam2.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I tmp/Lof0711_ind.GATK_no_dam2.beagle.gz.gzi -c >tmp/Lof0711_ind.GATK_no_dam2.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I tmp/Lof0711_ind.GATK_ex_inv.beagle.gz.gzi -c >tmp/Lof0711_ind.GATK_ex_inv.beagle.gz


# Lof 07-14 (n = 22 and 22, individuals #22-43 and #68-89)
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I tmp/Lof0714_ind.beagle.gz.gzi -c >tmp/Lof0714_ind.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I tmp/Lof0714_ind.GATK.beagle.gz.gzi -c >tmp/Lof0714_ind.GATK.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_no_dam2.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I tmp/Lof0714_ind.GATK_no_dam2.beagle.gz.gzi -c >tmp/Lof0714_ind.GATK_no_dam2.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I tmp/Lof0714_ind.GATK_ex_inv.beagle.gz.gzi -c >tmp/Lof0714_ind.GATK_ex_inv.beagle.gz


# Lof 11-14 (n = 24 and 22, individuals #44-67 and #68-89)
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I tmp/Lof1114_ind.beagle.gz.gzi -c >tmp/Lof1114_ind.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I tmp/Lof1114_ind.GATK.beagle.gz.gzi -c >tmp/Lof1114_ind.GATK.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_no_dam2.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I tmp/Lof1114_ind.GATK_no_dam2.beagle.gz.gzi -c >tmp/Lof1114_ind.GATK_no_dam2.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I tmp/Lof1114_ind.GATK_ex_inv.beagle.gz.gzi -c >tmp/Lof1114_ind.GATK_ex_inv.beagle.gz



###########################################
# Trim to unlinked loci (only for GATK loci)
###########################################

# make list of unlinked loci in correct format
zcat analysis/ld.unlinked.Can.gatk.csv.gz | tail -n +2 | awk --field-separator ',' '{print $2"_"$3}' > tmp/Can.gatk.unlink.beagle.pos
zcat analysis/ld.unlinked.Lof.gatk.csv.gz | tail -n +2 | awk --field-separator ',' '{print $2"_"$3}' > tmp/Lof.gatk.unlink.beagle.pos

# trim loci for Can. output has no header
zcat tmp/Can_ind.GATK.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Can.gatk.unlink.beagle.pos | bgzip -i -I tmp/Can_ind.GATK.unlink.beagle.gz.gzi >tmp/Can_ind.GATK.unlink.beagle.gz
zcat tmp/Can_ind.GATK_no_dam2.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Can.gatk.unlink.beagle.pos | bgzip -i -I tmp/Can_ind.GATK_no_dam2.unlink.beagle.gz.gzi >tmp/Can_ind.GATK_no_dam2.unlink.beagle.gz
zcat tmp/Can_ind.GATK_ex_inv.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Can.gatk.unlink.beagle.pos | bgzip -i -I tmp/Can_ind.GATK_ex_inv.unlink.beagle.gz.gzi >tmp/Can_ind.GATK_ex_inv.unlink.beagle.gz

# trim loci for Lof 0711
zcat tmp/Lof0711_ind.GATK.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof0711_ind.GATK.unlink.beagle.gz.gzi >tmp/Lof0711_ind.GATK.unlink.beagle.gz
zcat tmp/Lof0711_ind.GATK_no_dam2.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof0711_ind.GATK_no_dam2.unlink.beagle.gz.gzi >tmp/Lof0711_ind.GATK_no_dam2.unlink.beagle.gz
zcat tmp/Lof0711_ind.GATK_ex_inv.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof0711_ind.GATK_ex_inv.unlink.beagle.gz.gzi >tmp/Lof0711_ind.GATK_ex_inv.unlink.beagle.gz

# trim loci for Lof 0714
zcat tmp/Lof0714_ind.GATK.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof0714_ind.GATK.unlink.beagle.gz.gzi >tmp/Lof0714_ind.GATK.unlink.beagle.gz
zcat tmp/Lof0714_ind.GATK_no_dam2.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof0714_ind.GATK_no_dam2.unlink.beagle.gz.gzi >tmp/Lof0714_ind.GATK_no_dam2.unlink.beagle.gz
zcat tmp/Lof0714_ind.GATK_ex_inv.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof0714_ind.GATK_ex_inv.unlink.beagle.gz.gzi >tmp/Lof0714_ind.GATK_ex_inv.unlink.beagle.gz

# trim loci for Lof 1114
zcat tmp/Lof1114_ind.GATK.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof1114_ind.GATK.unlink.beagle.gz.gzi >tmp/Lof1114_ind.GATK.unlink.beagle.gz
zcat tmp/Lof1114_ind.GATK_no_dam2.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof1114_ind.GATK_no_dam2.unlink.beagle.gz.gzi >tmp/Lof1114_ind.GATK_no_dam2.unlink.beagle.gz
zcat tmp/Lof1114_ind.GATK_ex_inv.beagle.gz | java -jar scripts/filterlines.jar 1 tmp/Lof.gatk.unlink.beagle.pos | bgzip -i -I tmp/Lof1114_ind.GATK_ex_inv.unlink.beagle.gz.gzi >tmp/Lof1114_ind.GATK_ex_inv.unlink.beagle.gz



##############
# Run PCAngsd
##############
module --quiet purge  # Reset the modules to the system default
module load PCAngsd/200115-foss-2019a-Python-2.7.15 # 0.982


# Canada
#pcangsd.py -beagle tmp/Can_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can -sites_save # all loci
pcangsd.py -beagle tmp/Can_ind.GATK.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can_gatk -sites_save # gatk loci
pcangsd.py -beagle tmp/Can_ind.GATK_no_dam2.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can_gatk_no_dam2 -sites_save # gatk loci w/out aDNA damaged sites
pcangsd.py -beagle tmp/Can_ind.GATK_ex_inv.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can_gatk_ex_inv -sites_save # gatk loci w/out inversions
pcangsd.py -beagle tmp/Can_ind.GATK.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can_gatk_unlink -sites_save # unlinked gatk loci
pcangsd.py -beagle tmp/Can_ind.GATK_no_dam2.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can_gatk_no_dam2_unlink -sites_save # unlinked gatk no dam loci
pcangsd.py -beagle tmp/Can_ind.GATK_ex_inv.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_can_gatk_ex_inv_unlink -sites_save # unlinked gatk loci w/out inversions

# Lof 07-11
#pcangsd.py -beagle tmp/Lof0711_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711 -sites_save
pcangsd.py -beagle tmp/Lof0711_ind.GATK.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711_gatk -sites_save
pcangsd.py -beagle tmp/Lof0711_ind.GATK_no_dam2.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711_gatk_no_dam2 -sites_save
pcangsd.py -beagle tmp/Lof0711_ind.GATK_ex_inv.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711_gatk_ex_inv -sites_save
pcangsd.py -beagle tmp/Lof0711_ind.GATK.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711_gatk_unlink -sites_save
pcangsd.py -beagle tmp/Lof0711_ind.GATK_no_dam2.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711_gatk_no_dam2_unlink -sites_save
pcangsd.py -beagle tmp/Lof0711_ind.GATK_ex_inv.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0711_gatk_ex_inv_unlink -sites_save


# Lof 07-14
#pcangsd.py -beagle tmp/Lof0714_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714 -sites_save
pcangsd.py -beagle tmp/Lof0714_ind.GATK.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714_gatk -sites_save
pcangsd.py -beagle tmp/Lof0714_ind.GATK_no_dam2.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714_gatk_no_dam2 -sites_save
pcangsd.py -beagle tmp/Lof0714_ind.GATK_ex_inv.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714_gatk_ex_inv -sites_save
pcangsd.py -beagle tmp/Lof0714_ind.GATK.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714_gatk_unlink -sites_save
pcangsd.py -beagle tmp/Lof0714_ind.GATK_no_dam2.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714_gatk_no_dam2_unlink -sites_save
pcangsd.py -beagle tmp/Lof0714_ind.GATK_ex_inv.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof0714_gatk_ex_inv_unlink -sites_save

# Lof 11-14
#pcangsd.py -beagle tmp/Lof1114_ind.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114 -sites_save
pcangsd.py -beagle tmp/Lof1114_ind.GATK.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114_gatk -sites_save
pcangsd.py -beagle tmp/Lof1114_ind.GATK_no_dam2.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114_gatk_no_dam2 -sites_save
pcangsd.py -beagle tmp/Lof1114_ind.GATK_ex_inv.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114_gatk_ex_inv -sites_save
pcangsd.py -beagle tmp/Lof1114_ind.GATK.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114_gatk_unlink -sites_save
pcangsd.py -beagle tmp/Lof1114_ind.GATK_no_dam2.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114_gatk_no_dam2_unlink -sites_save
pcangsd.py -beagle tmp/Lof1114_ind.GATK_ex_inv.unlink.beagle.gz -selection -minMaf 0.05 -threads 2 -o analysis/pcangsd_lof1114_gatk_ex_inv_unlink -sites_save


####################################
# clean up the intermediate files
####################################
rm tmp/*.GATK.beagle.gz tmp/*.GATK.beagle.gz.gzi tmp/*.GATK_no_dam2.beagle.gz tmp/*.GATK_no_dam2.beagle.gz.gzi tmp/*.GATK_ex_inv.beagle.gz tmp/*.GATK_ex_inv.beagle.gz.gzi tmp/*.GATK.unlink.beagle.gz tmp/*.GATK.unlink.beagle.gz.gzi tmp/*.GATK_no_dam2.unlink.beagle.gz tmp/*.GATK_no_dam2.unlink.beagle.gz.gzi tmp/*.GATK_ex_inv.unlink.beagle.gz tmp/*.GATK_ex_inv.unlink.beagle.gz.gzi tmp/Can.gatk.unlink.beagle.pos tmp/Lof.gatk.unlink.beagle.pos
