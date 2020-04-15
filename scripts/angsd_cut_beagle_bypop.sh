#!/bin/bash

# Job name: Cut angsd Beagle file into populations
#SBATCH --job-name=cutBeagle
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-03:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=2G
#
# Number of cores:
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-intel-2018b

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Canada (n = 21 and 24, individuals #1-21 and #90-113
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,4-66,271-342 | bgzip -i -I data_31_01_20/Can_ind.beagle.gz.gzi -c >data_31_01_20/Can_ind.beagle.gz # all SNPs
zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,4-66,271-342 | bgzip -i -I data_31_01_20/Can_ind.GATK.beagle.gz.gzi -c >data_31_01_20/Can_ind.GATK.beagle.gz # GATK SNPs

zcat data_31_01_20/All_ind_beagle.GATK_no_dam.gz | cut -f 1-3,4-66,271-342 | bgzip -i -I data_31_01_20/Can_ind.GATK_no_dam.beagle.gz.gzi -c >data_31_01_20/Can_ind.GATK_no_dam.beagle.gz # GATK SNPs filtered for no aDNA damage

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,4-66,271-342 | bgzip -i -I data_31_01_20/Can_ind.GATK_ex_inv.beagle.gz.gzi -c >data_31_01_20/Can_ind.GATK_ex_inv.beagle.gz # GATK SNPs w/out inversions

# Lof 07-11 (n = 22 and 24, individuals #22-43 and #44-67)
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I data_31_01_20/Lof0711_ind.beagle.gz.gzi -c >data_31_01_20/Lof0711_ind.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I data_31_01_20/Lof0711_ind.GATK.beagle.gz.gzi -c >data_31_01_20/Lof0711_ind.GATK.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_no_dam.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I data_31_01_20/Lof0711_ind.GATK_no_dam.beagle.gz.gzi -c >data_31_01_20/Lof0711_ind.GATK_no_dam.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I data_31_01_20/Lof0711_ind.GATK_ex_inv.beagle.gz.gzi -c >data_31_01_20/Lof0711_ind.GATK_ex_inv.beagle.gz


# Lof 07-14 (n = 22 and 22, individuals #22-43 and #68-89)
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I data_31_01_20/Lof0714_ind.beagle.gz.gzi -c >data_31_01_20/Lof0714_ind.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I data_31_01_20/Lof0714_ind.GATK.beagle.gz.gzi -c >data_31_01_20/Lof0714_ind.GATK.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_no_dam.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I data_31_01_20/Lof0714_ind.GATK_no_dam.beagle.gz.gzi -c >data_31_01_20/Lof0714_ind.GATK_no_dam.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I data_31_01_20/Lof0714_ind.GATK_ex_inv.beagle.gz.gzi -c >data_31_01_20/Lof0714_ind.GATK_ex_inv.beagle.gz


# Lof 11-14 (n = 24 and 22, individuals #44-67 and #68-89)
#zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I data_31_01_20/Lof1114_ind.beagle.gz.gzi -c >data_31_01_20/Lof1114_ind.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I data_31_01_20/Lof1114_ind.GATK.beagle.gz.gzi -c >data_31_01_20/Lof1114_ind.GATK.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_no_dam.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I data_31_01_20/Lof1114_ind.GATK_no_dam.beagle.gz.gzi -c >data_31_01_20/Lof1114_ind.GATK_no_dam.beagle.gz

zcat data_31_01_20/All_ind_beagle.GATK_ex_inv.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I data_31_01_20/Lof1114_ind.GATK_ex_inv.beagle.gz.gzi -c >data_31_01_20/Lof1114_ind.GATK_ex_inv.beagle.gz

