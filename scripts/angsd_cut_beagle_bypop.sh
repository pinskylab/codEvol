#!/bin/bash

# Job name: Cut angsd Beagle file into populations
#SBATCH --job-name=cutBeagle
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-08:00:00
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

# Canada
zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,4-66,271-342 | bgzip -i -I data_31_01_20/Can_ind.beagle.gz.gzi -c >data_31_01_20/Can_ind.beagle.gz

# Lof 07-11
zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,67-132,133-204 | bgzip -i -I data_31_01_20/Lof0711_ind.beagle.gz.gzi -c >data_31_01_20/Lof0711_ind.beagle.gz

# Lof 07-14
zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,67-132,205-270 | bgzip -i -I data_31_01_20/Lof0714_ind.beagle.gz.gzi -c >data_31_01_20/Lof0714_ind.beagle.gz


# Lof 11-14
zcat data_31_01_20/All_ind.beagle.gz | cut -f 1-3,133-204,205-270 | bgzip -i -I data_31_01_20/Lof1114_ind.beagle.gz.gzi -c >data_31_01_20/Lof1114_ind.beagle.gz

