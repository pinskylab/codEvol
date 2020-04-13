#!/bin/bash

# Job name: Calculate pi and Tajima's D from ANGSD output
#SBATCH --job-name=calcTheta
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-02:00:00
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
module load angsd/0.931-GCC-8.2.0-2.31.1

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# calculate pi and theta_Watterson per site (all loci)
thetaStat print data_02.04.20/Can_40.thetas.idx | gzip >analysis/thetas.Can_40.pestPG.gz
thetaStat print data_02.04.20/Can_14.thetas.idx | gzip >analysis/thetas.Can_14.pestPG.gz
thetaStat print data_02.04.20/Lof_07.thetas.idx | gzip >analysis/thetas.Lof_07.pestPG.gz
thetaStat print data_02.04.20/Lof_11.thetas.idx | gzip >analysis/thetas.Lof_11.pestPG.gz
thetaStat print data_02.04.20/Lof_14.thetas.idx | gzip >analysis/thetas.Lof_14.pestPG.gz

# calculate windowed pi and Tajima's D by population (all loci)
thetaStat do_stat data_02.04.20/Can_40.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Can_40
thetaStat do_stat data_02.04.20/Can_14.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Can_14
thetaStat do_stat data_02.04.20/Lof_07.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Lof_07
thetaStat do_stat data_02.04.20/Lof_11.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Lof_11
thetaStat do_stat data_02.04.20/Lof_14.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Lof_14

# calculate pi and theta_Watterson per site (gatk loci)
thetaStat print data_02.04.20/GATK_Can_40.thetas.idx | gzip >analysis/thetas.Can_40.gatk.pestPG.gz
thetaStat print data_02.04.20/GATK_Can_14.thetas.idx | gzip >analysis/thetas.Can_14.gatk.pestPG.gz
thetaStat print data_02.04.20/GATK_Lof_07.thetas.idx | gzip >analysis/thetas.Lof_07.gatk.pestPG.gz
thetaStat print data_02.04.20/GATK_Lof_11.thetas.idx | gzip >analysis/thetas.Lof_11.gatk.pestPG.gz
thetaStat print data_02.04.20/GATK_Lof_14.thetas.idx | gzip >analysis/thetas.Lof_14.gatk.pestPG.gz

# windowed pi and D for only GATK loci
thetaStat do_stat data_02.04.20/GATK_Can_40.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Can_40.gatk
thetaStat do_stat data_02.04.20/GATK_Can_14.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Can_14.gatk
thetaStat do_stat data_02.04.20/GATK_Lof_07.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Lof_07.gatk
thetaStat do_stat data_02.04.20/GATK_Lof_11.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Lof_11.gatk
thetaStat do_stat data_02.04.20/GATK_Lof_14.thetas.idx -win 50000 -step 10000 -type 2 -outnames analysis/thetas.windowed.Lof_14.gatk
