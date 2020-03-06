#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Calculate sliding window FST for slim output with vcftools
#SBATCH --job-name=fstslim
#
# Project:
#SBATCH --account=nn9244k
#
# Wall time limit: DD-HH:MM:SS
#SBATCH --time=00-3:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=1G
#
# Number of cores:
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0


# run a set of burn-in sims
scripts/slim_makeburnin.py -n 100 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n100_1.vcf
scripts/slim_makeburnin.py -n 100 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n100_2.vcf
scripts/slim_makeburnin.py -n 100 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n100_3.vcf
scripts/slim_makeburnin.py -n 100 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n100_4.vcf
scripts/slim_makeburnin.py -n 100 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n100_5.vcf

scripts/slim_makeburnin.py -n 300 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n300_1.vcf
scripts/slim_makeburnin.py -n 300 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n300_2.vcf
scripts/slim_makeburnin.py -n 300 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n300_3.vcf
scripts/slim_makeburnin.py -n 300 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n300_4.vcf
scripts/slim_makeburnin.py -n 300 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n300_5.vcf

scripts/slim_makeburnin.py -n 1000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n1000_1.vcf
scripts/slim_makeburnin.py -n 1000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n1000_2.vcf
scripts/slim_makeburnin.py -n 1000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n1000_3.vcf
scripts/slim_makeburnin.py -n 1000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n1000_4.vcf
scripts/slim_makeburnin.py -n 1000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n1000_5.vcf

scripts/slim_makeburnin.py -n 3000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n3000_1.vcf
scripts/slim_makeburnin.py -n 3000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n3000_2.vcf
scripts/slim_makeburnin.py -n 3000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n3000_3.vcf
scripts/slim_makeburnin.py -n 3000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n3000_4.vcf
scripts/slim_makeburnin.py -n 3000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n3000_5.vcf

scripts/slim_makeburnin.py -n 10000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n10000_1.vcf
scripts/slim_makeburnin.py -n 10000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n10000_2.vcf
scripts/slim_makeburnin.py -n 10000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n10000_3.vcf
scripts/slim_makeburnin.py -n 10000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n10000_4.vcf
scripts/slim_makeburnin.py -n 10000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n10000_5.vcf
