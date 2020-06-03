#!/bin/bash

# note: run this from /cluster/projects/nn9244k/in_progress/historic_malin/

# Job name: Make coalescent simulations to use as burn-in for SLiM
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
# These are standard, but they cause conda activate to fail
#set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default

# get conda set up: https://documentation.sigma2.no/apps/userinstallsw.html
module load Anaconda3/2019.03 # for conda
export PS1=\$ # Set the ${PS1} (needed in the source of the Anaconda environment)
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh # Source the conda environment setup. The variable ${EBROOTANACONDA3} comes with the module load command

# Activate the environment by using the full path (not name)
# to the environment. The full path is listed if you do
# conda info --envs at the command prompt.
conda activate /cluster/home/mlpinsky/.conda/envs/slim

# run a set of burn-in sims
for i in {1..20}
do
	scripts/slim_makeburnin.py -n 100 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n100_$i.vcf
	scripts/slim_makeburnin.py -n 300 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n300_$i.vcf
	scripts/slim_makeburnin.py -n 1000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n1000_$i.vcf
	scripts/slim_makeburnin.py -n 3000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n3000_$i.vcf
	scripts/slim_makeburnin.py -n 10000 -L 30000000 -m 3.5e-9 -r 3.11e-8 -o analysis/slim_burnin_n10000_$i.vcf
done




