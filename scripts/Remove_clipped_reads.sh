#!/bin/bash

# Job name:
#SBATCH --job-name=2N_GNTP
#
# Project:
#SBATCH --account=nn9244k
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=45G
#
#if you need long time
#   --partition=long

#If you need hugemem...

#S --partition=hugemem

#SBATCH --ntasks-per-node=1

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi


##Removes clipped reads, including pairs for which one read is clipped from BAM files
## $1 feeds a bam file, pre clipped.

module load samtools


mkdir -p ../bam_files_no_clipped_reads/

samtools view $1 |awk '$6 ~ /H|S/ {print $1}' |sort -u > $1.exclude_list.txt

samtools view -h $1 | fgrep -v -f $1.exclude_list.txt | samtools view -b - -o ../bam_files_no_clipped_reads/$1.no_clip.bam


samtools index ../bam_files_no_clipped_reads/$1.no_clip.bam

rm $1.exclude_list.txt
