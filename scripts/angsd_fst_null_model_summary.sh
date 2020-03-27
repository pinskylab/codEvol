#!/bin/bash

# currently set up only to be run as an interactive bash session. only small tweaks would be needed to write files to disk, though.

###################
### CAN
maxvals=() # to hold the max Fst from each iteration
for i in {1..200} # for each iteration
do
	echo $i # print the iteration #
	maxvalsperlg=() # to hold the max Fst for each LG in this iteration
	for filename in analysis/Resampled_pairwise_ANGSD/analysis_Can/${i}_Can_40_13_LG??.slide # for each LG (not Unplaced)
	do
	#echo `sort -nrk5,5 $filename | head -1 | cut -f5`
	echo $filename # print filename
	maxvalsperlg+=(`sort -nrk5,5 $filename | head -1 | cut -f5`) # get highest fst from the LG. col 5 is the fst values
	done

	IFS=$'\n' # set the spacer character so that the next statement works
	maxvals+=(`echo "${maxvalsperlg[*]}" | sort -nr | head -1`) # find the max across the LGs in this iteration
	unset IFS # reset spacer char to default
done

echo "${maxvals[*]}" # whole thing

IFS=$'\n' 
echo "${maxvals[*]}" | sort -n  # whole thing sorted in a column
unset IFS # reset spacer char to default

echo "${#maxvals[@]}" # print the length: 200 iterations

IFS=$'\n'
echo "${maxvals[*]}" | sort -nr | head -1 # find the max in the array
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.95-0.5)]}' # 95th percentile
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.975-0.5)]}' # 95th percentile
unset IFS


###################
### LOF 07-11
maxvals=() # to hold the max Fst from each iteration
for i in {1..200} # for each iteration
do
echo $i # print the iteration #
maxvalsperlg=() # to hold the max Fst for each LG in this iteration
for filename in analysis/Resampled_pairwise_ANGSD/analysis_07_11/${i}_Lof_07_11_LG??.slide # for each LG (not Unplaced)
do
#echo `sort -nrk5,5 $filename | head -1 | cut -f5`
echo $filename # print filename
maxvalsperlg+=(`sort -nrk5,5 $filename | head -1 | cut -f5`) # get highest fst from the LG
done

IFS=$'\n' # set the spacer character so that the next statement works
maxvals+=(`echo "${maxvalsperlg[*]}" | sort -nr | head -1`) # find the max across the LGs in this iteration
unset IFS # reset spacer char to default
done

echo "${maxvals[*]}" # whole thing

IFS=$'\n' 
echo "${maxvals[*]}" | sort -n  # whole thing sorted in a column
unset IFS # reset spacer char to default

echo "${#maxvals[@]}" # print the length: 200 iterations

IFS=$'\n'
echo "${maxvals[*]}" | sort -nr | head -1 # find the max in the array
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.95-0.5)]}' # 95th percentile
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.975-0.5)]}' # 95th percentile
unset IFS

###################
### LOF 07-14
maxvals=() # to hold the max Fst from each iteration
for i in {1..200} # for each iteration
do
echo $i # print the iteration #
maxvalsperlg=() # to hold the max Fst for each LG in this iteration
for filename in analysis/Resampled_pairwise_ANGSD/analysis_07_14/${i}_Lof_07_14_LG??.slide # for each LG (not Unplaced)
do
#echo `sort -nrk5,5 $filename | head -1 | cut -f5`
echo $filename # print filename
maxvalsperlg+=(`sort -nrk5,5 $filename | head -1 | cut -f5`) # get highest fst from the LG
done

IFS=$'\n' # set the spacer character so that the next statement works
maxvals+=(`echo "${maxvalsperlg[*]}" | sort -nr | head -1`) # find the max across the LGs in this iteration
unset IFS # reset spacer char to default
done

echo "${maxvals[*]}" # whole thing

IFS=$'\n' 
echo "${maxvals[*]}" | sort -n  # whole thing sorted in a column
unset IFS # reset spacer char to default

echo "${#maxvals[@]}" # print the length: 200 iterations

IFS=$'\n'
echo "${maxvals[*]}" | sort -nr | head -1 # find the max in the array
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.95-0.5)]}' # 95th percentile
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.975-0.5)]}' # 95th percentile
unset IFS


###################
### LOF 11-14
maxvals=() # to hold the max Fst from each iteration
for i in {1..200} # for each iteration
do
echo $i # print the iteration #
maxvalsperlg=() # to hold the max Fst for each LG in this iteration
for filename in analysis/Resampled_pairwise_ANGSD/analysis_11_14/${i}_Lof_11_14_LG??.slide # for each LG (not Unplaced)
do
#echo `sort -nrk5,5 $filename | head -1 | cut -f5`
echo $filename # print filename
maxvalsperlg+=(`sort -nrk5,5 $filename | head -1 | cut -f5`) # get highest fst from the LG
done

IFS=$'\n' # set the spacer character so that the next statement works
maxvals+=(`echo "${maxvalsperlg[*]}" | sort -nr | head -1`) # find the max across the LGs in this iteration
unset IFS # reset spacer char to default
done

echo "${maxvals[*]}" # whole thing

IFS=$'\n' 
echo "${maxvals[*]}" | sort -n  # whole thing sorted in a column
unset IFS # reset spacer char to default

echo "${#maxvals[@]}" # print the length: 200 iterations

IFS=$'\n'
echo "${maxvals[*]}" | sort -nr | head -1 # find the max in the array
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.95-0.5)]}' # 95th percentile
unset IFS

IFS=$'\n'
echo "${maxvals[*]}" | sort -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.975-0.5)]}' # 95th percentile
unset IFS
