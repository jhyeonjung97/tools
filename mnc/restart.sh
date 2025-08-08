#!/bin/bash

rm *

IFS='/' read -r -a path <<< $PWD
metal=$(echo "${path[-3]}" | cut -d'_' -f2)
spin=$(echo "${path[-2]}" | cut -d'_' -f2)
dz=$(echo "${path[-1]}" | cut -d'_' -f1)

cp $1/restart.json .
python ~/bin/tools/mnc/dz.py $dz

cp $1/WAVECAR $1/CHGCAR .
ls -l WAVECAR CHGCAR

cp /scratch/x2755a09/3_MNC/3d/submit.sh .
grep --color=auto mnc-sol submit.sh

if [[ ${here} == 'nersc' ]]; then
    cp /pscratch/sd/j/jiuy97/6_MNC/scripts/submit.sh .
    sed -i -e "/#SBATCH -J/c\#SBATCH -J ${metal}${spin}r" submit.sh
    sbatch submit.sh
elif [[ ${here} == 'kisti' ]]; then
    cp /scratch/x2755a09/3_MNC/3d/submit.sh .
    sed -i -e "/#PBS -N/c\#PBS -N ${metal}${spin}${dz}" submit.sh
    qsub submit.sh
fi

# for dir in /pscratch/sd/j/jiuy97/6_MNC/revision/*_O2_*/*_*/*_/;
# do
#     cd $dir
#     IFS='/' read -r -a path <<< $PWD
#     spin=$(echo "${path[-1]}" | cut -d'_' -f1)
#     sed -i -e "s/nupdown0/nupdown${spin}/" submit.sh
# done

for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/*_*/most_stable/o/
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    subdir=${path[-3]}
    new_dir=/pscratch/sd/j/jiuy97/6_MNC/revision/5_DFT+U_O/$subdir
    if [[ ! -d $new_dir ]]; then
        mkdir -p $new_dir
    fi
    cp restart.json $new_dir
done