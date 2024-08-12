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

if [[ ${here} == 'nersc' ]]; then
    cp /pscratch/sd/j/jiuy97/6_MNC/scripts/submit.sh .
    sed -i -e "/#SBATCH -J/c\#SBATCH -J ${metal}${spin}r" submit.sh
    sbatch submit.sh
elif [[ ${here} == 'kisti' ]]; then
    cp /scratch/x2755a09/3_MNC/3d/submit.sh .
    sed -i -e "/#PBS -N/c\#PBS -N ${metal}${spin}${dz}" submit.sh
    qsub submit.sh
fi