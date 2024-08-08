#!/bin/bash

IFS='/' read -r -a path_components <<< $PWD
metal=$(echo "${path_components[-3]}" | cut -d'_' -f2)
spin=$(echo "${path_components[-2]}" | cut -d'_' -f2)
dz=$(echo "${path_components[-1]}" | cut -d'_' -f1)

mkdir nupdown
mv * nupdown

cp $1/restart.json .
python ~/bin/tools/mnc/dz.py $dz

cp $1/WAVECAR .
ls -l WAVECAR

if [[ ${here} == 'nersc' ]]; then
    cp /pscratch/sd/j/jiuy97/6_MNC/scripts/submit.sh .
    sed -i -e "/#SBATCH -J/c\#SBATCH -J ${metal}${spin}r" submit.sh
elif [[ ${here} == 'kisti' ]]; then
    cp /scratch/x2755a09/3_MNC/3d/submit.sh .
    sed -i -e "/#PBS -N/c\#PBS -N ${metal}${spin}${dz}n" submit.sh