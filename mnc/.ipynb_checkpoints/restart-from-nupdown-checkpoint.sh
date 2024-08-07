#!/bin/bash

IFS='/' read -r -a path_components <<< $PWD
metal=$(echo "${path_components[-2]}" | cut -d'_' -f2)
spin=$(echo "${path_components[-1]}" | cut -d'_' -f2)

cp $1/restart.json .
cp $1/WAVECAR .
ls -l WAVECAR
cp /pscratch/sd/j/jiuy97/6_MNC/scripts/submit.sh .
sed -i -e "/#SBATCH -J/c\#SBATCH -J $metal$spin" submit.sh