#!/bin/bash

for dir in /scratch/x2755a09/3_MNC/pre/4d/high_Ueff/*_*/*_*S/*_
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    sed -i -e "/#PBS -N/c\#PBS -N ${metal}${spin}_Ueff${dz}" submit.sh
    qsub submit.sh
done