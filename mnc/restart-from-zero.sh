#!/bin/bash

# metals=('2_Ti' '3_V' '4_Cr' '5_Mn' '6_Fe' '7_Co' '8_Ni' '9_Cu')
# spins=('1_LS' '2_IS' '2_HS' '3_HS')
# dzs=('1_' '2_' '3_' '4_' '5_' '6_')

mkdir "$1"_
mv * "$1"_

IFS='/' read -r -a path_components <<< $PWD
metal=$(echo "${path_components[-2]}" | cut -d'_' -f2)
spin=$(echo "${path_components[-1]}" | cut -d'_' -f2)

dzs=(0 1 2 3 4 5 6)
dir_now=$PWD
for dz in ${dzs[@]}; do
    if [[ $dz != $1 ]]; then
        mkdir "$dz"_/
        cp "$1"_/restart.json "$dz"_/
        cp "$1"_/WAVECAR "$dz"_/
        cd "$dz"_/
        cp /scratch/x2755a09/3_MNC/3d/submit.sh .
        sed -i -e "/#PBS -N/c\#PBS -N $metal$spin$dz" submit.sh
        python ~/bin/tools/mnc/dz.py $dz
        qsub submit.sh
        cd $dir_now
    fi
done