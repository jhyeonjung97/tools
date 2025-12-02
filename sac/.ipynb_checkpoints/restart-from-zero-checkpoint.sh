#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/2_OH/*d/*_*/*_*S/
do
    cd ${dir}; pwd
    IFS='/' read -r -a path <<< $PWD
    ads=$(echo "${path[-4]}" | cut -d'_' -f2)
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    for dz in {1..6}
    do
        mkdir ${dz}_
        cp ~/bin/tools/mnc/submit.sh ${dz}_
        cp 0_/restart.json 0_/WAVECAR ${dz}_
        sed -i "/#SBATCH -J/c\#SBATCH -J ${ads}${metal}${spin}${dz}" ${dz}_/submit.sh
        if [[ ${spin} == 'LS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' ${dz}_/submit.sh
        elif [[ ${spin} == 'IS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' ${dz}_/submit.sh
        elif [[ ${spin} == 'HS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' ${dz}_/submit.sh
        fi
        cd ${dz}_
        python ~/bin/tools/mnc/dz.py ${dz}
        cd ${dir}
    done
done