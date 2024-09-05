#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/4d/*_*/*_*S/
do
    cd ${dir}; pwd
    IFS='/' read -r -a path <<< $PWD
    ads=$(echo "${path[-4]}" | cut -d'_' -f2)
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    for dz in {0..6}
    do
        mkdir ${dz}_
        cp ~/bin/tools/mnc/submit.sh ${dz}_
        cp nupdown/restart.json nupdown/WAVECAR ${dz}_
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
        # sbatch submit.sh
        cd ${dir}
    done
done
# for dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/*d/*_*/*_*S/*_/
# do
#     cd $dir; pwd
#     IFS='/' read -r -a path <<< $PWD
#     ads=$(echo "${path[-5]}" | cut -d'_' -f2)
#     metal=$(echo "${path[-3]}" | cut -d'_' -f2)
#     spin=$(echo "${path[-2]}" | cut -d'_' -f2)
#     sed -i "/#SBATCH -J/c\#SBATCH -J ${ads}${metal}${spin}${dz}" submit.sh
#     sbatch submit.sh
# done
