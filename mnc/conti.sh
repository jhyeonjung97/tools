#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*S
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    if [[ -d nupdown ]]; then
        mkdir 1_ 2_ 3_ 4_ 5_ 6_
        sh ~/bin/verve/spread.sh nupdown/submit.sh
        sh ~/bin/verve/spread.sh nupdown/restart.json
        sh ~/bin/verve/spread.sh nupdown/WAVECAR
        mv nupdown 0_
    else
        echo -e "\e[36m$PWD\e[0m"
    fi
done

# for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*S/*_
# do
#     cd $dir; pwd
#     IFS='/' read -r -a path <<< $PWD
#     metal=$(echo "${path[-3]}" | cut -d'_' -f2)
#     spin=$(echo "${path[-2]}" | cut -d'_' -f2)
#     dz=$(echo "${path[-1]}" | cut -d'_' -f1)
#     python ~/bin/tools/mnc/dz.py $dz
#     sed -i "/#SBATCH -J/c\#SBATCH -J $metal$spin$dz" submit.sh
#     sbatch submit.sh
# done