#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/*_*/*d/*_*/*_*S
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    if [[ -s DONE ]] && [[ ! -d nupdown/ ]]; then
        mkdir 0_ 1_ 2_ 3_ 4_ 5_ 6_ nupdown
        sh ~/bin/verve/spread.sh submit.sh
        sh ~/bin/verve/spread.sh restart.json
        sh ~/bin/verve/spread.sh WAVECAR
        find . -maxdepth 1 -mindepth 1 ! -name opt -exec mv {} nupdown/ \;
    fi
done

# for dir in /pscratch/sd/j/jiuy97/6_MNC/*_*/*d/*_*/*_*S/*_
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