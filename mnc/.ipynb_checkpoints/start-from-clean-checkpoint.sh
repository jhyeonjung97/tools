#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/1_H/*d/*_*/*_*S
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    
    clean_path="/pscratch/sd/j/jiuy97/6_MNC/0_clean/${path[-3]}/${path[-2]}/1_LS"
    if [[ -s $clean_path/re_/CONTCAR ]]; then
        cp $clean_path/re_/CONTCAR POSCAR
    elif [[ -s $clean_path/CONTCAR ]]; then
        cp $clean_path/CONTCAR POSCAR
    fi
    sed -i -e "1c\C N $metal" POSCAR
    sed -i -e "6c\C N $metal" POSCAR
    python ~/bin/tools/mnc/add-h.py

    cp ~/bin/tools/mnc/submit.sh .
    sed -i -e "s/jobname/${metal}${spin}_H/" submit.sh
    if [[ $spin == 'LS' ]]; then
        sed -i -e "s/mnc-sol.py/mnc-sol-ls-nupdown.py/" submit.sh
    elif [[ $spin == 'IS' ]]; then
        sed -i -e "s/mnc-sol.py/mnc-sol-is-nupdown.py/" submit.sh
    elif [[ $spin == 'HS' ]]; then
        sed -i -e "s/mnc-sol.py/mnc-sol-hs-nupdown.py/" submit.sh
    fi
done
