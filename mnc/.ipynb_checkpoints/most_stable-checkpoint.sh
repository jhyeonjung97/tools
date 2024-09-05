#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/kisti/3_MNC/0_clean/*d/*_*
do
    cd $dir; pwd
    dzs=(0 1 2 3 4 5 6)
    all_done=true
    for i in {0..6}; do
        if [[ ! -f "${i}_/DONE" ]]; then
            echo "DONE file is missing in ${i}_/"
            all_done=false
        fi
    done
    if $all_done; then
        echo "All DONE files are present."
    else
        echo "Some DONE files are missing."
    fi
    
    # mkdir -p most_stable
    # cd most_stable
    # for dz in ${dzs[@]}
    # do
    #     mkdir -p ${dz}_
    # done
    # cd $dir_now
    
    # for dz in ${dzs[@]}
    # do
    #     lowest_dir=''
    #     lowest_energy=0
        
    #     for dir in *_*S/
    #     do
    #         if [ -d ${dir}${dz}_/ ]; then
    #             if [ -f ${dir}${dz}_/DONE ]; then
    #                 energy=$(grep -oP 'ENERGY\s+\K[-+]?[0-9]*\.?[0-9]+' "${dir}${dz}_/DONE")
    #                 if [[ $(echo "$energy < $lowest_energy" | bc) -eq 1 ]]; then
    #                     lowest_energy=$energy
    #                     lowest_dir=${dir}${dz}_/
    #                 fi
    #             fi
    #         fi
    #     done
        
    #     if [ -n "$lowest_dir" ]; then
    #         echo "${lowest_dir}restart.json" "most_stable/${dz}_/"
    #         cp "${lowest_dir}restart.json" "most_stable/${dz}_/"
    #         cp "${lowest_dir}WAVECAR" "most_stable/${dz}_/"
    #         # cp ~/bin/tools/mnc/submit.sh "most_stable/${dz}_/"
    #     else
    #         echo "No valid directory found for dz=${dz}"
    #     fi
    # done
done
