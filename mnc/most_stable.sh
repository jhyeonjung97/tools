#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/kisti/3_MNC/0_clean/*d/*_*
do
    cd $dir; pwd
    all_done=true
    for sub_sub_dir in ./*/*
    do
        if [[ ! -f $sub_sub_dir/DONE ]]; then
            echo "DONE file is missing in $sub_sub_dir/DONE"
            all_done=false
        fi
    done
    if $all_done; then
        mkdir -p most_stable
        cd most_stable
        dzs=(0 1 2 3 4 5 6)
        for dz in ${dzs[@]}
        do
            mkdir -p ${dz}_
        done
        cd $dir
        for dz in ${dzs[@]}
        do
            lowest_dir=''
            lowest_energy=0
            for sub_dir in *_*S/
            do
                energy=$(grep -oP 'ENERGY\s+\K[-+]?[0-9]*\.?[0-9]+' "${sub_dir}${dz}_/DONE")
                if [[ $(echo "$energy < $lowest_energy" | bc) -eq 1 ]]; then
                    lowest_energy=$energy
                    lowest_dir=${sub_dir}${dz}_/
                fi
            done
            if [ -n "$lowest_dir" ]; then
                echo "${lowest_dir}restart.json most_stable/${dz}_/"
                cp ${lowest_dir}restart.json most_stable/${dz}_/
                cp ${lowest_dir}WAVECAR most_stable/${dz}_/
                cp ~/bin/tools/mnc/submit.sh most_stable/${dz}_/
            else
                echo "No valid directory found for dz=${dz}"
            fi
        done
    fi
done
