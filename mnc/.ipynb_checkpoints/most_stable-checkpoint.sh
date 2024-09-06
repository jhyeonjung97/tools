#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/*_*/
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    echo -e "dz\tspin_state" > lowest.csv
    all_done=true
    for sub_sub_dir in ./*_*S/*_/
    do
        if [[ ! -s ${sub_sub_dir}DONE ]]; then
            # echo "DONE file is missing in $dir$sub_sub_dir/DONE"
            all_done=false
        fi
    done
    if [[ $all_done == true ]]; then
        mkdir -p most_stable
        cd most_stable
        dzs=(0 1 2 3 4 5 6)
        for dz in ${dzs[@]}
        do
            mkdir -p ${dz}_
        done
        cd $dir; pwd
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
                    lowest_sub_dir=${sub_dir}
                fi
            done
            if [[ -n "$lowest_dir" ]]; then
                product=$(echo "scale=2; $dz * 0.2" | bc)
                product=$(echo "$product" | awk '{printf "%.2f", $0}')
                clean_dir=$(echo "$lowest_sub_dir" | sed 's:/$::')
                spin=$(echo "$clean_dir" | cut -d'_' -f2)
                echo -e "$product\t$spin" >> lowest.csv
                if [[ ! -f most_stable/6_/restart.json ]]; then
                    echo "${lowest_dir}restart.json most_stable/${dz}_/"
                    cp ${lowest_dir}restart.json most_stable/${dz}_/
                    cp ${lowest_dir}WAVECAR most_stable/${dz}_/
                    cp ~/bin/tools/mnc/submit.sh most_stable/${dz}_/
                    sed -i -e "s/jobname/${metal}MS${dz}/" most_stable/${dz}_/submit.sh
                fi
            else
                echo "No valid directory found for dz=${dz}"
            fi
        done
    fi
done
