#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/kisti/3_MNC/0_clean/*d/*_*/
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    all_done=true
    for sub_sub_dir in ./*/*_/
    do
        if [[ ! -s $sub_sub_dir/DONE ]]; then
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
            if [[ -n "$lowest_dir" ]]; then
                echo "${lowest_dir}restart.json most_stable/${dz}_/"
                cp ${lowest_dir}restart.json most_stable/${dz}_/
                cp ${lowest_dir}WAVECAR most_stable/${dz}_/
                cp ~/bin/tools/mnc/submit.sh most_stable/${dz}_/
                sed -i -e "s/jobname/${metal}MS${dz}/" most_stable/${dz}_/submit.sh
            else
                echo "No valid directory found for dz=${dz}"
            fi
        done
        if [[ ! -n $(grep ${metal}MS${dz} ~/mystat.txt) ]]; then
            for sub_sub_dir in ${dir}most_stable/*_/
            do
                cd $sub_sub_dir
                pwd; sbatch submit.sh
            done
        fi
    fi
done
