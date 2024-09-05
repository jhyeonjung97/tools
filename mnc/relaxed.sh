#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/kisti/3_MNC/0_clean/*d/*_*/*_*S/
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    all_done=true
    for sub_dir in ./*_/
    do
        if [[ ! -s ${sub_dir}DONE ]]; then
            # echo "DONE file is missing in $dir$sub_sub_dir/DONE"
            all_done=false
        fi
    done
    
    if [[ $all_done == true ]]; then
        mkdir relaxed
        lowest_dir=''
        lowest_energy=0
        for sub_dir in ./*_/
        do
            energy=$(grep -oP 'ENERGY\s+\K[-+]?[0-9]*\.?[0-9]+' "${sub_dir}DONE")
            if [[ $(echo "$energy < $lowest_energy" | bc) -eq 1 ]]; then
                lowest_energy=$energy
                lowest_dir=$sub_dir
            fi
        done
        if [[ -n "$lowest_dir" ]] && [[ ! -n $(grep ${metal}${spin}rlx ~/mystat.txt) ]]; then
            cp ${lowest_dir}restart.json relaxed/
            sed -i -e "/constraints/d" relaxed/restart.json
            cp ${lowest_dir}WAVECAR relaxed/
            cp ~/bin/tools/mnc/submit.sh relaxed/
            sed -i -e "s/jobname/${metal}${spin}rlx/" relaxed/submit.sh
            cd relaxed
            pwd; sbatch submit.sh
        fi
    fi
done
