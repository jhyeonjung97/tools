#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*/
do
    cd $dir; # pwd
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    if [[ $spin == *stable ]]; then
        spin='MS'
    fi
    all_done=true
    for sub_dir in ./*_/
    do
        if [[ ! -s ${sub_dir}DONE ]]; then
            # echo "DONE file is missing in $dir$sub_sub_dir/DONE"
            all_done=false
        fi
    done
    if [[ $all_done == true ]] && [[ ! -f relaxed/restart.json ]]; then
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
        if [[ -n "$lowest_dir" ]]; then
            # echo $lowest_energy $lowest_dir
            cp ${lowest_dir}restart.json relaxed/
            sed -i -e "/constraints/d" relaxed/restart.json
            cp ${lowest_dir}WAVECAR relaxed/
            cp ~/bin/tools/mnc/submit.sh relaxed/
            if [[ $spin == 'LS' ]]; then
                sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' relaxed/submit.sh
            elif [[ $spin == 'IS' ]]; then
                sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' relaxed/submit.sh
            elif [[ $spin == 'HS' ]]; then
                sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' relaxed/submit.sh
            elif [[ $spin == 'stable' ]]; then
                spin='MS'
            fi
            sed -i "/#SBATCH -J/c\#SBATCH -J ${metal}${spin}r" relaxed/submit.sh
        fi
    fi
    if [[ -d relaxed ]]; then
        cd relaxed
        if [[ ! -s DONE ]] && [[ ! -n $(grep " ${metal}${spin}r" ~/mystat.txt) ]]; then
            pwd; sbatch submit.sh
        fi
    fi
done

for dir in /pscratch/sd/j/jiuy97/6_MNC/*_O*/*_*/*_*/
do
    cd $dir; # pwd
    IFS='/' read -r -a path <<< $PWD
    ads=$(echo "${path[-3]}" | cut -d'_' -f2)
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
    if [[ $all_done == true ]] && [[ ! -f relaxed/restart.json ]]; then
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
        if [[ -n "$lowest_dir" ]]; then
            # echo $lowest_energy $lowest_dir
            cp ${lowest_dir}restart.json relaxed/
            sed -i -e "/constraints/d" relaxed/restart.json
            cp ${lowest_dir}WAVECAR relaxed/
            cp ~/bin/tools/mnc/submit.sh relaxed/
            if [[ $spin == 'LS' ]]; then
                sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' relaxed/submit.sh
            elif [[ $spin == 'IS' ]]; then
                sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' relaxed/submit.sh
            elif [[ $spin == 'HS' ]]; then
                sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' relaxed/submit.sh
            elif [[ $spin == 'stable' ]]; then
                spin='MS'
            fi
            sed -i "/#SBATCH -J/c\#SBATCH -J ${ads}${metal}${spin}r" relaxed/submit.sh
        fi
    fi
    if [[ -d relaxed ]]; then
        cd relaxed
        if [[ ! -s DONE ]] && [[ ! -n $(grep "${ads}${metal}${spin}r" ~/mystat.txt) ]]; then
            pwd; sbatch submit.sh
        fi
    fi
done