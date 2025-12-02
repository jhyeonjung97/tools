#!/bin/bash

squeue --me > ~/mystat.txt
# cat ~/kisti.txt ~/nersc.txt > ~/mystat.txt

for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/*_*/*_*/*
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ $spin == 'stable' ]]; then
        spin='MS'
    fi
    if [[ $dz == 'relaxed' ]]; then
        dz='r'
    elif [[ $dz == 'x' ]]; then
        continue
    fi
    if [[ ! -f submit.sh ]]; then
        cp ~/bin/tools/mnc/submit.sh .
        sed -i "/#SBATCH -J/c\#SBATCH -J $metal$spin$dz" submit.sh
        if [[ $spin == 'LS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' submit.sh
        elif [[ $spin == 'IS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' submit.sh
        elif [[ $spin == 'HS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' submit.sh
        fi
    fi
    if [[ -n $(grep $metal$spin$dz ~/mystat.txt) ]]; then
        :
    elif [[ -s vasp.out ]]; then
        if [[ -s DONE ]]; then
            :
        else
            sh ~/bin/verve/correct-contcar.sh; python ~/bin/get_restart3
            # ~/bin/shoulder/rm_mv *.e* *.o* *.log
            pwd; sbatch submit.sh
            # if [[ -n $(grep 'please rerun with smaller EDIFF' vasp.out) ]]; then
            # elif [[ -n $(grep 'Call to ZHEGV failed' vasp.out) ]]; then
            # elif [[ -n $(grep 'exceeded limit' *.e*) ]]; then
            # else
            #     echo -e "\e[32m$PWD\e[0m"
            # fi
        fi
    else
        python ~/bin/tools/mnc/dz.py $dz
        echo -e "\e[35m$PWD\e[0m"; sbatch submit.sh
    fi
done

for dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/3d/*_*/*_*/*_
do
    cd $dir
    sbatch submit.sh
done