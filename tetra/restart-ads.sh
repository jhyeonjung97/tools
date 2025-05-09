#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/8_V_slab/*_*_*/*/*_*
do
    cd $dir    
    if [[ ! -f 'unmatched' ]]; then
        for subdir in $dir/o*/
        do
            if [[ -d $subdir ]]; then
                cd $subdir
                IFS='/' read -r -a path <<< $subdir
                coord=$(echo "${path[-4]}" | cut -d'_' -f3)
                row=$(echo "${path[-3]}" | cut -d'_' -f1)
                numb=$(echo "${path[-2]}" | cut -d'_' -f1)
                metal=$(echo "${path[-2]}" | cut -d'_' -f2)
                ads=$(echo "${path[-1]}" | cut -d'_' -f1)
                jobname=${coord}${row}${numb}${ads}
                if [[ ! -f 'DONE' ]] && [[ ! -n $(squeue --me | grep $jobname) ]]; then
                    python ~/bin/get_restart3.py
                    sed -i "/#SBATCH -t/c\#SBATCH -t 04:00:00" submit.sh
                    sed -i "/#SBATCH -q/c\#SBATCH -q regular" submit.sh
                    pwd; sbatch submit.sh
                fi
            fi
        done
    fi
done