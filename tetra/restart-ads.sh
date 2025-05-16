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
                if [[ ! -f 'DONE' ]] && [[ ! -f 'unmatched' ]]; then #&& [[ -z $(squeue -o "%.12j" --me | grep $jobname) ]]; then
                    echo "\e[31m$jobname\e[0m" $subdir
                    python ~/bin/get_restart3.py
                    sed -i "/#SBATCH -G/d" submit.sh
                    sed -i "/#SBATCH -N/c\#SBATCH -N 1" submit.sh
                    sed -i "/#SBATCH -C/c\#SBATCH -C cpu" submit.sh
                    sed -i "/#SBATCH -t/c\#SBATCH -t 12:00:00" submit.sh
                    sed -i "s/run_vasp_gpu2.py/run_vasp_cpu.py/" submit.sh
                    pwd; sbatch submit.sh
                fi
            fi
        done
    fi
done