#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/*/*_*
do    
    IFS='/' read -r -a path <<< "$dir"
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname="${coord}${row}${numb}"
    
    if [[ -f "$dir/DONE" ]] && [[ ! -n $(grep 'PROFILE, used timers' "$dir/vasp.out") ]]; then
        cd $dir; rm DONE
        pwd; python ~/bin/tools/tetra/get_restart3.py; sbatch submit.sh
    fi
    
    if [[ ! -f "$dir/DONE" ]] && [[ -n $(grep 'PROFILE, used timers' "$dir/vasp.out") ]]; then
        cd $dir; python ~/bin/tools/tetra/get_restart3.py
        if [[ -f "$dir/DONE" ]]; then
            pwd; squeue --me | grep --color=auto $jobname
        fi
    fi
done