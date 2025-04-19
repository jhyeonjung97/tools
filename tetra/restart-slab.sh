#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/8_V_slab/*_*_*/*/*_*
do
    cd $dir
    IFS='/' read -r -a path <<< $dir
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname=${coord}${row}${numb}
    
    if [[ -n $(squeue --me | grep $jobname) ]] || [[ -f "DONE" ]]; then
        continue
    elif [[ -f "submit.sh" ]]; then
        cp ~/bin/tools/tetra/submit_slab.sh ./submit.sh
        sed -i -e "s/jobname/${jobname}/" submit.sh
        pwd; #sbatch submit.sh
    fi
done