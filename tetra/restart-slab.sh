#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/8_V_slab/8_Tetrahedral_WZ/*/*_* /pscratch/sd/j/jiuy97/8_V_slab/9_Tetrahedral_ZB/*/*_*
do
    cd $dir
    IFS='/' read -r -a path <<< $dir
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname=${coord}${row}${numb}
    
    if [[ -f 'unmatched' ]]; then
        rm *.*
    else
        sed -i -e "s/jobname/${jobname}/" submit.sh
        mv slab.traj start.traj
        pwd; sbatch submit.sh
    fi
done