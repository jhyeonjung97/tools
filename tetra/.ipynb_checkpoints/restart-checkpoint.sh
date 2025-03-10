#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/7_V_bulk/7_Pyramidal_LT/*/*_* /pscratch/sd/j/jiuy97/7_V_bulk/8_Tetrahedral_AQ/*/*_* /pscratch/sd/j/jiuy97/7_V_bulk/9_SquarePlanar_AU/*/*_*
do
    IFS='/' read -r -a path <<< "$dir"
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname="${coord}${row}${numb}"
    
    if [[ $row == '3d' ]] && [[ -f "$dir/DONE" ]] && [[ -f "$dir/restart.json" ]] && [[ ! -f "$dir_fm/start.traj" ]]; then
        cd $dir_fm; cp $dir/CONTCAR $dir/submit.sh .
        echo -e "\e[32mcp CONTCAR submit.sh $dir_fm\e[0m"
        ase convert CONTCAR start.traj; rm CONTCAR
        sh ~/bin/verve/minute.sh 30
        sed -i -e "s/$jobname/${coord}fm${numb}/" -e 's/afm/fm/' submit.sh
        pwd; sbatch submit.sh
    fi
    
    if [[ ! -f "$dir/DONE" ]] && [[ -f submit.sh ]]; then
        cd $dir; sbatch submit.sh
    fi
done