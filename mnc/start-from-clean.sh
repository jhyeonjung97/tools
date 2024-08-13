#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/*d/*_*/*_*S; then
do
    cd $dir
    pwd
    
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    
    clean_path="/scratch/x2755a09/3_MNC/0_clean/${path[-3]}/${path[-2]}/1_LS"
    ase converf -f $clean_path/CONTCAR start.traj
done