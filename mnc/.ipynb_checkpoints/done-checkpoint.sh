#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*S/*_
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ ! -s DONE ]] && [[ ! -n $(grep ${metal}${spin}${dz} ~/mystat.txt) ]]; then
        pwd
    fi
done

for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*S/relaxed
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    if [[ ! -s DONE ]] && [[ ! -n $(grep ${metal}${spin}r ~/mystat.txt) ]]; then
        pwd
    fi
done

for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/most_stable/*_
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ ! -s DONE ]] && [[ ! -n $(grep ${metal}MS${dz} ~/mystat.txt) ]]; then
        pwd
    fi
done