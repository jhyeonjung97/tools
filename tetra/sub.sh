#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/3_V_bulk/isif8/*_*_*/*/*_*/
do
    IFS='/' read -r -a path <<< $dir
    path1=${path[-1]}
    path2=${path[-2]}
    path3=${path[-3]}
    numb=$(echo $path1 | cut -d'_' -f1)
    metal=$(echo $path1 | cut -d'_' -f2)
    row=$path2
    coord=$(echo $path3 | cut -d'_' -f3)

    if [[ $row == fm ]]; then
        sed -i -e 's/afm.py/fm.py/' submit.sh
    fi
    sed -i -e "s/jobname/${coord}${row}${numb}/" submit.sh
done