#!/bin/bash/

for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< "$PWD"
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row="${path[-2]}"
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    
    ls -t ${coord}${row}${numb}*.o* | tail -n +2 | xargs ~/bin/shoulder/rm_mv
    ls -t ${coord}${row}${numb}*.e* | tail -n +2 | xargs ~/bin/shoulder/rm_mv
done
