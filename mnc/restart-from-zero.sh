#!/bin/bash

# Step 1: Run spread.sh on directories matching the first pattern
for dir in /scratch/x2755a09/3_MNC/*d/*_*/*_*S/; do
    cd $dir
    if [[ -s 'DONE' ]] && [[ -s 'restart.json' ]]; then
        sh ~/bin/verve/spread.sh 0_/restart.json
        sh ~/bin/verve/spread.sh 0_/WAVECAR
    else
        echo "Calculation in $dir is not completed"
    fi
done

# # Step 2: Submit jobs for specific subdirectories
# for subdir in /scratch/x2755a09/3_MNC/*d/*_*/*_*S/*_/; do
#     cd $subdir
#     IFS='/' read -r -a path_components <<< $PWD
#     dz=$(echo ${path_components[-1]} | cut -d'_' -f1)
#     if [[ $dz != 0 ]]; then
#         qsub submit.sh
#     fi
# done