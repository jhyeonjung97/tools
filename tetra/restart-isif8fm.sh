#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/*d/*_*
do    
    if [[ "$row" == '3d' ]] && [[ -f 'DONE' ]] && [[ -f 'restart.json' ]]; then
        pwd
        dir_fm="/pscratch/sd/j/jiuy97/7_V_bulk/${path[-3]}/fm/${path[-1]}"
        mkdir -p "$dir_fm"  # Ensure destination directory exists
        cp ./restart.json ./submit.sh "$dir_fm"
        sed -i -e "s/${coord}3d${numb}/${coord}fm${numb}/" -e 's/afm/fm/' "$dir_fm/submit.sh"
    fi
done