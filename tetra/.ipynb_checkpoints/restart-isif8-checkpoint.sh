#!/bin/bash

# i=$(squeue --me | grep 'gpu' | wc -l)

# for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/*d/*_*
# do
#     if(( i > 4 )); then
#         exit
#     fi
    
#     cd "$dir"
#     IFS='/' read -r -a path <<< $dir
#     coord=$(echo "${path[-3]}" | cut -d'_' -f3)
#     row=$(echo "${path[-2]}" | cut -d'_' -f1)
#     numb=$(echo "${path[-1]}" | cut -d'_' -f1)
#     metal=$(echo "${path[-1]}" | cut -d'_' -f2)
#     jobname=${coord}${row}${numb}

#     if [[ -n $(squeue --me | grep $jobname) ]]; then
#         continue
#     elif [[ -z "$(find . -maxdepth 1 -type f ! -name 'start.traj' ! -name 'submit.sh' ! -name '.*')" ]]; then
#         pwd; sbatch submit.sh; ((i+=1))
#     fi
# done

for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/*d/*_*
do    
    cd "$dir" || continue  # Ensure script doesn't fail if cd is unsuccessful
    IFS='/' read -r -a path <<< "$dir"
    
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname="${coord}${row}${numb}"

    if [[ -n "$(squeue --me | grep "$jobname")" ]] || [[ -f 'DONE' ]]; then
        continue
    else
        pwd; python ~/bin/get_restart3; sbatch submit.sh
    fi
    
    if [[ "$row" == '3d' ]] && [[ -f 'DONE' ]] && [[ -f 'restart.json' ]]; then
        dir_fm="/pscratch/sd/j/jiuy97/7_V_bulk/${path[-3]}/fm/${path[-1]}"
        mkdir -p "$dir_fm"  # Ensure destination directory exists
        cp ./restart.json ./submit.sh "$dir_fm"
        sed -i -e 's/3d/fm/' -e 's/afm/fm/' "$dir_fm/submit.sh"
    fi
done