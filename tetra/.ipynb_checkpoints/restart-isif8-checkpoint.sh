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
    IFS='/' read -r -a path <<< "$dir"
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname="${coord}${row}${numb}"
    dir_fm="/pscratch/sd/j/jiuy97/7_V_bulk/${path[-3]}/fm/${path[-1]}"
    
    if [[ $row == '3d' ]] && [[ -f "$dir/DONE" ]] && [[ -f "$dir/restart.json" ]] && [[ ! -f "${dir_fm}/start.traj" ]]; then
        cp $dir/CONTCAR $dir/submit.sh .
        echo "cp CONTCAR submit.sh $dir_fm"
        ase convert CONTCAR start.traj; rm CONTCAR
        sed -i -e "s/$jobname/${coord}fm${numb}/" -e 's/afm/fm/' submit.sh
        pwd; sbatch submit.sh
    fi
    
    cd $dir
    if [[ -n $(squeue --me | grep "$jobname") ]]; then
        continue
    elif [[ -s "$dir/DONE" ]]; then
        echo 'please'; continue
    else
        pwd; python ~/bin/get_restart3; sbatch submit.sh
    fi
done