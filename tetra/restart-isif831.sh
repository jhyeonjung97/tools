#!/bin/bash

i=0
for dir in /pscratch/sd/j/jiuy97/7_V_bulk/6_Octahedral_RS/3d/*_*
do
    if(( i == 5 )); then
        exit
    fi
    
    cd "$dir"
    IFS='/' read -r -a path <<< $dir
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    jobname=${coord}${row}${numb}

    if [[ -n $(squeue --me | grep $jobname) ]]; then
        ((i+=1)); continue
    elif [[ -z "$(find . -maxdepth 1 -type f ! -name 'CHGCAR' ! -name 'restart.json' ! -name 'submit.sh' ! -name 'WAVECAR' ! -name '.*')" ]]; then
        pwd; sbatch submit.sh; ((i+=1))
    elif [[ ! -f "DONE" ]]; then
        pwd; echo -e "\e[31mCheck this directory!\e[0m"
    elif [[ -d "isif3" ]]; then
        continue
    elif [[ -d "isif8" ]]; then
        mkdir isif3; find . -maxdepth 1 -type f -exec mv {} isif3/ \;
        cd isif3; cp CHGCAR WAVECAR restart.json submit.sh $dir
        cd $dir; sed -i -e 's/opt_bulk3/opt_bulk2/' submit.sh
        echo 'python ~/bin/verve/bader.py' >> submit.sh
        pwd; sbatch submit.sh; ((i+=1))
    else
        mkdir isif8; find . -maxdepth 1 -type f -exec mv {} isif8/ \;
        cd isif8; cp CHGCAR WAVECAR restart.json submit.sh $dir
        cd $dir; sed -i -e 's/opt_bulk8/opt_bulk3/' submit.sh
        pwd; sbatch submit.sh; ((i+=1))
    fi
done