#!/bin/bash

squeue --me > ~/mystat.txt
# cat ~/kisti.txt ~/nersc.txt > ~/mystat.txt

for dir in /pscratch/sd/j/jiuy97/4_V_slab/kisti/6_V_slab/*_*_*/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    
    if [[ -n $(grep ${coord}${row}${numb}s ~/mystat.txt) ]]; then
        :
    elif [[ $coord == 'AU' ]] || [[ $coord == 'AQ' ]]; then
        :
    elif [[ $numb == *x ]] || [[ $numb == *z ]] || [[ $numb == *s ]]; then
        :
    elif [[ -s DONE ]]; then
        :
        # mkdir full_relaxed
        # python ~/bin/get_restart3
        # cp restart.json full_relaxed/
        # cp /global/homes/j/jiuy97/bin/tools/tetra/submit.sh full_relaxed/
        # cd full_relaxed/
        # sed -i -e '/constraints/d' restart.json
        # sed -i -e "/#SBATCH -J/c\#SBATCH -J ${coord}${row}${numb}f" submit.sh
        # pwd; sbatch submit.sh
    elif [[ -s vasp.out ]]; then
        if [[ -n $(grep 'WARNING: random wavefunctions but no delay for mixing, default for NELMD' vasp.out) ]] || [[ -n $(grep 'please rerun with smaller EDIFF, or copy CONTCAR' vasp.out) ]] || [[ -n $(grep 'exceeded limit' *.e*) ]]; then
            python ~/bin/get_restart3
            cp /global/homes/j/jiuy97/bin/tools/tetra/submit.sh .
            sed -i -e "/#SBATCH -J/c\#SBATCH -J ${coord}${row}${numb}s" submit.sh
            if [[ $row == 'fm' ]]; then
                sed -i -e "s/opt_slab2_afm.py/opt_slab2_fm.py/" submit.sh
            elif [[ $coord == 'NB' ]]; then
                sed -i -e "s/opt_slab2_afm.py/opt_slab2_NB.py/" submit.sh
            fi
            pwd; sbatch submit.sh
        else
            echo -e "\e[32m$PWD\e[0m"
        fi
    else
        echo -e "\e[32m$PWD\e[0m"
    fi
done