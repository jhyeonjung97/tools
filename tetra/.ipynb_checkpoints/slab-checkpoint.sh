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
    
    if [[ -n $(grep ${coord}${row}${numb}s ~/mystat.txt) ]] || [[ -s DONE ]]; then
        :
    elif [[ $coord == 'AU' ]] || [[ $coord == 'AQ' ]]; then
        :
    elif [[ $numb == *x_* ]] || [[ $numb == *z_* ]] || [[ $numb == *s_* ]]; then
        :
    elif [[ -s vasp.out ]]; then
        if [[ -n $(grep 'WARNING: random wavefunctions but no delay for mixing, default for NELMD' vasp.out) ]] || [[ -n $(grep 'please rerun with smaller EDIFF, or copy CONTCAR' vasp.out) ]]; then
            python ~/bin/get_restart3
            cp /pscratch/sd/j/jiuy97/4_V_slab/kisti/6_V_slab/submit.sh .
            sed -i -e "/#SBATCH -J/c\#SBATCH -J ${coord}${row}${numb}s" submit.sh
            if [[ $row == 'fm' ]]; then
                sed -i -e "s/opt_slab2_afm.py/opt_slab2_fm.py/" submit.sh
            elif [[ $coord == 'NB' ]]; then
                sed -i -e "s/opt_slab2_afm.py/opt_slab2_NB.py/" submit.sh
            fi
            pwd; sbatch submit.sh
        elif [[ -n $(grep 'exceeded limit' *.e*) ]]; then
            python ~/bin/get_restart3
            cp /pscratch/sd/j/jiuy97/4_V_slab/kisti/6_V_slab/submit.sh .
            sed -i -e "/#SBATCH -J/c\#SBATCH -J ${coord}${row}${numb}s" submit.sh
            if [[ $row == 'fm' ]]; then
                sed -i -e "s/opt_slab2_afm.py/opt_slab2_fm.py/" submit.sh
            elif [[ $coord == 'NB' ]]; then
                sed -i -e "s/opt_slab2_afm.py/opt_slab2_NB.py/" submit.sh
            fi
            rm *.e*
            pwd; sbatch submit.sh
        else
            echo -e "\e[32m$PWD\e[0m"
        fi
    else
        echo -e "\e[32m$PWD\e[0m"
    fi
done