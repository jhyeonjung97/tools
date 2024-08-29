#!/bin/bash

qstat -u x2755a09 > ~/mystat.txt
for dir in /scratch/x2755a09/3_MNC/0_clean/*d/*_*/*_*S/*_
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ -n $(grep $metal$spin$dz ~/mystat.txt) ]]; then
        :
    elif [[ -s vasp.out ]]; then
        if [[ -s DONE ]]; then
            :
        else
            python ~/bin/get_restart3
            if [[ -n $(grep 'please rerun with smaller EDIFF' vasp.out) ]]; then
                pwd; qsub submit.sh
            elif [[ -n $(grep 'Call to ZHEGV failed' vasp.out) ]]; then
                sed -i 's/nupdown.py/nupdown-fast.py/' submit.sh
                pwd; qsub submit.sh
            elif [[ -n $(grep 'exceeded limit' *.e*) ]]; then
                pwd; qsub submit.sh
            else
                echo -e "\e[32m$PWD\e[0m"
        fi
    else
        python ~/bin/tools/mnc/dz.py $dz
        sed -i "/#PBS -N/c\#PBS -N $metal$spin$dz" submit.sh
        if [[ $spin == 'LS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' submit.sh
        elif [[ $spin == 'IS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' submit.sh
        elif [[ $spin == 'HS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' submit.sh
        fi
        sed -i -e "s/ncpus=40/ncpus=64/" submit.sh
        sed -i -e "s/mpiprocs=40/mpiprocs=64/" submit.sh
        sed -i -e "s/run_vasp_flat.py/run_vasp.py/" submit.sh
        sed -i -e "/#PBS -l walltime/c\#PBS -l walltime=48:00:01" submit.sh
        sed -i -e "/#PBS -q/c\#PBS -q long" submit.sh
        pwd; qsub submit.sh
    fi
done