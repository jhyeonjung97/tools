#!/bin/bash

# for dir in /scratch/x2755a09/3_MNC/0_clean/*d/*_*/*_*S
# do
#     cd $dir; pwd
#     IFS='/' read -r -a path <<< $PWD
#     metal=$(echo "${path[-2]}" | cut -d'_' -f2)
#     spin=$(echo "${path[-1]}" | cut -d'_' -f2)
#     if [[ -d nupdown ]]; then
#         mkdir 1_ 2_ 3_ 4_ 5_ 6_
#         sh ~/bin/verve/spread.sh /scratch/x2755a09/3_MNC/0_clean/submit.sh
#         sh ~/bin/verve/spread.sh nupdown/restart.json
#         sh ~/bin/verve/spread.sh nupdown/WAVECAR
#         # mv nupdown 0_
#     else
#         echo -e "\e[36m$PWD\e[0m"
#     fi
# done

qstat -u x2755a09 > ~/mystat-mnc.txt
for dir in /scratch/x2755a09/3_MNC/0_clean/*d/*_*/*_*S/*_
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ -n $(grep $metal$spin$dz ~/mystat-mnc.txt) ]]; then
        :
    elif [[ -s vasp.out ]]; then
        if [[ -n $(grep 'Sub-Space-Matrix is not hermitian in DAV' vasp.out) ]]; then
            sh ~/bin/verve/correct-contcar.sh; python ~/bin/get_restart3
            sed -i -e 's/nupdown.py/nupdown-fast.py/' submit.sh
            pwd; qsub submit.sh
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
        pwd; qsub submit.sh
    fi
done