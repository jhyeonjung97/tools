#!/bin/bash/

for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    if [[ -d opt ]] && [[ -s DONE ]]; then
        if [[ -d conti_1 ]]; then
            rm -r conti*/
        fi
        cp ~/bin/tools/tetra/lobsterin .
        sed -i -e "s/X/$metal/g" lobsterin
        cp ~/bin/tools/tetra/lobster.sh .
        sed -i -e "s/jobname/${coord}${row}${numb}_lobster/g" lobster.sh
        pwd; qsub lobster.sh
    elif [[ ! -d opt ]] && [[ -s DONE ]]; then
        if [[ -d conti_1 ]]; then
            cp conti_1/start.traj .
            rm -r conti*/
        fi
        mkdir opt; mv * opt
        cp opt/restart.json .
        cp opt/submit.sh .
        cp opt/WAVECAR .
        cp ~/bin/tools/tetra/static.sh
        sed -i -e "s/jobname/${coord}${row}${numb}_static/g" static.sh
        pwd; qsub static.sh
    elif [[ ! -s vasp.out ]]; then
        pwd; mystat | grep --color=auto ${coord}${row}${numb}
    else
        pwd | grep --color=auto $metal
    fi
done
