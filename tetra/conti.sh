#!/bin/bash/

for dir in /scratch/x2755a09/5_V_bulk/6_Octahedral_RS/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    if [[ -d opt ]] && [[ -n $(grep sec OUTCAR) ]]; then
        if [[ -d conti_1 ]]; then
            rm -r conti*/
        fi
        cp ~/bin/tools/tetra/lobsterin .
        sed -i -e "s/X/$metal/g" lobsterin
        cp ~/bin/tools/tetra/lobster.sh .
        sed -i -e "s/jobname/${coord}${row}${numb}lob/g" lobster.sh
        pwd; #qsub lobster.sh
    elif [[ ! -d opt ]] && [[ -s DONE ]] && [[ ! -n $(grep CONTCAR vasp.out) ]]; then
        if [[ -d conti_1 ]]; then
            cp conti_1/start.traj .
            rm -r conti*/
        fi
        mkdir opt; find . -maxdepth 1 -mindepth 1 ! -name opt -exec mv {} opt/ \;
        cp opt/restart.json .
        cp opt/WAVECAR .
        cp ~/bin/tools/tetra/static.sh .
        sed -i -e "s/jobname/${coord}${row}${numb}stc/g" static.sh
        pwd; #qsub static.sh
    elif [[ ! -s vasp.out ]]; then
        pwd; qstat -u x2755a09 | grep --color=auto ${coord}${row}${numb}
    else
        echo -e "\e[31m$PWD\e[0m"
    fi
done
