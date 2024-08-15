#!/bin/bash/

for dir in /scratch/x2755a09/5_V_bulk/7_Pyramidal_LT/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    if [[ -n $(qstat -u x2755a09 | grep --color=auto ${coord}${row}${numb}) ]]; then
        :
    elif [[ -d opt ]] && [[ -s OUTCAR ]]; then
        if [[ -n $(grep sec OUTCAR) ]]; then
            if [[ -d conti_1 ]]; then
                rm -r conti*/
            fi
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/$metal/g" lobsterin
            cp ~/bin/tools/tetra/lobster.sh .
            sed -i -e "s/jobname/${coord}${row}${numb}lob/g" lobster.sh
            pwd; #qsub lobster.sh
        else
            echo -e "\e[31m$PWD\e[0m"
        fi
    elif [[ ! -d opt ]] && [[ -s DONE ]] && [[ -s vasp.out ]]; then
        if [[ -d conti_1 ]]; then
            cp conti_1/start.traj .
            rm -r conti*/
        fi
        if [[ -n $(grep CONTCAR vasp.out) ]]; then
            sed -i -e 's/m.py/m_ediff.py/' submit.sh
            pwd; qsub submit.sh
        elif [[ -n $(grep Sub-Space-Matrix vasp.out) ]]; then
            sed -i -e 's/m.py/m_fast.py/' submit.sh
            pwd; qsub submit.sh
        elif [[ -n $(grep WARNING vasp.out) ]]; then
            echo -e "\e[31m$PWD\e[0m"
        else
            mkdir opt; find . -maxdepth 1 -mindepth 1 ! -name opt -exec mv {} opt/ \;
            cp opt/restart.json .
            cp opt/WAVECAR .
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "s/jobname/${coord}${row}${numb}stc/g" static.sh
            pwd; #qsub static.sh
        fi
    else
        echo -e "\e[31m$PWD\e[0m"
    fi
done
