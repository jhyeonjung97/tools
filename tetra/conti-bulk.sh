#!/bin/bash

squeue --me > ~/mystat.txt
# cat ~/kisti.txt ~/nersc.txt > ~/mystat.txt

for dir in /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    
    if [[ -n $(grep "${coord}${row}${numb} " ~/mystat.txt) ]] || [[ -n $(grep "${coord}${row}${numb}t" ~/mystat.txt) ]]; then
        :
    elif [[ -d opt ]] && [[ -s icohp.txt ]]; then
        if [[ -n $(grep 'BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES' vasp.out) ]] || [[ -n $(grep 'while reading WAVECAR, plane wave coefficients changed' vasp.out) ]]; then
            cp opt/WAVECAR opt/restart.json .
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/${metal}/g" lobsterin
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "/#SBATCH -J/c\#SBATCH -J ${coord}${row}${numb}t" static.sh
            rm vasp.out
            pwd; # sbatch static.sh
        elif [[ ! -s opt/DONE ]]; then
            echo -e "\e[32m$PWD\e[0m"
        elif [[ ! -s full_relaxed.json ]]; then
            pwd; ase convert -f OUTCAR full_relaxed.json
        fi

    elif [[ -d opt ]] && [[ ! -s icohp.txt ]]; then
        if [[ ! -s vasp.out ]]; then
            cp opt/restart.json opt/WAVECAR .
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/${metal}/g" lobsterin
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "/#SBATCH -J/c\#SBATCH -J ${coord}${row}${numb}t" static.sh
            echo -e "\e[35m$PWD\e[0m"; # sbatch static.sh
        elif [[ -s static.sh ]] && [[ -n $(grep 'BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES' vasp.out) ]]; then
            # ~/bin/shoulder/rm_mv vasp.out
            echo -e "\e[35m$PWD\e[0m"; # sbatch static.sh
        elif [[ -s static.sh ]] && [[ -n $(grep 'Call to ZHEGV failed' vasp.out) ]]; then
            # sed -i -e "s/static_bulk2.py/static_bulk2_fast.py/" static.sh
            echo -e "\e[35m$PWD\e[0m"; # sbatch static.sh
        else
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/${metal}/g" lobsterin
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "s/jobname/${coord}${row}${numb}t/" static.sh
            echo -e "\e[35m$PWD\e[0m"; # sbatch static.sh
        fi
    elif [[ ! -d opt ]] && [[ -s vasp.out ]]; then
        if [[ -n $(grep CONTCAR vasp.out) ]]; then
            sed -i -e 's/m.py/m_ediff.py/' submit.sh
            sed -i -e 's/m_fast.py/m_ediff.py/' submit.sh
            sed -i -e 's/m_symprec.py/m_ediff.py/' submit.sh
            echo -e "\e[35m$PWD\e[0m"; # sbatch submit.sh
        elif [[ -n $(grep Sub-Space-Matrix vasp.out) ]] || [[ -n $(grep EDDDAV vasp.out) ]]; then
            sed -i -e 's/m.py/m_fast.py/' submit.sh
            sed -i -e 's/m_ediff.py/m_fast.py/' submit.sh
            sed -i -e 's/m_symprec.py/m_fast.py/' submit.sh
            echo -e "\e[35m$PWD\e[0m"; # sbatch submit.sh
        elif [[ -s DONE ]]; then
            mkdir opt; find . -maxdepth 1 -mindepth 1 ! -name opt -exec mv {} opt/ \;
            cp opt/restart.json opt/WAVECAR .
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/$metal/g" lobsterin
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "s/jobname/${coord}${row}${numb}t/" static.sh 
            echo -e "\e[35m$PWD\e[0m"; # sbatch static.sh
        else
            python ~/bin/get_restart3
            cp ~/bin/tools/tetra/submit-bulk.sh submit.sh
            sed -i -e "s/jobname/${coord}${row}${numb}/" submit.sh 
            echo -e "\e[35m$PWD\e[0m"; # sbatch submit.sh
        fi
    else
        echo -e "\e[32m$PWD\e[0m"
    fi
done