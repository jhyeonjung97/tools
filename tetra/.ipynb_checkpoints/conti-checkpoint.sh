#!/bin/bash

squeue --me > ~/mystat.txt
# cat ~/kisti.txt ~/nersc.txt > ~/mystat.txt

for dir in /pscratch/sd/j/jiuy97/3_V_shape/kisti/5_V_bulk/*_*_*/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    
    if [[ -n $(grep "${coord}${row}${numb} " ~/mystat.txt) ]] || [[ -n $(grep "${coord}${row}${numb}t" ~/mystat.txt) ]]; then
        :
    elif [[ -d opt ]] && [[ -s icohp.txt ]]; then
        :
    elif [[ -d opt ]] && [[ ! -s icohp.txt ]]; then
        if [[ ! -s vasp.out ]]; then
            cp opt/restart.json opt/WAVECAR .
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/${metal}/g" lobsterin
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "/#SBATCH -J/c\#SBATCH -J ${coord}${row}${numb}t" static.sh
            pwd; sbatch static.sh
        elif [[ -s static.sh ]] && [[ -n $(grep 'BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES' vasp.out) ]]; then
            # ~/bin/shoulder/rm_mv vasp.out
            pwd; # sbatch static.sh
        elif [[ -s static.sh ]] && [[ -n $(grep 'Call to ZHEGV failed' vasp.out) ]]; then
            # sed -i -e "s/static_bulk2.py/static_bulk2_fast.py/" static.sh
            pwd; # sbatch static.sh
        # elif [[ -n $(grep SYMPREC vasp.out) ]]; then
        #     mv opt/* .; rm -r opt; rm WAVECAR
        #     sed -i -e '/walltime/c\#PBS -l walltime=06:00:00' submit.sh
        #     pwd; qsub submit.sh
        # elif [[ -n $(grep 'plane wave coefficients changed' vasp.out) ]]; then
        #     mv opt/* .; rm -r opt; rm WAVECAR
        #     sed -i -e '/walltime/c\#PBS -l walltime=06:00:00' submit.sh
        #     pwd; qsub submit.sh
        # elif [[ -n $(grep 'Inconsistent Bravais lattice types' vasp.out) ]]; then
        #     mv opt/* .; rm -r opt; rm WAVECAR
        #     sed -i -e '/walltime/c\#PBS -l walltime=06:00:00' submit.sh
        #     pwd; qsub submit.sh
        else
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/${metal}/g" lobsterin
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "s/jobname/${coord}${row}${numb}t/" static.sh
            pwd; # sbatch static.sh
        fi
    elif [[ ! -d opt ]] && [[ -s vasp.out ]]; then
        if [[ -n $(grep CONTCAR vasp.out) ]]; then
            sed -i -e 's/m.py/m_ediff.py/' submit.sh
            sed -i -e 's/m_fast.py/m_ediff.py/' submit.sh
            sed -i -e 's/m_symprec.py/m_ediff.py/' submit.sh
            sed -i -e '/walltime/c\#PBS -l walltime=06:00:00' submit.sh
            echo -e "\e[35m$PWD\e[0m"; # sbatch submit.sh
        elif [[ -n $(grep Sub-Space-Matrix vasp.out) ]] || [[ -n $(grep EDDDAV vasp.out) ]]; then
            sed -i -e 's/m.py/m_fast.py/' submit.sh
            sed -i -e 's/m_ediff.py/m_fast.py/' submit.sh
            sed -i -e 's/m_symprec.py/m_fast.py/' submit.sh
            sed -i -e '/walltime/c\#PBS -l walltime=06:00:00' submit.sh
            echo -e "\e[35m$PWD\e[0m"; # sbatch submit.sh
        # elif [[ -n $(grep SYMPREC vasp.out) ]]; then
        #     sed -i -e 's/m.py/m_symprec.py/' submit.sh
        #     sed -i -e 's/m_ediff.py/m_symprec.py/' submit.sh
        #     sed -i -e 's/m_fast.py/m_symprec.py/' submit.sh
        #     sed -i -e '/walltime/c\#PBS -l walltime=06:00:00' submit.sh
        #     pwd; qsub submit.sh
        elif [[ -s DONE ]]; then
            mkdir opt; find . -maxdepth 1 -mindepth 1 ! -name opt -exec mv {} opt/ \;
            cp opt/restart.json opt/WAVECAR .
            cp ~/bin/tools/tetra/lobsterin .
            sed -i -e "s/X/$metal/g" lobsterin
            cp ~/bin/tools/tetra/static.sh .
            sed -i -e "s/jobname/${coord}${row}${numb}stc/g" static.sh 
            echo -e "\e[35m$PWD\e[0m"; # sbatch static.sh
        else
            echo -e "\e[35m$PWD\e[0m"
        fi
    else
        echo -e "\e[32m$PWD\e[0m"
    fi
done