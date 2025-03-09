#!/bin/bash

# i=$(squeue --me | grep 'RS' | wc -l)
for dir in /pscratch/sd/j/jiuy97/7_V_bulk/6_Octahedral_RS/*/*_*
do
    # if(( i > 4 )); then
    #     exit
    # fi
    
    cd $dir
    IFS='/' read -r -a path <<< $dir
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname=${coord}${row}${numb}

    if [[ -n $(squeue --me | grep $jobname) ]]; then
        continue
    elif [[ -z "$(find . -maxdepth 1 -type f ! -name 'start.traj' ! -name 'submit.sh' ! -name '.*')" ]]; then
        pwd; sbatch submit.sh #; ((i+=1))
    elif [[ -z "$(find . -maxdepth 1 -type f ! -name 'restart.json' ! -name 'submit.sh' ! -name 'WAVECAR' ! -name '.*')" ]]; then
        if [[ -d "isif3" ]] && [[ ! -n $(grep 'bader' submit.sh) ]]; then
            echo 'python ~/bin/verve/bader.py' >> submit.sh
        fi
        pwd; sbatch submit.sh #; ((i+=1))
    elif [[ -z "$(find . -maxdepth 1 -type f ! -name 'restart.json' ! -name 'submit.sh' ! -name 'WAVECAR' ! -name 'lobsterin' ! -name '.*')" ]]; then
        if [[ -d "isif2" ]] && [[ -n $(grep 'bader' submit.sh) ]]; then
            sed -i -e 's/opt_bulk3/opt_bulk2/' submit.sh
        fi
        pwd; sbatch submit.sh #; ((i+=1))
    elif [[ ! -f "DONE" ]]; then
        pwd; echo -e "\e[31mCheck this directory!\e[0m"
    elif [[ -d "isif2" ]]; then
        continue
    elif [[ -d "isif3" ]]; then
        mkdir isif2; find . -maxdepth 1 -type f -exec mv {} isif2/ \;
        cd isif2; cp WAVECAR restart.json submit.sh $dir
        cd $dir; sed -i -e 's/opt_bulk2_afm/static_bulk/' submit.sh
        sed -i -e '/bader/d' submit.sh
        echo '~/bin/lobster-5.0.0/lobster-5.0.0' >> submit.sh
        echo 'python ~/bin/tools/tetra/cohp.py > icohp.txt' >> submit.sh
        echo 'python ~/bin/tools/tetra/cobi.py > icobi.txt' >> submit.sh
        echo 'python ~/bin/tools/tetra/coop.py > icoop.txt' >> submit.sh
        echo 'python ~/bin/tools/tetra/spilling.py' >> submit.sh
        cp ~/bin/tools/tetra/lobsterin .
        sed -i -e "s/X/${metal}/" lobsterin
        pwd; sbatch submit.sh #; ((i+=1))
    elif [[ -d "isif8" ]]; then
        mkdir isif3; find . -maxdepth 1 -type f -exec mv {} isif3/ \;
        cd isif3; cp WAVECAR restart.json submit.sh $dir
        cd $dir; sed -i -e 's/opt_bulk3/opt_bulk2/' submit.sh
        echo '' >> submit.sh
        echo 'python ~/bin/verve/bader.py' >> submit.sh
        pwd; sbatch submit.sh #; ((i+=1))
    else
        mkdir isif8; find . -maxdepth 1 -type f -exec mv {} isif8/ \;
        cd isif8; cp WAVECAR restart.json submit.sh $dir
        cd $dir; sed -i -e 's/opt_bulk8/opt_bulk3/' submit.sh
        pwd; sbatch submit.sh #; ((i+=1))
    fi
done