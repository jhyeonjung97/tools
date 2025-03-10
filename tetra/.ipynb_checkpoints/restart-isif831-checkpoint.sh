#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/*/*_*
do
    cd $dir
    IFS='/' read -r -a path <<< $dir
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    jobname=${coord}${row}${numb}
    
    dir_fm="/pscratch/sd/j/jiuy97/7_V_bulk/${path[-3]}/fm/${path[-1]}"
    if [[ $row == '3d' ]] && [[ -f "$dir/isif8/DONE" ]] && [[ -f "$dir/isif8/restart.json" ]] && [[ ! -f "$dir_fm/start.traj" ]]; then
        cd $dir_fm; cp $dir/isif8/CONTCAR $dir/isif8/submit.sh .
        echo -e "\e[32mcp CONTCAR submit.sh $dir_fm\e[0m"
        ase convert CONTCAR start.traj; rm CONTCAR
        sh ~/bin/verve/minute.sh 30
        sed -i -e "s/$jobname/${coord}fm${numb}/" -e 's/afm/fm/' submit.sh
        pwd; sbatch submit.sh
    fi
    
    if [[ -n $(squeue --me | grep $jobname) ]] || [[ -z $(find . -maxdepth 1 -type f) ]]; then
        continue
    elif [[ $coord == 'RS' ]] || [[ -n $(grep '12:00:00' $dir/submit.sh) ]]; then
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
        continue
        # if [[ -f "start.traj" ]]; then
        #     python ~/bin/tools/tetra/get_restart3.py
        #     if [[ -f "$dir/DONE" ]]; then
        #         pwd; echo -e "\e[31mCheck this directory!\e[0m"
        #     else
        #         pwd; sbatch submit.sh
        #     fi
        # fi
    elif [[ -d "isif2" ]]; then
        continue
    elif [[ -d "isif3" ]]; then
        mkdir isif2; find . -maxdepth 1 -type f -exec mv {} isif2/ \;
        cd isif2; cp WAVECAR restart.json submit.sh $dir; cd $dir
        sed -i -e 's/\.\/opt/~\/bin\/tools\/tetra\/opt/' submit.sh
        sed -i -e 's/_symprec.py/.py/' submit.sh
        sed -i -e '/bader/d' submit.sh
        sed -i -e 's/opt_bulk2_afm/static_bulk/' submit.sh
        sed -i -e 's/opt_bulk2_fm/static_bulk/' submit.sh
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
        cd isif3; cp WAVECAR restart.json submit.sh $dir; cd $dir
        sed -i -e 's/\.\/opt/~\/bin\/tools\/tetra\/opt/' submit.sh
        sed -i -e 's/_symprec.py/.py/' submit.sh
        sed -i -e 's/opt_bulk3/opt_bulk2/' submit.sh
        echo '' >> submit.sh
        echo 'python ~/bin/verve/bader.py' >> submit.sh
        pwd; sbatch submit.sh #; ((i+=1))
    else
        mkdir isif8; find . -maxdepth 1 -type f -exec mv {} isif8/ \;
        cd isif8; cp WAVECAR restart.json submit.sh $dir; cd $dir
        sed -i -e 's/\.\/opt/~\/bin\/tools\/tetra\/opt/' submit.sh
        sed -i -e 's/_symprec.py/.py/' submit.sh
        sed -i -e 's/opt_bulk8/opt_bulk3/' submit.sh
        pwd; sbatch submit.sh #; ((i+=1))
    fi
done