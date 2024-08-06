#!/bin/bash/

for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*/*_*/; do
    cd $dir
    if [[ -n $(grep CONTCAR vasp.out) ]] || [[ ! -s vasp.out ]]; then
        echo $PWD
        if [[ -d conti_2 ]]; then
            rm -r conti_3 conti_4 conti_5 conti_6
            cp conti_2/* .
            sed -i -e 's/fm.py/fm_ediff.py/' submit.sh
        fi
        # sh ~/bin/verve/conti.sh
    elif [[ -s DONE ]]; then
        mkdir opt
        mv * opt
        if [[ -d conti_1 ]]; then
            cp conti_1/start.traj .
        fi
        cp opt/restart.json .
        cp opt/submit.sh .
        cp opt/WAVECAR .
        sed -i -e 's/fm_ediff.py/fm.py/' submit.sh
        sed -i -e 's/opt_bulk3/static_bulk2/' submit.sh
        qsub submit.sh
    fi
done
