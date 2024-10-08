#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*/*/
do
    if [[ ! -f ${dir}final_with_calculator.json ]] && [[ -s ${dir}DONE ]]; then
        cd ${dir}; pwd
        ase convert -n -1 OUTCAR final_with_calculator.json
    fi
done
