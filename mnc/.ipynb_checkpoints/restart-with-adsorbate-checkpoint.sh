#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/*_*/*d/*_*/*_*S
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    adsorbate=$(echo "${path[-4]}" | cut -d'_' -f2)
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)

    if [[ adsorbate != 'clean' ]]; then
        if [[ -s final_with_calculator.json ]] && [[ -n $(grep nupdown submit.sh) ]]; then
            mkdir nupdown; mv * nupdown
            cp nupdown/restart.json nupdown/WAVECAR .
            cp ~/bin/tools/mnc/submit.sh .
            sed -i -e "s/jobname/${metal}${spin}_H/" submit.sh
            if [[ $spin == 'LS' ]]; then
                sed -i -e "s/mnc-sol.py/mnc-sol-ls-nupdown.py/" submit.sh
            elif [[ $spin == 'IS' ]]; then
                sed -i -e "s/mnc-sol.py/mnc-sol-is-nupdown.py/" submit.sh
            elif [[ $spin == 'HS' ]]; then
                sed -i -e "s/mnc-sol.py/mnc-sol-hs-nupdown.py/" submit.sh
            fi
            echo "$PWD without_nupdown"
            sbatch submit.sh
        elif [[ -n $(grep CONTCAR vasp.out) ]] && [[ ! -d conti ]]; then
            mkdir conti; mv * conti
            cp conti/restart.json conti/WAVECAR conti/submit.sh .
            echo "$PWD continue"
            sbatch submit.sh
        else
            echo "$PWD err"
        fi
    fi
done
