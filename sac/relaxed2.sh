#!/bin/bash

root=$PWD

# Loop through the specified directories
for dir in {HS1,HS5,IS1,IS5,LS1,LS5}
do
    cd $root/$dir || continue # Ensure the directory exists

    mkdir -p relaxed
    lowest_dir=''
    lowest_energy=999999 # Set an initial high value for energy comparison

    # Determine the spin type from the directory name
    spin=${dir:0:2}

    # Loop through subdirectories ending with '/'
    for sub_dir in */; do
        energy=$(grep -oP 'ENERGY\s+\K[-+]?[0-9]*\.?[0-9]+' "${sub_dir}DONE")
        if [[ $(echo "$energy < $lowest_energy" | bc) -eq 1 ]]; then
            lowest_energy=$energy
            lowest_dir=$sub_dir
        fi
    done

    if [[ -n "$lowest_dir" ]]; then
        cp "${lowest_dir}restart.json" relaxed/
        sed -i -e "/constraints/d" relaxed/restart.json
        cp "${lowest_dir}WAVECAR" relaxed/
        cp ~/bin/tools/mnc/submit.sh relaxed/

        # Modify the submit.sh file based on the spin type
        if [[ $spin == 'LS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' relaxed/submit.sh
        elif [[ $spin == 'IS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' relaxed/submit.sh
        elif [[ $spin == 'HS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' relaxed/submit.sh
        fi

        # Update the job name in submit.sh
        sed -i "/#SBATCH -J/c\#SBATCH -J ${metal}${spin}r" relaxed/submit.sh
    fi

    # If the relaxed directory exists, submit the job
    if [[ -d relaxed ]]; then
        cd relaxed || continue
        sbatch submit.sh
    fi
done
