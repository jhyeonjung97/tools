#!/bin/bash

# Navigate through specified directories to find and submit jobs based on energy values.
for dir in /pscratch/sd/j/jiuy97/6_MNC/pourbaix/*_/; do
    cd "$dir" || { echo "Failed to change directory to $dir"; continue; }
    pwd
    IFS='/' read -r -a path <<< "$PWD"
    ads=${path[-1]}
    numb=$(echo "${path[-2]}" | cut -d'_' -f1)
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    mkdir -p most_stable
    
    lowest_dir=''
    lowest_energy=9999999  # Initialize to a high value

    # Define subdirectories to check
    sub_dirs=(HS1 HS5 IS1 IS5 LS1 LS5)
    
    for sub_dir in "${sub_dirs[@]}"; do
        if [[ -f "${sub_dir}/DONE" ]]; then
            energy=$(grep -oP 'ENERGY\s+\K[-+]?[0-9]*\.?[0-9]+' "${sub_dir}/DONE")
            if [[ $(echo "$energy < $lowest_energy" | bc) -eq 1 ]]; then
                lowest_energy=$energy
                lowest_dir=$sub_dir
            fi
        fi
    done
    
    if [[ -n "$lowest_dir" ]]; then
        echo "Copying files from ${dir}/${lowest_dir} to most_stable/"
        cp "${lowest_dir}/restart.json" most_stable/ || echo "Failed to copy restart.json"
        cp "${lowest_dir}/WAVECAR" most_stable/ || echo "Failed to copy WAVECAR"
        cp ~/bin/tools/mnc/submit.sh most_stable/ || echo "Failed to copy submit.sh"
        
        # Modify submit.sh with the appropriate job name and time limit
        sed -i -e "/#SBATCH -t/c\#SBATCH -t 00:30:00" most_stable/submit.sh
        sed -i -e "/#SBATCH -J/c\#SBATCH -J ${numb}${metal}MS${ads}" most_stable/submit.sh
    else
        echo "No valid directory found with lower energy."
    fi
    
    if [[ -f most_stable/submit.sh ]]; then
        cd most_stable || exit
        sbatch submit.sh || echo "Failed to submit job"
    fi
done
