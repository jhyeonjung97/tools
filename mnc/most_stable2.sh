#!/bin/bash

# Navigate through specified directories to find and submit jobs based on energy values.
for dir in /pscratch/sd/j/jiuy97/6_MNC/pourbaix/*_*/*; do
    # Check if "most_stable" directory already exists, and skip if it does
    if [[ -d "$dir/most_stable" ]]; then
        echo "most_stable directory already exists in $dir, skipping..."
        continue
    fi

    cd "$dir" || { echo "Failed to change directory to $dir"; continue; }
    pwd
    IFS='/' read -r -a path <<< "$PWD"
    ads=${path[-1]}
    numb=$(echo "${path[-2]}" | cut -d'_' -f1)
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    
    lowest_dir=''
    lowest_energy=0  # Initialize to 0 as requested

    # Check if all required subdirectories contain a DONE file
    all_done=true
    for sub_dir in ${dir}/*S*/; do
        if [[ ! -f "${sub_dir}/DONE" ]]; then
            echo "Skipping ${dir}, ${sub_dir} is missing DONE file."
            all_done=false
            break
        fi
    done
    
    # If any DONE file is missing, skip to the next dir
    if [[ "$all_done" = false ]]; then
        continue
    fi

    # Loop through each subdirectory to find the one with the lowest energy
    for sub_dir in ${dir}/*S*/; do
        energy=$(grep -oP 'ENERGY\s+\K[-+]?[0-9]*\.?[0-9]+' "${sub_dir}/DONE")
        # Proceed only if energy value was successfully extracted
        if [[ -n "$energy" && $(echo "$energy < $lowest_energy" | bc -l) -eq 1 ]]; then
            lowest_energy=$energy
            lowest_dir=$sub_dir
        fi
    done

    # Copy files if a valid lowest energy directory was found
    mkdir -p most_stable
    if [[ -n "$lowest_dir" ]]; then
        echo "Copying files from ${lowest_dir} to most_stable/"
        cp "${lowest_dir}/restart.json" most_stable/ || echo "Failed to copy restart.json"
        cp "${lowest_dir}/WAVECAR" most_stable/ || echo "Failed to copy WAVECAR"
        cp ~/bin/tools/mnc/submit.sh most_stable/ || echo "Failed to copy submit.sh"
        
        # Modify submit.sh with the appropriate job name and time limit
        sed -i -e "/#SBATCH -t/c\#SBATCH -t 00:30:00" most_stable/submit.sh
        sed -i -e "/#SBATCH -J/c\#SBATCH -J ${numb}${metal}MS${ads}" most_stable/submit.sh
    else
        echo "No valid directory found with lower energy."
    fi
    
    # Submit the job if submit.sh exists in most_stable
    if [[ -f most_stable/submit.sh ]]; then
        cd most_stable || { echo "Failed to change directory to most_stable"; continue; }
        sbatch submit.sh || echo "Failed to submit job"
    fi
done
