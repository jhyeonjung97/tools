#!/bin/bash

dir=$PWD
IFS='/' read -r -a path <<< "$PWD"
ads=${path[-1]}
numb=$(echo "${path[-2]}" | cut -d'_' -f1)
metal=$(echo "${path[-2]}" | cut -d'_' -f2)
mkdir -p most_stable

# Initialize lowest_energy to a high value for comparison
lowest_dir=''
lowest_energy=0

# Loop through the specified subdirectories
for sub_dir in {HS1,HS5,IS1,IS5,LS1,LS5}; do
    energy=$(grep -oP 'ENERGY\s+\K[-+]?[0-9]*\.?[0-9]+' "${sub_dir}/DONE")
    if [[ $(echo "$energy < $lowest_energy" | bc) -eq 1 ]]; then
        lowest_energy=$energy
        lowest_dir=$sub_dir
    fi
done

# If a lowest_dir was found, copy the necessary files
if [[ -n "$lowest_dir" ]]; then
    echo "Copying files from ${dir}/${lowest_dir} to most_stable/"
    cp "${lowest_dir}/restart.json" most_stable/
    cp "${lowest_dir}/WAVECAR" most_stable/
    cp ~/bin/tools/mnc/submit.sh most_stable/
    
    # Modify submit.sh file with the appropriate job name and time limit
    sed -i -e "/#SBATCH -t/c\#SBATCH -t 00:30:00" most_stable/submit.sh
    sed -i -e "/#SBATCH -J/c\#SBATCH -J ${numb}${metal}MS${ads}" most_stable/submit.sh
else
    echo "No valid directory found with lower energy."
fi
