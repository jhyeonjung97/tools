#!/bin/bash

declare -A rows
rows=(
    ["3d"]="Ti V Cr Mn Fe Co Ni Cu"
    ["4d"]="Zr Nb Mo Tc Ru Rh Pd"
    ["5d"]="Hf Ta W Re Os Ir Pt"
)

spins=("LS" "IS" "HS")
small_spins=("ls" "is" "hs")
specific_metals=("Fe" "Ni" "Mo" "Ru" "Pd" "W" "Pt")
adsorbates=("clean" "h" "o-o" "o-oh" "o" "oh-o" "oh-oh" "oh" "oho" "ohoh" "oo" "ooh")

i=0

# Loop through rows and metals
for row in "${!rows[@]}"; do
    metals=${rows[$row]}
    m=1
    for metal in $metals; do
        ((m++))

        # Check if the metal is in specific_metals
        if [[ " ${specific_metals[@]} " =~ " ${metal} " ]]; then
            ((i++))

            for s in $(seq 0 2); do
                spin=${spins[$s]}
                small_spin=${small_spins[$s]}

                for dz in "1" "5"; do
                    echo $row $metal $m $i $s
                    # Find the first matching directory
                    base_path="/pscratch/sd/j/jiuy97/6_MNC/0_clean/${row}/${m}_${metal}"
                    path=$(find "$base_path" -maxdepth 1 -type d -name "*_${spin}" 2>/dev/null)

                    # Check if path exists
                    if [ -n "$path" ]; then
                        target_dir="${path}/${dz}_"
                        echo "Processing: $target_dir"

                        if [ -d "$target_dir" ]; then
                            cd "$target_dir" || continue

                            # Loop through add*.py scripts in ~/bin/tools/mnc/ and execute them
                            for file in ~/bin/tools/mnc/add*.py; do
                                python "$file"
                            done
                        else
                            echo "Target directory does not exist: $target_dir"
                            continue
                        fi
                    else
                        echo "Path does not exist: $base_path/*_${spin}"
                        continue
                    fi
                    
                    pourbaix="/pscratch/sd/j/jiuy97/6_MNC/pourbaix/${i}_${metal}"
                    mkdir -p "$pourbaix"

                    for adsorbate in "${adsorbates[@]}"; do
                        adsorbate_dir="${pourbaix}/${adsorbate}/${spin}${dz}"
                        mkdir -p "$adsorbate_dir"
                        cd "$adsorbate_dir" || continue

                        if [[ "$adsorbate" == "clean" ]]; then
                            restart_file="${target_dir}/restart.json"
                        else
                            restart_file="${target_dir}/restart-${adsorbate}.json"
                        fi

                        if [ -f "$restart_file" ]; then
                            cp "$restart_file" .
                        else
                            echo "File not found: $restart_file"
                            continue
                        fi

                        cp ~/bin/tools/mnc/submit.sh .
                        sed -i -e "s/jobname/${metal}_${adsorbate}/" submit.sh
                        sed -i -e "s/mnc-sol/mnc-sol-${small_spin}/" submit.sh
                    done
                done
            done
        fi
    done
done
