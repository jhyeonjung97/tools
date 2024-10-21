#!/bin/bash

declare -A rows
rows=(
    ["3d"]="Ti V Cr Mn Fe Co Ni Cu"
    ["4d"]="Zr Nb Mo Tc Ru Rh Pd"
    ["5d"]="Hf Ta W Re Os Ir Pt"
)

spins=("LS" "IS" "HS")
small_spins=("ls" "is" "hs")
specific_metals=("Fe" "Ni" "Mo" "W" "Ru" "Pd" "Pt")
adsorbates=("clean" "h" "o-o" "o-oh" "o" "oh-o" "oh-oh" "oh" "oho" "ohoh" "oo" "ooh")

# 순서를 보장하기 위해 명시적으로 행을 지정
ordered_rows=("3d" "4d" "5d")

# Loop through rows in defined order
i=0
for row in "${ordered_rows[@]}"; do
    metals=${rows[$row]}
    metal_index=1  # Reset index for each row

    for metal in $metals; do
        # Check if the metal is in specific_metals
        if [[ " ${specific_metals[@]} " =~ " ${metal} " ]]; then
            ((i++))

            for s in $(seq 0 2); do
                spin=${spins[$s]}
                small_spin=${small_spins[$s]}

                for dz in "1" "5"; do
                    base_path="/pscratch/sd/j/jiuy97/6_MNC/0_clean/${row}/${metal_index}_${metal}"
                    path=$(find "$base_path" -maxdepth 1 -type d -name "*_${spin}" 2>/dev/null)

                    # Check if path exists
                    if [ -n "$path" ]; then
                        target_dir="${path}/${dz}_"
                        echo "Processing: $target_dir"

                        if [ -d "$target_dir" ]; then
                            cd "$target_dir" || continue

                            # Execute add*.py scripts
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
                    
                    pourbaix="/pscratch/sd/j/jiuy97/6_MNC/pourbaix/${metal_index}_${metal}"

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
        ((metal_index++))  # Increment metal index
    done
done
