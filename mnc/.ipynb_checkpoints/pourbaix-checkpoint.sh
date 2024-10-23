#!/bin/bash

declare -A rows
rows=(
    ["3d"]="Ti V Cr Mn Fe Co Ni Cu"
    ["4d"]="Zr Nb Mo Tc Ru Rh Pd"
    ["5d"]="Hf Ta W Re Os Ir Pt"
)

spins=("LS" "IS" "HS")
small_spins=("ls" "is" "hs")
# specific_metals=("Fe" "Co" "Mo" "W" "Ru" "Pd" "Pt")
specific_metals=("Co")
adsorbates=("clean" "h" "o-o" "o-oh" "o" "oh-o" "oh-oh" "oh" "oho" "ohoh" "oo" "ooh")
# adsorbates=("ooh-oh" "oh-ooh" "ooh-o" "o-ooh" "oo-oh" "oh-oo" "oo-o" "o-oo")

# 순서를 보장하기 위해 명시적으로 행을 지정
ordered_rows=("3d" "4d" "5d")

# Loop through rows in defined order
index=2
for row in "${ordered_rows[@]}"; do
    metals=${rows[$row]}
    metal_index=2  # Reset index for each row

    for metal in $metals; do
        # Check if the metal is in specific_metals
        if [[ " ${specific_metals[@]} " =~ " ${metal} " ]]; then
            # ((index++))

            for s in $(seq 0 2); do
                spin=${spins[$s]}
                small_spin=${small_spins[$s]}

                for dz in "1" "5"; do
                    base_path="/pscratch/sd/j/jiuy97/6_MNC/0_clean/${row}/${metal_index}_${metal}"
                    path=$(find "$base_path" -maxdepth 1 -type d -name "*_${spin}" 2>/dev/null)
                    pourbaix="/pscratch/sd/j/jiuy97/6_MNC/pourbaix/${index}_${metal}"

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
                        echo "Processing: $target_dir"
                        spin='IS'
                    fi
                    
                    # Create pourbaix directories
                    for adsorbate in "${adsorbates[@]}"; do
                        adsorbate_dir="${pourbaix}/${adsorbate}/${spin}${dz}"
                        echo "$adsorbate_dir"
                        mkdir -p "$adsorbate_dir"
                        cd "$adsorbate_dir" || continue

                        if [[ "$adsorbate" == "clean" ]]; then
                            restart_file="${target_dir}/restart.json"
                        else
                            restart_file="${target_dir}/restart-${adsorbate}.json"
                        fi

                        if [ -f "$restart_file" ]; then
                            cp "$restart_file" ./restart.json
                            sed -i -e '/constraints/d' 'restart.json'
                            ase convert -f 'restart.json' 'POSCAR'
                            ase convert -f 'POSCAR' 'start.traj'
                            rm 'restart.json'
                        else
                            echo "File not found: $restart_file"
                            continue
                        fi
                        
                        cp ~/bin/tools/mnc/submit.sh .
                        sed -i -e "s/jobname/${index}${metal}${spin}${adsorbate}${dz}/" submit.sh
                        sed -i -e "s/mnc-sol/mnc-sol-${small_spin}-nupdown/" submit.sh
                    done
                done
            done
        fi
        ((metal_index++))  # Increment metal index
    done
done
