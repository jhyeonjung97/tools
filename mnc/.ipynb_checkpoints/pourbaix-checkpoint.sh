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
                    path="/pscratch/sd/j/jiuy97/6_MNC/0_clean/${row}/${m}_${metal}/*_${spin}/${dz}_"
                    pourbaix="/pscratch/sd/j/jiuy97/6_MNC/pourbaix/${i}_${metal}"

                    # Check if path exists
                    if [ -d "$path" ]; then
                        cd "$path" || continue

                        # Loop through add*.py scripts in ~/bin/tools/mnc/ and execute them
                        for file in ~/bin/tools/mnc/add*.py; do
                            python "$file"
                        done
                    else
                        echo "Path does not exist: $path"
                    fi
                    
                    for adsorbate in "${adsorbates[@]}"; do
                        adsorbate_dir="${pourbaix}/${adsorbate}/${m}_${small_spin}${dz}"
                        mkdir -p "$adsorbate_dir"
                        cd "$adsorbate_dir" || continue

                        if [[ "$adsorbate" == "clean" ]]; then
                            cp "${path}/restart.json" .
                        else
                            cp "${path}/restart-${adsorbate}.json" .
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