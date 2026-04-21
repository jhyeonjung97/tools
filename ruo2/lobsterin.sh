#!/usr/bin/env bash

set -euo pipefail

dopants_3d=(Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge)
dopants_4d=(Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn)
dopants_5d=(La Hf Ta W Re Os Ir Pt Au Hg Tl Pb)

in_array() {
    local value="$1"
    shift
    local item
    for item in "$@"; do
        if [[ "$item" == "$value" ]]; then
            return 0
        fi
    done
    return 1
}

for dir in */; do
    base_dir="${dir%/}"
    IFS='_' read -r _ dopant <<< "$base_dir"
    target_file="${dir}lobsterin"

    if [[ ! -f "$target_file" ]]; then
        continue
    fi

    if in_array "$dopant" "${dopants_3d[@]}"; then
        sed -i -e "/basisFunctions X/c\\basisFunctions ${dopant} 3d 4s" "$target_file"
    elif in_array "$dopant" "${dopants_4d[@]}"; then
        sed -i -e "/basisFunctions X/c\\basisFunctions ${dopant} 4d 5s" "$target_file"
    elif in_array "$dopant" "${dopants_5d[@]}"; then
        sed -i -e "/basisFunctions X/c\\basisFunctions ${dopant} 5d 6s" "$target_file"
    fi

    sed -i -e "s/X/${dopant}/g" "$target_file"
done