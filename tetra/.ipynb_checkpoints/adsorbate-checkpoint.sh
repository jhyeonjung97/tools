for dir in /pscratch/sd/j/jiuy97/4_V_slab/*_*_*/*d/*_*
do
    IFS='/' read -r -a path <<< "$dir"
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    row=${path[-2]}
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)

    cd "$dir" || { echo "Failed to change directory to $dir"; continue; }
    pwd

    # Skip directories ending with 'x', 'z', or 's'
    if [[ $numb == *x ]] || [[ $numb == *z ]] || [[ $numb == *s ]]; then
        continue

    # Case when restart-o.json exists
    elif [[ -f restart-o.json ]]; then
        mkdir -p o oh
        cp ./restart-o.json o/
        cp ./restart-oh.json oh/
        cp ~/bin/tools/tetra/submit-slab.sh o/submit.sh
        cp ~/bin/tools/tetra/submit-slab.sh oh/submit.sh
        # Update jobname in both directories' submit.sh
        sed -i -e "s/jobname/$coord$row$numb/" o/submit.sh
        sed -i -e "s/jobname/$coord$row$numb/" oh/submit.sh

    # Case when restart-o1.json and restart-o2.json exist
    elif [[ -f restart-o1.json ]] && [[ -f restart-o2.json ]]; then
        mkdir -p o1 o2 oh1 oh2
        cp ./restart-o1.json o1/
        cp ./restart-o2.json o2/
        cp ./restart-oh1.json oh1/
        cp ./restart-oh2.json oh2/
        
        # Copy the submit script to all directories in one go
        for d in o1 o2 oh1 oh2; do
            mkdir -p $d
            cp ./restart-$d.json $d/
            cp ~/bin/tools/tetra/submit-slab.sh "$d/submit.sh"
            sed -i -e "s/jobname/$coord$row$numb/" "$d/submit.sh"
        done

    # If no valid conditions are met, print error
    else
        echo "error: missing files in $dir"
    fi
done
