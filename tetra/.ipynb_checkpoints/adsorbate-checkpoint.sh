for dir in /pscratch/sd/j/jiuy97/4_V_slab/*_*_*/*d/*_*/
do
    # Break directory path into parts
    IFS='/' read -r -a path <<< "$dir"
    
    # Extract parts of the path
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    row=${path[-2]}
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)

    # Move to the current directory, skip if it fails
    cd "$dir" || { echo "Failed to change directory to $dir"; continue; }
    pwd

    # Skip directories where 'numb' ends with 'x', 'z', or 's'
    if [[ $numb == *x ]] || [[ $numb == *z ]] || [[ $numb == *s ]]; then
        continue

    # Case where 'restart-o.json' exists (single file case)
    elif [[ -f restart-o.json ]]; then
        mkdir -p o oh
        cp ./restart-o.json o/
        cp ./restart-oh.json oh/
        cp ~/bin/tools/tetra/submit-slab.sh o/submit.sh
        cp ~/bin/tools/tetra/submit-slab.sh oh/submit.sh
        # Update jobname in both directories' submit.sh
        sed -i -e "s/jobname/$coord$row$numb/" o/submit.sh
        sed -i -e "s/jobname/$coord$row$numb/" oh/submit.sh

    # Case where 'restart-o1.json' and 'restart-o2.json' exist (multi-file case)
    elif [[ -f restart-o1.json ]] && [[ -f restart-o2.json ]]; then
        # Create directories
        mkdir -p o1 o2 oh1 oh2
        
        # Copy restart files to the appropriate directories
        cp ./restart-o1.json o1/
        cp ./restart-o2.json o2/
        cp ./restart-oh1.json oh1/
        cp ./restart-oh2.json oh2/
        
        # Loop through directories and copy submit scripts, then update jobname
        for d in o1 o2 oh1 oh2; do
            cp ~/bin/tools/tetra/submit-slab.sh "$d/submit.sh"
            sed -i -e "s/jobname/$coord$row$numb/" "$d/submit.sh"
        done

    # If no restart files are found, report an error
    else
        echo "error: missing files in $dir"
    fi
done
