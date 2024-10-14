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

    # Skip directories where 'numb' ends with 'x', 'z', or 's'
    if [[ $numb == *x ]] || [[ $numb == *z ]] || [[ $numb == *s ]]; then
        continue
    elif [[ $coord == 'AQ' ]] || [[ $coord == 'AU' ]]; then
        continue
    else
        pwd
    fi

    # Case where 'restart-o.json' exists (single file case)
    if [[ -f restart-o.json ]]; then
        mkdir -p o oh
        cp ./restart-o.json o/restart.json
        cp ./restart-oh.json oh/restart.json
        cp ~/bin/tools/tetra/submit-slab.sh o/submit.sh
        cp ~/bin/tools/tetra/submit-slab.sh oh/submit.sh
        # Update jobname in both directories' submit.sh
        sed -i -e "s/jobname/${coord}${row}${numb}o/" o/submit.sh
        sed -i -e "s/jobname/${coord}${row}${numb}oh/" oh/submit.sh

    # Case where 'restart-o1.json' and 'restart-o2.json' exist (multi-file case)
    elif [[ -f restart-o1.json ]] && [[ -f restart-o2.json ]]; then
        # Create directories
        mkdir -p o1 o2 oh1 oh2
        
        # Copy restart files to the appropriate directories
        cp ./restart-o1.json o1/restart.json
        cp ./restart-o2.json o2/restart.json
        cp ./restart-oh1.json oh1/restart.json
        cp ./restart-oh2.json oh2/restart.json
        
        # Loop through directories and copy submit scripts, then update jobname
        for d in o1 o2 oh1 oh2; do
            cp ~/bin/tools/tetra/submit-slab.sh "$d/submit.sh"
            sed -i -e "s/jobname/${coord}${row}${numb}${d}/" "$d/submit.sh"
        done

    # If no restart files are found, report an error
    else
        echo "error: missing files in $dir"
    fi
done
