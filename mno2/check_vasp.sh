#!/bin/bash

THRESHOLD=900   # vasp.out update threshold (15 min)
MIN_RUNTIME=900 # minimum runtime before checking (15 min)

for jobid in $(squeue -u $USER -t RUNNING -h -o %A); do

    runtime=$(squeue -j $jobid -h -o %M)

    runtime_sec=$(echo $runtime | awk -F: '
    NF==3 {print $1*3600 + $2*60 + $3}
    NF==2 {print $1*60 + $2}
    ')

    if [ "$runtime_sec" -lt "$MIN_RUNTIME" ]; then
        continue
    fi

    workdir=$(scontrol show job $jobid | awk -F= '/WorkDir/ {print $2}' | awk '{print $1}')

    vaspout="$workdir/vasp.out"

    if [ -f "$vaspout" ]; then

        last_update=$(stat -c %Y "$vaspout")
        current_time=$(date +%s)

        diff=$((current_time - last_update))

        if [ $diff -gt $THRESHOLD ]; then
            echo "Cancelling RUNNING job $jobid : vasp.out not updated for $diff seconds"
            scancel "$jobid"
            if cd "$workdir"; then
                python "$HOME/bin/get_restart3"
                sh "$HOME/bin/verve/sub.sh"
            else
                echo "Failed to cd to $workdir; skip get_restart3 and sub.sh"
            fi
        fi

    fi

done