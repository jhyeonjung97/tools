#!/bin/bash

for dir in /pscratch/sd/j/jiuy97/6_MNC/*_*/*d/*_*/*_*S
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    mkdir nupdown
    mv * nupdown
    mkdir 0_ 1_ 2_ 3_ 4_ 5_ 6_
    sh ~/bin/verve/spread.sh ~/bin/tools/mnc/submit.sh
    sh ~/bin/verve/spread.sh nupdown/restart.json
    sh ~/bin/verve/spread.sh nupdown/WAVECAR
    sh ~/bin/verve/jobname.sh -r $metal$spin
    sh ~/bin/verve/sub.sh -r
done