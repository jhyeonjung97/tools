#!/bin/bash

squeue --me > ~/mystat.txt
for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*/*
do
    cd $dir
    if [[ ! -s DONE ]] && [[ ! -n $(grep $metal$spin$dz ~/mystat.txt) ]]; then
        pwd
    fi
done