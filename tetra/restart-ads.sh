#!/bin/bash

for row in 3d 4d; do
    for ads in o oh; do
        for dir in /pscratch/sd/j/jiuy97/8_V_slab/8_Tetrahedral_WZ/${row}/*_*/${ads}
        do
            cd $dir
            sbatch submit.sh
        done
    done
done