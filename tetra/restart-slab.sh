#!/bin/bash

for coord in 1_Tetrahedral_WZ 2_Tetrahedral_ZB 3_SquarePlanar_TN 4_SquarePlanar_PD 5_SquarePlanar_NB 6_Octahedral_RS 7_Pyramidal_LT
do
    for row in 3d 4d 5d fm
    do
        dir=/pscratch/sd/j/jiuy97/8_V_slab/${coord}/${row}
        cd ${dir}
        if [[ -f 'DONE' ]] || [[ -f 'unmatched' ]]; then
            mkdir -p clean
            mv * clean
            if [[ -d 'clean/o' ]]; then
                mv clean/o/ old_o
                mv clean/oh/ old_oh
            elif [[ -d 'clean/o1' ]]; then
                mv clean/o1/ old_o1
                mv clean/o2/ old_o2
                mv clean/oh1 old_oh1
                mv clean/oh2 old_oh2
            fi
        else
            pwd
        fi
    done
done