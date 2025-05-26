#!/bin/bash

for coord in 1_Tetrahedral_WZ 2_Tetrahedral_ZB 3_SquarePlanar_TN 4_SquarePlanar_PD 5_SquarePlanar_NB 6_Octahedral_RS 7_Pyramidal_LT 8_Tetrahedral_WZ 9_Tetrahedral_ZB
do
    for row in 3d 4d 5d
    do
        for dir in /pscratch/sd/j/jiuy97/8_V_slab/${coord}/${row}/*_*/
        do
            cd ${dir}
            if [[ -f 'unmatched' ]] || [[ -f 'DONE' ]]; then
                continue
            elif [[ -f 'submit.sh' ]]; then
                python ~/bin/get_restart3
            fi
        done
    done
done