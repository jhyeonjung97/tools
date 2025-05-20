#!/bin/bash

for coord in 1_Tetrahedral_WZ 2_Tetrahedral_ZB 3_SquarePlanar_TN 4_SquarePlanar_PD 5_SquarePlanar_NB 6_Octahedral_RS 7_Pyramidal_LT 8_Tetrahedral_WZ 9_Tetrahedral_ZB
do
    for row in 3d 4d 5d fm
    do
        dir=/pscratch/sd/j/jiuy97/8_V_slab/${coord}/${row}
        cd ${dir}
        for subdir in *_*/
        do
            cd ${dir}/${subdir}
            if [[ -f submit.sh ]]; then
                if [[ ! -f 'DONE' ]]; then
                    if [[ ! -f 'unmatched' ]]; then
                        continue
                    fi
                fi
            fi
            mkdir -p static
            cd static
            cp ../POSCAR .
            cp ../POTCAR .
            cp ../KPOINTS .
            cp ../INCAR .
            cp ../submit.sh .
            cd ..
        done
    done
done