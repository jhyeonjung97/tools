#!/bin/bash

# metals=('2_Ti' '3_V' '4_Cr' '5_Mn' '6_Fe' '7_Co' '8_Ni' '9_Cu')
# spins=('1_LS' '2_IS' '2_HS' '3_HS')
# dzs=('1_' '2_' '3_' '4_' '5_' '6_')

# for metal in ${metals[@]}; do
#     for spin in ${spins[@]}; do
#         for dz in ${dzs[@]}; do
#             path="/scratch/x2755a09/3_MNC/3d/$metal/$spin/$dz"
#             if [ -d $path ]; then
#                 cd $path
#                 echo $PWD
#                 # if [ -
#                 # qsub submit.sh
#             fi
#         done
#     done
# done

# sh ~/bin/verve/spread.sh 0_/restart.json
# sh ~/bin/verve/spread.sh 0_/WAVECAR

dzs=('1_' '2_' '3_' '4_' '5_' '6_')
dir_now=$PWD
for dz in ${dzs[@]}; do
    if [ -d $dzs ]; then
        cd "$dir_now"/"$dzs"
        echo $PWD
    fi
done