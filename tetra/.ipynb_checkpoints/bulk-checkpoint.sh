#!/bin/bash

metals_3d=('Ca' 'Sc' 'Ti' 'V' 'Cr' 'Mn' 'Fe' 'Co' 'Ni' 'Cu' 'Zn' 'Ga' 'Ge')
metals_4d=('Sr' 'Y' 'Zr' 'Nb' 'Mo' 'Tc' 'Ru' 'Rh' 'Pd' 'Ag' 'Cd' 'In' 'Sn')
metals_5d=('Ba' 'La' 'Hf' 'Ta' 'W' 'Re' 'Os' 'Ir' 'Pt' 'Au' 'Hg' 'Tl' 'Pb')

dir_now=$PWD

mkdir -p 3d 4d 5d

# Process 3d metals
cd 3d
for i in "${!metals_3d[@]}"; do
    if [ $i -lt 10 ]; then
        mkdir -p 0"$i"_"${metals_3d[$i]}"
    else
        mkdir -p "$i"_"${metals_3d[$i]}"
    fi
done
cd $dir_now

# Process 4d metals
cd 4d
for i in "${!metals_4d[@]}"; do
    if [ $i -lt 10 ]; then
        mkdir -p 0"$i"_"${metals_4d[$i]}"
    else
        mkdir -p "$i"_"${metals_4d[$i]}"
    fi
done
cd $dir_now

# Process 5d metals
cd 5d
for i in "${!metals_5d[@]}"; do
    if [ $i -lt 10 ]; then
        mkdir -p 0"$i"_"${metals_5d[$i]}"
    else
        mkdir -p "$i"_"${metals_5d[$i]}"
    fi
done
cd $dir_now

cp $1 POSCAR
sed -i -e "s/$2/XX/" POSCAR
coord=$(basename $PWD | cut -d'_' -f3)
sed -i -e "s/XX/$coord/" submit.sh
sh ~/bin/verve/spread.sh -rr POSCAR submit.sh

for dir in *d/*_*/; do
    cd $dir
    metal=$(basename $PWD | cut -d'_' -f2)
    sed -i -e "s/XX/$metal/" POSCAR
    ase convert -f POSCAR start.traj
    cd $dir_now
done

for dir in *d/; do
    cd $dir
    clean_dir=$(basename $dir)
    sh ~/bin/verve/jobname.sh -rc "$coord$clean_dir"
    cd $dir_now
done