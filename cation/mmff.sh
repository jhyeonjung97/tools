#!/bin/bash

if [[ $1 =~ '-h' ]]; then
    echo 'usage: pack [NAME.inp] (lattice A, B, C)'
    exit 1
fi

name="${1%.*}"
seed=$2

a=$3
b=$4
c=$5
if [[ -z $a ]]; then
    echo 'use default lattice parameter 30 A, 30 A, 40 A...'
    a=30.
    b=30.
    c=40.
else
    if [[ $a != '*.*' ]]; then
        a=$a.0
    fi
    if [[ $b != '*.*' ]]; then
        b=$b.0
    fi
    if [[ $c != '*.*' ]]; then
        c=$c.0
    fi
fi

for i in {0..29}
do
    j=$(printf "%02d" $i)
    sed -i "/output/c\output $name$j.xyz" $name.inp
    sed -i "/seed/c\seed $i" $name.inp
    ~/bin/packmol/packmol < $name.inp
    echo "obabel $name$j.xyz -O $name$j.mol2"
    obabel $name$j.xyz -O $name$j.mol2
done

for i in {0..29}
do
    j=$(printf "%02d" $i)    
    echo "obminimize -n 100000000 -sd -c 1e-8 -ff MMFF94s $name$j.mol2 > $name$j.pdb"
    obminimize -n 100000000 -sd -c 1e-8 -ff MMFF94s $name$j.mol2 > $name$j.pdb
done

sh ~/bin/verve/mmff-result.sh
python3 ~/bin/orange/convert.py pdb traj $a $b $c

for i in {0..29}
do
      if [[ $i -lt 10 ]]; then
              dirname="0$i"
      else
              dirname="$i"
      fi
      mkdir "$dirname"
      mv "cation$dirname."* "$dirname/"
done

dir_now=$PWD
for dir in */
do
    cd $dir
    python3 ~/bin/verve/hopping.py -n 8 cation*.traj
    python3 ~/bin/orange/convert.py traj vasp $a $b $c
    python3 ~/bin/orange/convert.py vasp traj $a $b $c
    rm *.vasp
    cd $dir_now
done