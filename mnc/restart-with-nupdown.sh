#!/bin/bash

mkdir nupdown
mv * nupdown
ase convert -f nupdown/CONTCAR start.traj
cp nupdown/submit.sh .
sed -i -e "s/mnc-sol/mnc-sol-$1/" submit.sh