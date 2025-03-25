for f in chg mag volume l_bond n_bond grosspop madelung ICOHPm ICOHPn ICOBIm ICOBIn ICOOPm ICOOPn pauling ion1 ion2 ion12 ion3 Natom mass density Vatom dipole Rcoval Rmetal Rvdw Tboil Tmelt Hevap Hfus Hform
do
    keep=$(echo "chg mag volume l_bond n_bond grosspop madelung ICOHPm ICOHPn ICOBIm ICOBIn ICOOPm ICOOPn pauling ion1 ion2 ion12 ion3 Natom mass density Vatom dipole Rcoval Rmetal Rvdw Tboil Tmelt Hevap Hfus Hform" | sed "s/\b$f\b//g")
    python ~/bin/tools/tetra/lr.py --Y form --X $keep --output leaveout_$f
done

for f in chg mag volume l_bond n_bond grosspop madelung ICOHPm ICOHPn ICOBIm ICOBIn ICOOPm ICOOPn pauling ion1 ion2 ion12 ion3 Natom mass density Vatom dipole Rcoval Rmetal Rvdw Tboil Tmelt Hevap Hfus Hform
do
    python ~/bin/tools/tetra/lr.py --Y form --X $f --output addon_$f
done
