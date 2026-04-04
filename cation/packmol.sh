for i in {0..99}; do
    j=$(printf "%02d" $i)
    sed -i "/output/c\output cation$j.xyz" packmol.inp
    sed -i "/seed/c\seed $i" packmol.inp
    packmol < packmol.inp
done