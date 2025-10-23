fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
~/bin/tools/ruo2/pydos-ruo2 -z $fermi -y -12 12 --notot -p \
--elem='Ru' --spd='d' --lc 'C2' \
--elem='O' --spd='p' --lc 'C1'

# fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
# ~/bin/tools/ruo2/pydos-reruo2 -z $fermi -y -12 12 --notot \
# --elem='Re' --spd='d' --lc 'C0' \
# --elem='Ru' --spd='d' --lc 'C2' \
# --elem='O' --spd='p' --lc 'C1'

sumo-dosplot --format png --prefix sumo --orbitals Ru.d --elements Ru.d

~/bin/tools/ruo2/pydos-ruo2 --notot -p "12" -p "12" -p "12" -p "12" -p "12" \
-l "dxy" -l "dyz" -l "dz2" -l "dxz" -l "dx2" \
--spd "dxy" --spd "dyz" --spd "dz2" --spd "dxz" --spd "dx2" \
-o atom12_d_orbitals.png

~/bin/tools/ruo2/pydos-ruo2 --notot -p "13" -p "13" -p "13" -p "13" -p "13" \
-l "dxy" -l "dyz" -l "dz2" -l "dxz" -l "dx2" \
--spd "dxy" --spd "dyz" --spd "dz2" --spd "dxz" --spd "dx2" \
-o atom13_d_orbitals.png

~/bin/tools/ruo2/pydos-ruo2 --notot -p "11" -p "11" -p "11" -p "11" -p "11" \
-l "dxy" -l "dyz" -l "dz2" -l "dxz" -l "dx2" \
--spd "dxy" --spd "dyz" --spd "dz2" --spd "dxz" --spd "dx2" \
-o atom11_d_orbitals.png

~/bin/tools/ruo2/pydos-ruo2 --notot -p "0"
-l "dxy" -l "dyz" -l "dz2" -l "dxz" -l "dx2" \
--spd "dxy" --spd "dyz" --spd "dz2" --spd "dxz" --spd "dx2" \
-o atom15_d_orbitals.png

for i in 12 13 14 15; do
    for d in dxy dyz dz2 dxz dx2; do
        ~/bin/tools/ruo2/pydos-ruo2 --notot -p $i -l $d --spd $d \
        -ymin 1.5 -ymax 1.5 -o atom${i}_${d}_orbitals.png
    done
done