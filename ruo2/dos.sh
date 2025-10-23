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

~/bin/tools/ruo2/pydos-ruo2 -notot -p "12" -p "12" -p "12" -p "12" -p "12" \
-l "dxy" -l "dyz" -l "dz2" -l "dxz" -l "dx2" \
--spd "dxy" --spd "dyz" --spd "dz2" --spd "dxz" --spd "dx2" \
-o atom12_d_orbitals.png