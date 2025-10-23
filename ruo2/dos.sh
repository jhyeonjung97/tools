# fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
# ~/bin/tools/ruo2/pydos-ruo2 -z $fermi -y -12 12 --notot \
# --elem='Ru' --spd='d' --lc 'C2' \
# --elem='O' --spd='p' --lc 'C1'

fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
~/bin/tools/ruo2/pydos-ruo2 -z $fermi -y -12 12 --notot \
--elem='Re' --spd='d' --lc 'C0' \
--elem='Ru' --spd='d' --lc 'C2' \
--elem='O' --spd='p' --lc 'C1'

# sumo-dosplot --format png --prefix sumo --elements Ru.d,O.p