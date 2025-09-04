fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
pydos -z $fermi --elem='Ru' --spd='d' --elem='O' --spd='p'
sumo-dosplot --format png --prefix sumo --elements Ru.d,O.p