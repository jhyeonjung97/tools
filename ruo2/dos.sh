<<<<<<< HEAD
fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
pydos -z $fermi --elem='Ru' --spd='d' --elem='O' --spd='p'
sumo-dosplot --format png --prefix sumo --elements Ru.d,O.p
=======
#!/bin/bash

sumo-dosplot --format png --prefix sumo --elements Ru.d

fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
pydos -z $fermi --elem='Ru' --spd='d'
>>>>>>> f81fdce7a77ac8b27f3fa988a235f5ba060a3047
