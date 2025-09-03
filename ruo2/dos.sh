<<<<<<< HEAD
#!/bin/bash

sumo-dosplot --elements Ru.d --format png --prefix py

fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
echo $fermi
=======
sumo-dosplot --format png --prefix sumo --elements Ru.d

fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
pydos -z $fermi --elem='Ru' --spd='d'
>>>>>>> 3bd78e0a8761a4b10a1a6b332c0c51781d39307e
