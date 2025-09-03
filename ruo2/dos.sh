#!/bin/bash

sumo-dosplot --elements Ru.d --format png --prefix py

fermi=$(grep "E-fermi" OUTCAR | tail -n 1 | awk '{print $3}')
echo $fermi
