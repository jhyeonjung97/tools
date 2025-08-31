#!/bin/bash

for i in {19..27}; do
    sumo-dosplot \
        --elements Cu.d \
        --atoms Cu.${i} \
        --prefix host${i}
done

sumo-dosplot \
    --elements Sc.d \
    --prefix fermi \
    --no


sumo-dosplot \
    --elements Co.d \
    --atoms Co.27 \
    --prefix host27