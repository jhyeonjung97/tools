#!/bin/bash

for i in {19..26}; do
    sumo-dosplot \
        --elements Co.d \
        --atoms Co.${i} \
        --prefix host${i}
done

sumo-dosplot \
    --elements Sc.d \
    --prefix fermi \
    --no