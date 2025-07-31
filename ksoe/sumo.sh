#!/bin/bash

sumo-dosplot \
    --elements Co.d \
    --atoms Co.19.20.22.23.25.26 \
    --prefix surround

sumo-dosplot \
    --elements Co.d \
    --atoms Co.19.20.21.22.23.24.25.26 \
    --prefix top

sumo-dosplot \
    --elements Co.d \
    --prefix host