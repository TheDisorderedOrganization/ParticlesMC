#!/bin/bash

for temperature in 3.0 2.0 1.5 1.0
do
    cd $temperature
    particlesMC params.toml
    python3 ../plotE.py
    python3 ../plotg.py 
    cd ..
done
