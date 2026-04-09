#!/bin/bash

for temperature in 3.0 2.0 1.6 1.4 1.25 1.2 1.15 1.1 1.05 1.0
do
    cd $temperature
    particlesmc params.toml

    cd ..
    
done
