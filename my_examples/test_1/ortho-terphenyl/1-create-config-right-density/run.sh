#!/bin/bash

uv run --script create-initial-config.py

density=0.2

while (( $(echo "$density < 1.2" | bc -l ) ))
do
    echo "Density" $density

    sed "s/DENSITY/$density/" params-template.toml > params.toml
    particlesmc params.toml
    cp trajectories/1/lastframe.xyz inputframe.xyz

    density=$(echo "$density" | awk '{printf "%f", $1 * 1.1}')
done
