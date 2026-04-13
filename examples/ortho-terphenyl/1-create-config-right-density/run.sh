#!/bin/bash

uv run --script create-initial-config.py

densities=(0.2 0.22 0.24 0.27 0.30 0.33 0.36 0.40 0.44 0.48 0.53 0.58 0.64 0.70 0.77 0.85 0.93 1.02 1.12 1.2)

for density in "${densities[@]}"
do
    echo "Density" $density
    sed "s/DENSITY/$density/" params-template.toml > params.toml
    particlesmc params.toml
    cp chains/1/lastframe.xyz inputframe.xyz
    echo "Done density: $density"
done

echo "Final density: $density"