#!/bin/bash

for temperature in 3.0 2.0 1.5 1.0
do
    mkdir -p $temperature
    sed "s/TEMPERATURE/$temperature/" params-template.toml > $temperature/params.toml
done
