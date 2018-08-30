#!/bin/bash

make --directory pipeline/data/dcr-source/dcr
chmod +x pipeline/data/dcr-source/dcr/dcr
mkdir dcrbin
cp -v  pipeline/data/dcr-source/dcr/dcr dcrbin