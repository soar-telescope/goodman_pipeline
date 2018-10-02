#!/bin/bash

make --directory goodman_pipeline/data/dcr-source/dcr
chmod +x goodman_pipeline/data/dcr-source/dcr/dcr
mkdir dcrbin
cp -v  goodman_pipeline/data/dcr-source/dcr/dcr dcrbin