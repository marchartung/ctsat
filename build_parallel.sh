#!/bin/bash

cd parallel
make -j
cd ..
mkdir -p bin
cp parallel/ctsat_parallel bin/
