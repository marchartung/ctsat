#!/bin/bash

cd mpi
make -j
cd ..
mkdir -p bin
cp mpi/ctsat_mpi bin/
