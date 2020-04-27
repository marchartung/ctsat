#!/bin/bash

cd mpi
make -j
dmesg | grep cc1plus
cd ..
mkdir -p bin
cp mpi/ctsat_mpi bin/
