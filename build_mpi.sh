#!/bin/bash

cd mpi
make r -j
cd ..
mkdir -p bin
cp mpi/ctsat_mpi_release bin/
