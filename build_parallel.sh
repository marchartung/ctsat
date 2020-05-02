#!/bin/bash

cd parallel
make r -j
./ctsat_parallel_release -nthreads=1 ./sokoban-p09.sas.cr.25.cnf -cpu-lim=20
cd ..
mkdir -p bin
cp parallel/ctsat_parallel_release bin/
