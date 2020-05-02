#!/bin/bash

cd parallel
make r -j
cd ..
mkdir -p bin
cp parallel/ctsat_parallel_release bin/
