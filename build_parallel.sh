#!/bin/bash

cd parallel
make r -j
cd ..

cp parallel/ctsat_parallel_release bin/
