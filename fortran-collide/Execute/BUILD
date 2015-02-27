#!/bin/bash

make clobber
make
make clean

mpirun -np $1 ./COLLIDE
