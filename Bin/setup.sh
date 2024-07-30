#!/usr/bin/env bash

module purge

module load intel/2024.0
module load gcc/11
module load openmpi/4.1

export OMP_STACKSIZE=16000
export OMP_NUM_THREADS=1

#cmake  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran ../
#make -j