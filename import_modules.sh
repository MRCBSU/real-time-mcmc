#!/bin/bash

module load gcc/gcc-9.4.0 R openmp intel_compiler OpenBLAS openmpi slurm

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/packages/OpenBLAS/0.3.17/lib64"
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/packages/openmp/12.0.1/lib"
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/packages/gsl/2.7/lib"
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/packages/gcc/9.4.0/lib64"
