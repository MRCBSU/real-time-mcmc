#!/bin/bash

source /mnt/azfiles/apps/hpc/hpc.env
module load Anaconda3
source /mnt/azfiles/apps/hpc/conda_init_hpc.sh
conda activate R_env_rtm
module load GCC GSL intel-compilers impi OpenBLAS