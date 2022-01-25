#!/bin/bash

source /mnt/azfiles/apps/htc/htc.env
module load Anaconda3
source /mnt/azfiles/apps/htc/conda_init_htc.sh
conda activate R_env_rtm
module load GCC GSL intel-compilers impi OpenBLAS