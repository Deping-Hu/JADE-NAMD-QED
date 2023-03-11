#!/bin/bash
#SBATCH -p debug
#SBATCH -o output.log
#SBATCH --mem-per-cpu=1GB
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -N 1

export PATH="/home/dhu13/Program/Anaconda341/bin:$PATH"
export PATH=/home/dhu13/Program/JADE_NAMD_mu_learn/src/krr_compoment/:$PATH

1_pre.py
#2_select_train_test_all_in_one.py
#3_kernel_regression_ini.py
