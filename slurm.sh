#!/bin/bash

#SBATCH --qos=gp_bsccs
# Expected runtime. Format: DD-HH:MM:SS
#SBATCH --time=00-00:20:00
# Numero de PROCESOS
#SBATCH --ntasks=1
#SBATCH --nodes=1

# Enable or disable SMT
#SBATCH --hint=nomultithread
#SBATCH --exclusive

#SBATCH --job-name=prueba1
# Working directory
#SBATCH --chdir=/home/bsc/bsc348002/QuickAffine
#SBATCH --output=rops.out
#SBATCH --error=profiling_2.stderr


./run_multiple_tests.sh


