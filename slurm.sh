#!/bin/bash

#SBATCH --qos=gp_bsccs
# Expected runtime. Format: DD-HH:MM:SS
#SBATCH --time=00-00:05:00
# Numero de PROCESOS
#SBATCH --ntasks=1
#SBATCH --nodes=1

# Enable or disbale SMT
#SBATCH --hint=nomultithread
#SBATCH --exclusive

#SBATCH --job-name=prueba1
# Working directory
#SBATCH --chdir=/home/bsc/bsc348002/QuickedAffine
#SBATCH --output=ress.out
#SBATCH --error=profiling_2.stderr


./quickedaffine_align test_datasets/real/Illumina_100s.seq res.out 1000 64 16 0 6 5 3 3 Illumina


