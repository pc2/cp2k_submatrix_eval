#!/bin/bash

#SBATCH -p short
#SBATCH -A hpc-prf-ldft
#SBATCH -t 00:30:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=1

##SBATCH -n 80
##SBATCH -m cyclic

source ../env.sh

ulimit -s unlimited
export OMP_NUM_THREADS=40
/usr/bin/time -v mpirun -genv I_MPI_DEBUG=4 -outfile-pattern=output.log-%r-%g-%h ../exe/Linux-x86-64-intel/cp2k.psmp H2O-dft-ls.inp > out 2> err
