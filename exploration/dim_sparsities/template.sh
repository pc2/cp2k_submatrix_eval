#!/bin/bash

#SBATCH -p fpga
#SBATCH -A hpc-prf-ldft
#SBATCH -t 02:00:00
#SBATCH -N _NODES_
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=emul

##SBATCH -n _n_
##SBATCH -m cyclic

source ../env.sh

ulimit -s unlimited
export OMP_NUM_THREADS=40
/usr/bin/time -v mpirun -genv I_MPI_DEBUG=4 -outfile-pattern=output.log-%r-%g-%h ../exe/Linux-x86-64-intel/cp2k.psmp H2O-dft-ls.inp > out 2> err
#/usr/bin/time -v ../exe/local/cp2k.popt H2O-dft-ls.inp > out 2> err
#mpirun -genv I_MPI_PIN=1 -genv I_MPI_PIN_DOMAIN=socket -genv I_MPI_PIN_ORDER=spread -genv I_MPI_DEBUG=4 -outfile-pattern=output.log-%r-%g-%h ../exe/local/cp2k.popt H2O-dft-ls.inp > out 2> err
#/usr/bin/time -v mpirun -outfile-pattern=output.log-%r-%g-%h cp2k.popt H2O-dft-ls.inp > out 2> err
