#!/bin/bash
#SBATCH -Ahpc-prf-ldft
#SBATCH -pbatch
#SBATCH -N4
#SBATCH -n160
#SBATCH -t0:10:00

ENV=../env.sh
EXE=../../exe/Linux-x86-64-intel/cp2k.psmp

source $ENV
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=10g
export SM_THREADING=HYBRID

cp H2O-dft-ls_SZV-SM.inp input.inp
sed -i "s/###FILTER###/1.0E-5/g" input.inp
sed -i "s/###NREP###/7/g" input.inp
srun $EXE input.inp | tee 4-SM.out
