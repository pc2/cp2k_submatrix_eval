#!/bin/bash
#SBATCH -Ahpc-prf-ldft
#SBATCH -pbatch
#SBATCH -N16
#SBATCH -n128
#SBATCH -t0:20:00
#SBATCH --cpu-bind=sockets

LINREP=80

ENV=../../env.sh
EXE=../../../exe/Linux-x86-64-intel/cp2k.psmp

source $ENV
export OMP_NUM_THREADS=5
export OMP_STACKSIZE=10g
export SM_THREADING=HYBRID

cp H2O-dft-ls_SZV-NS.inp ${LINREP}-input.inp
sed -i "s/###FILTER###/1.0E-5/g" ${LINREP}-input.inp
sed -i "s/###LINREP###/${LINREP}/g" ${LINREP}-input.inp
srun $EXE ${LINREP}-input.inp | tee ${LINREP}-NS.out
