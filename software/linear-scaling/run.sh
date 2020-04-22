#!/bin/bash
#SBATCH -Ahpc-prf-ldft
#SBATCH -pbatch
#SBATCH -N2
#SBATCH -n80
#SBATCH -t2:00:00

ENV=../env.sh
EXE=../../exe/Linux-x86-64-intel/cp2k.psmp

source $ENV
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=10g
export SM_THREADING=HYBRID

for nrep in $(seq 1 10); do
	cp H2O-dft-ls_SZV-SM.inp input.inp
	sed -i "s/###FILTER###/1.0E-5/g" input.inp
	sed -i "s/###NREP###/${nrep}/g" input.inp
	srun $EXE input.inp | tee ${nrep}-SM.out
done
