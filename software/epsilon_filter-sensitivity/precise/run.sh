#!/bin/bash
#SBATCH -Ahpc-prf-ldft
#SBATCH -pbatch
#SBATCH -N6
#SBATCH -n240
#SBATCH -t1:00:00

ENV=../../env.sh
EXE=../../../exe/Linux-x86-64-intel/cp2k.psmp

source $ENV
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=10g
export SM_THREADING=HYBRID

EPS_0=1.0               # 10^(0/4)
EPS_1=1.778279410038923 # 10^(1/4)
EPS_2=3.162277660168379 # 10^(2/4)
EPS_3=5.623413251903491 # 10^(3/4)

for EXP in -15; do
	for eps in 1.0 ; do
		cp H2O-dft-ls_SZV-NS.inp input.inp
		sed -i "s/###FILTER###/${eps}E${EXP}/g" input.inp
		srun $EXE input.inp | tee ${eps:0:3}E${EXP}-NS.out
	done
done
