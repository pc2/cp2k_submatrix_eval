#!/bin/bash
#SBATCH -Ahpc-prf-ldft
#SBATCH -pbatch
#SBATCH -N4
#SBATCH -n160
#SBATCH -t0:30:00

ENV=env.sh
EXE=../../cp2k/exe/Linux-x86-64-intel/cp2k.psmp

source $ENV
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=10g
export SM_THREADING=HYBRID

for nrep in $(seq 1 6); do
	cp H2O-dft-ls_SZV-NS.inp input.inp
	sed -i "s/###FILTER###/1.0E-12/g" input.inp
	sed -i "s/###NREP###/${nrep}/g" input.inp
	srun $EXE input.inp | tee ${nrep}-precise.out

	cp H2O-dft-ls_SZV-NS.inp input.inp
	sed -i "s/###FILTER###/1.0E-4/g" input.inp
	sed -i "s/###NREP###/${nrep}/g" input.inp
	srun $EXE input.inp | tee ${nrep}-truncated-4.out

	cp H2O-dft-ls_SZV-NS.inp input.inp
	sed -i "s/###FILTER###/1.0E-5/g" input.inp
	sed -i "s/###NREP###/${nrep}/g" input.inp
	srun $EXE input.inp | tee ${nrep}-truncated-5.out

	cp H2O-dft-ls_SZV-NS.inp input.inp
	sed -i "s/###FILTER###/1.0E-6/g" input.inp
	sed -i "s/###NREP###/${nrep}/g" input.inp
	srun $EXE input.inp | tee ${nrep}-truncated-6.out

	cp H2O-dft-ls_SZV-NS.inp input.inp
	sed -i "s/###FILTER###/1.0E-7/g" input.inp
	sed -i "s/###NREP###/${nrep}/g" input.inp
	srun $EXE input.inp | tee ${nrep}-truncated-7.out

	cp H2O-dft-ls_SZV-NS.inp input.inp
	sed -i "s/###FILTER###/1.0E-8/g" input.inp
	sed -i "s/###NREP###/${nrep}/g" input.inp
	srun $EXE input.inp | tee ${nrep}-truncated-8.out
done
