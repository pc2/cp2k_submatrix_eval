#!/bin/bash
#SBATCH -t 1-0
#SBATCH -n 1
#SBATCH -J mm_16_16_8_32_32_16_1_1_7
#SBATCH -A hpc-prf-ldft
#SBATCH -p long
source env.sh
export OMP_NUM_THREADS=40

mkdir syn
aoc -seed=7 device/matrix_mult.cl -o syn/matrix_mult.aocx -global-ring -duplicate-ring -fp-relaxed -DPE_ROWS=16 -DPE_COLS=16 -DDOT_PROD_VECTOR_SIZE=8 -DROWS_INTERLEAVED=32 -DCOLUMNS_INTERLEAVED=32 -DSCALING_FACTOR=16 -board=p520_hpc_sg280l > syn/out 2> syn/err
