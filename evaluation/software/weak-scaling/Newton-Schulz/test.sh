#!/bin/bash
#SBATCH -Ahpc-prf-ldft
#SBATCH -pshort
#SBATCH -N2
#SBATCH -n16
#SBATCH -t0:1:00
#SBATCH --cpu-bind=sockets

hostname
