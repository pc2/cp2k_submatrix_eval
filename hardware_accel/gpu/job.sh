#!/bin/bash
#CCS -t 10:00:00
#CCS -N "gpusign"
#CCS -g pc2-mitarbeiter
#CCS --res=rset=1:ncpus=16:rtx2080=1,place=free:excl

export OMP_DISPLAY_AFFINITY=TRUE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=16
export OMP_STACKSIZE=1g

ulimit -s unlimited

source env.sh

bash collect_environment.sh > gpu_node.txt

for f in half halfsingle single double;
do
  for iter in `seq 0 20`;
  do
    n=`echo "3972" | bc`
    bs=192
    conv=0.0
    ./test_interface_${f}.x $n $bs $conv $iter > ${n}_${bs}_${f}_${conv}_${iter}.log
  done
done
