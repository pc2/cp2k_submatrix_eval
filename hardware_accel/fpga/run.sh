source env.sh
for i in `seq 1 20`;
do
  export OMP_NUM_THREADS=20
  numactl --membind 0 --cpubind 0 ./bin/hostf 3972 192 0 $i > 3972_192_single_0.0_${i}.log
done
