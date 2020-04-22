#!/bin/bash

for f in single;
do
  rm ${f}.dat
  for iter in `seq 0 20`;
  do
    n=`echo "3972" | bc`
    bs=192
    conv=0.0
    e=`grep "E=" ${n}_${bs}_${f}_${conv}_${iter}.log | awk 'BEGIN { FS = " " } ; { print $3 }'`
    inv=`grep "INV=" ${n}_${bs}_${f}_${conv}_${iter}.log | awk 'BEGIN { FS = " " } ; { print $3, $4 }'`
    echo "$n $bs $conv $iter $e $inv" >> ${f}.dat
  done
done
gnuplot plot.gnu
