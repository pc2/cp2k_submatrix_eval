export OMP_NUM_THREADS=40
module reset
module load compiler/GCC/8.3.0
module load lang/Python/3.7.0-foss-2018b

gfortran -fcheck=all -march=native -O3 -fopenmp sm.f90 kmeans.f -g -o sm.x
cat $1/fort.* > $1/blocks.dat
./sm.x $1 $2 $3
