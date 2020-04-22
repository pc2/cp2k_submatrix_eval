source env.sh
set -x
for f in half halfsingle single double;
do
  cp cuda_sign_tc.cpp tmp.cpp 
  sed -i "s/WHATTYPE/$f/g" tmp.cpp
  nvcc tmp.cpp -O3 -c -o cuda_sign_tc.o -g -Xlinker -lopenblas -lpthread -lm -ldl -lcublas -Xcompiler "-fopenmp -march=native"
  gfortran -O3 test_interface.f90 -c -o test_interface.o -fopenmp -march=native
  gfortran -L/usr/local/cuda-10.2/lib64 test_interface.o cuda_sign_tc.o  -O3 -o test_interface_$f.x -g -lopenblas -lpthread -lm -ldl -lcudart -lcublas -fopenmp -march=native
done
