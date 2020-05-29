source env.sh
set -x
opt=""
#opt="-fstack-protector -fsanitize=address"
rm *.o

outs=""
touch host/src/multiply.o
for f in host/src/multiply.cpp ../common/src/AOCLUtils/opencl.cpp ../common/src/AOCLUtils/options.cpp host/src/main.cpp;
do
  out=`echo "$f" | awk 'BEGIN { FS = "/" } ; { print $NF }' | sed "s/cpp/o/g"`
  outs="$outs $out"
  icc  -O3 -fopenmp -g -xHost  -L/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64 -Wl,--no-as-needed -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl  -m64 -I/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/include $opt -D_FORTIFY_SOURCE=2 -Wformat -Wformat-security -fPIE -fPIC  -Ihost/inc -I../common/inc -I/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/include   $f -c -o $out  -L/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/linux64/lib -z noexecstack -Wl,-z,relro,-z,now -Wl,-Bsymbolic -pie -lOpenCL -lrt -lpthread
done

icc  -O3 -fopenmp -g -xHost -L/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64 -Wl,--no-as-needed -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl  -m64 -I/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/include -fstack-protector $opt -D_FORTIFY_SOURCE=2 -Wformat -Wformat-security -fPIE -fPIC -Ihost/inc -I../common/inc -I/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/include  $outs -L/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/linux64/lib -z noexecstack -Wl,-z,relro,-z,now -Wl,-Bsymbolic -pie -lOpenCL -lrt -lpthread -o bin/host

for f in main.f90;
do
  out=`echo "$f" | awk 'BEGIN { FS = "/" } ; { print $NF }' | sed "s/f90/o/g"`
#  outs="$outs $out"
  ifort -fcheck=all -ffree-form -ffixed-line-length-none  -O3 -fopenmp -g -xHost  -L/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64 -Wl,--no-as-needed -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl -m64 -I/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/include $opt -D_FORTIFY_SOURCE=2 -fPIE -fPIC -Ihost/inc -I../common/inc -I/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/include   $f -c -o $out  -L/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/linux64/lib -z noexecstack -Wl,-z,relro,-z,now -Wl,-Bsymbolic -pie -lOpenCL -lrt -lpthread
done

ifort -O3 -fopenmp -g -xHost  -L/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64 -Wl,--no-as-needed -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl  -m64 -I/cm/shared/opt/intel/19.04/compilers_and_libraries_2019.4.243/linux/mkl/include $opt -D_FORTIFY_SOURCE=2 -Wformat -Wformat-security -fPIE -fPIC  -Ihost/inc -I../common/inc -I/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/include  $outs -L/cm/shared/opt/intelFPGA_pro/19.2.0/hld/host/linux64/lib -z noexecstack -Wl,-z,relro,-z,now -Wl,-Bsymbolic -pie -lOpenCL -lrt -lpthread -o bin/hostf -lstdc++
