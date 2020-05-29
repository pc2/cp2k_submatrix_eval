  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <omp.h>
#include <vector>
#include <mkl.h>

#ifdef _WIN32
#include <time.h>
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "CL/opencl.h"
#include "AOCLUtils/aocl_utils.h"

  using namespace aocl_utils;
  
void randomize_array(float* array, const int size)
  {
      for (int i = 0; i < size; ++i)
      {
          array[i] = ((float)rand() / (float)RAND_MAX);
  //        array[i] = i;
      }
      printf("first three elements %f %f %f\n", array[0],array[1],array[2]);
  }
  void unit_matrix(float* array, const int dim)
  {
      for (int i = 0; i < dim; ++i)
      {
        for (int j = 0; j < dim; ++j)
        {
          if(i==j){
            array[i*dim+j] = 1.0;
          }else{
            array[i*dim+j] = 0.0;
          }
        }
      }
  }

  void randomize_array(float* array, const int size, int fake)
  {
      randomize_array(array, size);
  }


extern "C" void multiply(int dimm_np,int dimn_np,int dimk_np,float* At,float* Bt,float * Ct,int test);
  float* At_in =NULL;
  float* Bt_in =NULL;
  float* Ct_fpga =NULL;


int main(int argc, char** argv) {
      int test=0;
      int seed=atoi(argv[1]);
      srand(seed);
      printf("seed=%d\n",seed);
      //build actual input matrices
      int dimm=atoi(argv[2]);
      int dimn=atoi(argv[3]);
      int dimk=atoi(argv[4]);

      int elema=dimm*dimk;
      int elemb=dimk*dimn;
      int elemc=dimm*dimn;
      printf("dimm=%d\n",dimm);
      printf("dimn=%d\n",dimn);
      printf("dimk=%d\n",dimk);
      //row major!!!
      At_in = (float*)malloc(elema*sizeof(float));
      Bt_in = (float*)malloc(elemb*sizeof(float));
      Ct_fpga = (float*)malloc(elemc*sizeof(float));
      randomize_array(At_in, elema);
      //unit_matrix(At,dimm);
      randomize_array(Bt_in, elemb);
      
      for(int i=0;i<2;i++)
      { 
        double t0=omp_get_wtime();
        multiply(dimm,dimn,dimk,At_in,Bt_in,Ct_fpga,test); 
        double t1=omp_get_wtime();
        printf("multiply GFLOPS=%f %f \n",(double)1.0e-9*((double)2.0 * dimm * dimn * dimk)/(t1-t0),t1-t0);
      }
  }

