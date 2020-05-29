  // Copyright (C) 2013-2019 Altera Corporation, San Jose, California, USA. All rights reserved.
  // Permission is hereby granted, free of charge, to any person obtaining a copy of this
  // software and associated documentation files (the "Software"), to deal in the Software
  // without restriction, including without limitation the rights to use, copy, modify, merge,
  // publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to
  // whom the Software is furnished to do so, subject to the following conditions:
  // The above copyright notice and this permission notice shall be included in all copies or
  // substantial portions of the Software.
  // 
  // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  // EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
  // OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  // NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
  // HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
  // WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  // FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
  // OTHER DEALINGS IN THE SOFTWARE.
  // 
  // This agreement shall be governed in all respects by the laws of the State of California and
  // by the laws of the United States of America.

  ///////////////////////////////////////////////////////////////////////////////////
  // This host program executes a matrix multiplication kernel to perform:
  //  C = A * B
  // where A is a N x K matrix, B is a K x M matrix and C is a N x M matrix.
  // All dimensions must be a multiple of BLOCK_SIZE, which affects the
  // underlying kernel.
  //
  // Verification is performed against the same computation on the host CPU.
  ///////////////////////////////////////////////////////////////////////////////////

  // System includes
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



  // Constants
#define HOST

#define ACL_ALIGNMENT 64

#ifdef _WIN32
  void* acl_aligned_malloc (size_t size) {
      return _aligned_malloc (size, ACL_ALIGNMENT);
  }
  void acl_aligned_free (void *ptr) {
      _aligned_free (ptr);
  }
#else
  void* acl_aligned_malloc (size_t size) {
      void *result = NULL;
      if (posix_memalign(&result, ACL_ALIGNMENT, size) != 0)
          printf("acl_aligned_malloc() failed.\n");
      return result;
  }
  void acl_aligned_free (void *ptr) {
      free (ptr);
  }
#endif // LINUX


#define NUM_NON_AUTORUN_KERNELS 6

#include "config.h"
#include "matrixMul_large.h"

#define AOCX_FILE "bin/matrix_mult.aocx"

#define NUM_KERNELS   ( PE_ROWS*PE_COLS + PE_ROWS + PE_COLS + (PE_COLS-1) + (1 + 1 + 1) )
#define NUM_QUEUES    NUM_KERNELS

#define   KID_FEED_MAT_A1    0
#define   KID_FEED_MAT_B1    1
#define   KID_DRAIN_MAT_C1   2
#define   KID_FEED_MAT_A2    3
#define   KID_FEED_MAT_B2    4
#define   KID_DRAIN_MAT_C2   5
#define   KID_HOST_COMM     6

  // only create non-autorun kernels for the HW run
#define NUM_KERNELS_TO_CREATE   NUM_NON_AUTORUN_KERNELS
#define NUM_QUEUES_TO_CREATE    NUM_NON_AUTORUN_KERNELS
#define NUM_QUEUES_TO_FINISH    NUM_NON_AUTORUN_KERNELS

#ifndef min
#define min(a,b) ((a<b) ? (a) : (b))
#endif

  // OpenCL runtime configuration
  cl_kernel kernel[NUM_KERNELS_TO_CREATE];
  cl_command_queue cmdQueue[NUM_QUEUES_TO_CREATE+1]; // extra queue for host communication
  cl_event kernel_exec_event[NUM_QUEUES];

  std::vector<cl_mem> d_matrix_mul_outputC                   ;
  std::vector<cl_mem> d_matrix_mul_inputA                   ;
  std::vector<cl_mem> d_matrix_mul_inputB                  ; 

  int fpga_init_done=0;
  cl_program program                            = NULL;
  cl_context context                            = NULL;

  cl_platform_id platform                       = NULL;
  cl_device_id* devices                         = NULL;

  float* Ct_chunked =NULL;
  float* Ct_ref =NULL;
      
  unsigned int num_elem_A;
  unsigned int num_elem_B;
  unsigned int num_elem_C;
      size_t globalWorkSize[1];
      size_t localWorkSize[1];

  float* matrix_mul_inputA1                      = NULL;
  float* matrix_mul_inputA2                      = NULL;
  float* matrix_mul_inputB1                      = NULL; // transposed
  float* matrix_mul_inputB2                      = NULL; // transposed
  float* matrix_mul_outputC1                     = NULL;
  float* matrix_mul_outputC2                     = NULL;

  float* matrix_mul_inputB_transposed           = NULL; // non-transposed

  // The design can either generate input data on-chip (running with -i option)
  // or load input data from DDR.
  // With only 1 DDR bank, the bandwidth becomes saturated when loading the input data

  // Use real data
  // By flipping this to 1, tell the kernels to inject dummy data instead
  unsigned char disableA = 0;
  unsigned char disableB = 0;

  const char *kernel_name[] = {
      "loadA1",
      "loadB1",
      "store1",
      "loadA2",
      "loadB2",
      "store2",
  };

  // Function prototypes
  void cleanup_host_side_resources();
  void cleanup();

  // Check the status returned by the OpenCL API functions
#define CHECK(status)                               \
  if (status != CL_SUCCESS)                       \
  {                                   \
      fprintf(stderr, "error %d in line %d.\n", status, __LINE__);    \
      exit(1);                            \
  }                                   \

  // Check the status returned by the OpenCL API functions, don't exit on error
#define CHECK_NO_EXIT(status)                               \
  if (status != CL_SUCCESS)                       \
  {                                   \
      fprintf(stderr, "error %d in line %d.\n", status, __LINE__);    \
  }                                   \


  void printkerneltype()
  {
      printf("\n===== Host-CPU checking the systolic array matrix multiplication parameters ======\n\n");

      printf("HA: \t\t%d\n", HA);
      printf("WA: \t\t%d\n\n", WA);

      printf("HB: \t\t%d\n", HB);
      printf("WB: \t\t%d\n\n", WB);

      printf("HC: \t\t%d\n", HC);
      printf("WC: \t\t%d\n\n", WC);

      printf("PE_ROWS: \t\t%d\n", PE_ROWS);
      if (PE_ROWS<1) {
          printf("--->ERROR, PE_ROWS must be larger than 0\n");
      }
      printf("PE_COLS: \t\t%d\n", PE_COLS);
      if (PE_COLS<1l) {
          printf("--->ERROR, PE_COLS must be larger than 0\n");
          exit(1);
      }
      printf("DOT_PROD_VECTOR_SIZE: \t\t%d\n", DOT_PROD_VECTOR_SIZE);
      if (DOT_PROD_VECTOR_SIZE!=4 &&
          DOT_PROD_VECTOR_SIZE!=8 &&
          DOT_PROD_VECTOR_SIZE!=16) {
          printf("Illegal DOT_PROD_VECTOR_SIZE, supported: 4,8,16\n");
          exit(1);
      }
      printf("\n");

      printf("ACCUM_SHIFT_REG_SIZE: \t\t%d\n", ACCUM_SHIFT_REG_SIZE);

      printf("\n");

      printf("ROWS_INTERLEAVED: \t\t%d\n",    ROWS_INTERLEAVED);
      printf("COLUMNS_INTERLEAVED: \t\t%d\n", COLUMNS_INTERLEAVED);
      printf("\n");

      printf("MAT_A_BLOCK_HEIGHT: \t\t%d\n",  MAT_A_BLOCK_HEIGHT);
      printf("MAT_A_BLOCK_WIDTH: \t\t%d\n",   MAT_A_BLOCK_WIDTH);
      printf("MAT_A_BLOCK_SIZE: \t\t%d\n",    MAT_A_BLOCK_SIZE);

      printf("MAT_A_BLOCK_NUM_VECTORS: \t\t%d\n",   MAT_A_BLOCK_NUM_VECTORS);
      if (MAT_A_BLOCK_SIZE % DOT_PROD_VECTOR_SIZE) {
          printf("MAT_A_BLOCK_SIZE must be a multiple of DOT_PROD_VECTOR_SIZE\b");
      }
      printf("MAT_A_NUM_BLOCKS_IN_ROW: \t\t%d\n",   MAT_A_NUM_BLOCKS_IN_ROW);
      if (WA % MAT_A_BLOCK_WIDTH) {
          printf("WA must be a multiple of MAT_A_BLOCK_WIDTH\n");
      }
      printf("MAT_A_NUM_BLOCKS_IN_COL: \t\t%d\n",   MAT_A_NUM_BLOCKS_IN_COL);
      if (HA % MAT_A_BLOCK_HEIGHT) {
          printf("HA must be a multiple of MAT_A_BLOCK_HEIGHT\n");
      }
      printf("MAT_A_NUM_VECTORS_IN_ROW_OF_BLOCKS: \t\t%d\n",   MAT_A_NUM_VECTORS_IN_ROW_OF_BLOCKS);
      printf("\n");

      printf("MAT_B_BLOCK_HEIGHT: \t\t%d\n",  MAT_B_BLOCK_HEIGHT);
      printf("MAT_B_BLOCK_WIDTH: \t\t%d\n",   MAT_B_BLOCK_WIDTH);
      printf("MAT_B_BLOCK_SIZE: \t\t%d\n",    MAT_B_BLOCK_SIZE);

      printf("MAT_B_BLOCK_NUM_VECTORS: \t\t%d\n",   MAT_B_BLOCK_NUM_VECTORS);
      if (MAT_B_BLOCK_SIZE % DOT_PROD_VECTOR_SIZE) {
          printf("MAT_B_BLOCK_SIZE must be a multiple of DOT_PROD_VECTOR_SIZE\b");
      }
      printf("MAT_B_NUM_BLOCKS_IN_ROW: \t\t%d\n",   MAT_B_NUM_BLOCKS_IN_ROW);
      if (WB % MAT_B_BLOCK_WIDTH) {
          printf("WB must be a multiple of MAT_B_BLOCK_WIDTH\n");
      }
      printf("MAT_B_NUM_BLOCKS_IN_COL: \t\t%d\n",   MAT_B_NUM_BLOCKS_IN_COL);
      if (HB % MAT_B_BLOCK_HEIGHT) {
          printf("HB must be a multiple of MAT_B_BLOCK_HEIGHT\n");
      }
      printf("MAT_B_NUM_VECTORS_IN_COL_OF_BLOCKS: \t\t%d\n",  MAT_B_NUM_VECTORS_IN_COL_OF_BLOCKS);
      printf("MAT_B_NUM_VECTORS_IN_MATRIX: \t\t%d\n",         MAT_B_NUM_VECTORS_IN_MATRIX);
      printf("\n");

      printf("MAT_C_BLOCK_HEIGHT: \t\t%d\n",  MAT_C_BLOCK_HEIGHT);
      printf("MAT_C_BLOCK_WIDTH: \t\t%d\n",   MAT_C_BLOCK_WIDTH);


  }

  bool compare_L2_norm(const float* ref_array, const float* output_array, const unsigned int size, const float epsilon)
  {
      // compute the L^2-Norm of the difference between the output array and reference array
      // and compare it against the L^2-Norm of the reference.
      float diff = 0.0f;
      float ref = 0.0f;
      for (unsigned int i = 0; i < size; ++i) {

          const float o = output_array[i];
          const float r = ref_array[i];
          const float d = o - r;
          diff += d * d;
          ref += r * r;
      }

      const float diff_l2norm = sqrtf(diff);
      const float ref_l2norm = sqrtf(ref);
      const float error = diff_l2norm / ref_l2norm;
      const bool pass = error < epsilon;
      printf("%f < %f\n", error, epsilon);
      printf("diff %f %f\n", ref_array[0], output_array[0]);

      return pass;
  }

  void transpose_matrix( float * B, int hB, int wB)
  {
    double t1=omp_get_wtime();
      const char row_major = 'R';
      const char transpose = 'T';
      const float alpha = 1.0f;
      mkl_simatcopy (row_major, transpose, hB, wB, alpha, B, wB, hB);
      double t2=omp_get_wtime();
      //printf("transpose_matrix inplace %f\n",t2-t1);
  }

  void transpose_matrix( float * B_orig, float * B_transposed, int hB, int wB)
  {
    double t1=omp_get_wtime();
    if(0==1){
      double t1=omp_get_wtime();
      for(int i=0; i < wB; ++i) {
          for(int j=0; j < hB; ++j) {
              B_transposed[i*hB + j] = B_orig[j*wB + i];
          }
      }
      double t2=omp_get_wtime();
      //printf("transpose_matrix %f\n",t2-t1);
    }else{
      const char row_major = 'R';
      const char transpose = 'T';
      const float alpha = 1.0f;
      mkl_somatcopy (row_major, transpose, hB, wB, alpha, B_orig, wB,B_transposed, hB);
    }
      double t2=omp_get_wtime();
      //printf("transpose_matrix oplace %f\n",t2-t1);
  }

  void block_wise_reformat_matrix( float * A_orig, float * A_block_wise, int mat_height, int mat_width, int block_height, int block_width)
  {
      double t1=omp_get_wtime();
#pragma omp parallel for schedule(static)
      for(int i=0; i < mat_height/block_height; i++) {
          for(int j=0; j < mat_width/block_width; j++) {
              for(int k=0; k < block_height; k++) {
                  for(int l=0; l < block_width; l++) {
                      int word_id=l+block_width*(k+block_height*(j+mat_width/block_width*i));
                      A_block_wise[word_id] = A_orig[(i*block_height+k)*mat_width + (j*block_width+l)];
                  }
              }
          }
      }
      double t2=omp_get_wtime();
      //printf("block_wise_reformat_matrix %f\n",t2-t1);
  }

  double compute_kernel_execution_time(cl_event &event, double &start_d, double &end_d)
  {
      cl_ulong start, end;

      clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END,      sizeof(cl_ulong), &end,     NULL);
      clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START,    sizeof(cl_ulong), &start,   NULL);

      start_d = (double)1.0e-9 * start;
      end_d   = (double)1.0e-9 * end;

      return    (double)1.0e-9 * (end - start); // nanoseconds to seconds
  }

  int get_block_id(std::vector<std::vector<int> > blocks,int a,int b)
  {
    for(int i=0;i<blocks.size();i++){
      if(blocks[i][0]==a && blocks[i][1]==b){
        return i;
      }
    }
    return -1;
  }

  void get_chunk_from_matrix(int dimr,int dimc,int h,int w,std::vector<int> blocks,float* A,float* matrix)
  {
    double t1=omp_get_wtime();
    int Ai1=blocks[2];
    int Aj1=blocks[4];
#pragma omp parallel for schedule(static)
    for(int i=0;i<h;i++){
      for(int j=0;j<w;j++){
       if(Ai1+i<dimr && Aj1+j<dimc){
          matrix[i*w+j]=A[(Ai1+i)+(Aj1+j)*dimr];
        }else{
          matrix[i*w+j]=0.0;
        }
      }
    }
    double t2=omp_get_wtime();
    //printf("get_chunk_from_matrix %f %d %d, %d %d, %d %d\n",t2-t1,dimr,dimc,Ai1,Aj1,h,w);
  }

  void block_wise_reformat_chunk(int dimr,int dimc,float*A,int h,int w,std::vector<int> blocks,float * A_block_wise, int block_height, int block_width)
  {
      double t1=omp_get_wtime();
      int Ai1=blocks[2];
      int Aj1=blocks[4];
#pragma omp parallel for schedule(static)
      for(int j=0; j < w/block_width; j++) {
        for(int i=0; i < h/block_height; i++) {
            for(int l=0; l < block_width; l++) {
              for(int k=0; k < block_height; k++) {
                      int word_id=l+block_width*(k+block_height*(j+w/block_width*i));
                      if(Ai1+i*block_height+k<dimr && Aj1+j*block_width+l<dimc){
                        A_block_wise[word_id] =A[(Ai1+i*block_height+k)+dimr*(Aj1+j*block_width+l)];
                      }else{
                        A_block_wise[word_id] =0.0;
                      }
                  }
              }
          }
      }
      double t2=omp_get_wtime();
//      printf("block_wise_reformat_chunk %f\n",t2-t1);
  }
  void set_matrix_from_chunk_add(int dimr,int dimc,int h,int w,std::vector<int> blocks,float* A,float* matrix,int tadd)
  {
    int Ai1=blocks[2];
    int Aj1=blocks[4];
    int h2=min(h,dimr-Ai1);
    int w2=min(w,dimc-Aj1);
    printf("set_matrix_from_chunk_add %d %d, %d %d, %d %d,  %d\n",dimr,dimc,Ai1,Aj1,h2,w2,tadd);
    if(tadd==1){
      #pragma omp parallel for schedule(static)
      for(int j=0;j<w2;j++){
        for(int i=0;i<h2;i++){
              A[(Ai1+i)+(Aj1+j)*dimr]+=matrix[i*w+j];
            }
          }
    }else{
      #pragma omp parallel for schedule(static)
      for(int j=0;j<w2;j++){
        for(int i=0;i<h2;i++){
              A[(Ai1+i)+(Aj1+j)*dimr]=matrix[i*w+j];
            }
          }
    }
  }
  void set_matrix_from_chunk_reorder_add(int dimr,int dimc,std::vector<int> blocks,float* A,float* C0,int mat_height, int mat_width,int block_height, int block_width,int num_sys_arr_columns,int tadd)
  {
    int Ai1=blocks[2];
    int Aj1=blocks[4];
    
    double t1=omp_get_wtime();
    int num_elems = HC*WC;
    int column_interleaving = MAT_C_BLOCK_WIDTH / PE_COLS;
    
    if(tadd==1){
      #pragma omp parallel for schedule(static)
      for(int i=0; i < num_elems/MAT_C_BLOCK_WIDTH; i++) {
          for(int j=0; j < column_interleaving; j++) {
        for(int k=0; k < PE_COLS ; k++) {
                  int word_id=i*MAT_C_BLOCK_WIDTH+j+k*column_interleaving;
                  int word_id0=k+PE_COLS*(j+column_interleaving*i);

                  int l2=word_id%MAT_C_BLOCK_WIDTH;
                  int k2=((word_id)/MAT_C_BLOCK_WIDTH)%MAT_C_BLOCK_HEIGHT;
                  int j2=(((word_id)/MAT_C_BLOCK_WIDTH)/MAT_C_BLOCK_HEIGHT)%(WC/MAT_C_BLOCK_WIDTH);
                  int i2=(((word_id)/MAT_C_BLOCK_WIDTH)/MAT_C_BLOCK_HEIGHT)/(WC/MAT_C_BLOCK_WIDTH);
                  if((Ai1+(i2*MAT_C_BLOCK_HEIGHT+k2))<dimr && (Aj1+(j2*MAT_C_BLOCK_WIDTH+l2))<dimc){           
                    A[(Ai1+(i2*MAT_C_BLOCK_HEIGHT+k2))+dimr*(Aj1+(j2*MAT_C_BLOCK_WIDTH+l2))]+=C0[word_id0];
                  }
              }
          }
      }
    }else{
      #pragma omp parallel for schedule(static)
          for(int i=0; i < num_elems/MAT_C_BLOCK_WIDTH; i++) {
      for(int j=0; j < column_interleaving; j++) {
        for(int k=0; k < PE_COLS ; k++) {
                  int word_id=i*MAT_C_BLOCK_WIDTH+j+k*column_interleaving;
                  int word_id0=k+PE_COLS*(j+column_interleaving*i);

                  int l2=word_id%MAT_C_BLOCK_WIDTH;
                  int k2=((word_id)/MAT_C_BLOCK_WIDTH)%MAT_C_BLOCK_HEIGHT;
                  int j2=(((word_id)/MAT_C_BLOCK_WIDTH)/MAT_C_BLOCK_HEIGHT)%(WC/MAT_C_BLOCK_WIDTH);
                  int i2=(((word_id)/MAT_C_BLOCK_WIDTH)/MAT_C_BLOCK_HEIGHT)/(WC/MAT_C_BLOCK_WIDTH);
                  
                  if((Ai1+(i2*MAT_C_BLOCK_HEIGHT+k2))<dimr && (Aj1+(j2*MAT_C_BLOCK_WIDTH+l2))<dimc){           
                    A[(Ai1+(i2*MAT_C_BLOCK_HEIGHT+k2))+dimr*(Aj1+(j2*MAT_C_BLOCK_WIDTH+l2))]=C0[word_id0];
                  }
              }
          }
       }

    }
    
  }
  void initfpga(){
      if(fpga_init_done==1){
        return;
      }
      printkerneltype();
//exit(1);
      FILE *f_out = stdout;
      cl_int status;

      //----------------------------------------------
      // Get the OpenCL platform
      //----------------------------------------------
      platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
      if(platform == NULL) {
          printf("ERROR: Unable to find Intel(R) FPGA OpenCL platform\n");
          cleanup_host_side_resources();
          exit(1);
      }

      //----------------------------------------------
      // Discover and initialize the devices
      //----------------------------------------------

      cl_uint numDevices = 0;

      // Device info
      char buffer[4096];
      unsigned int buf_uint;
      int device_found = 0;

      printf("Initializing IDs\n");
      status = clGetDeviceIDs(platform,
                      CL_DEVICE_TYPE_ALL,
                      0,
                      NULL,
                      &numDevices);

      if(status == CL_SUCCESS){
          clGetPlatformInfo(platform,
                          CL_PLATFORM_VENDOR,
                          4096,
                          buffer,
                          NULL);

          if(strstr(buffer, "Intel(R)") != NULL){
                  device_found = 1;
          }
          printf("%s\n", buffer);

          if(device_found){
              // Allocate enough space for each device
              devices = (cl_device_id*)
              acl_aligned_malloc (numDevices * sizeof(cl_device_id));

              // Fill in devices with clGetDeviceIDs()
              status = clGetDeviceIDs(platform,
                              CL_DEVICE_TYPE_ALL,
                              numDevices,
                              devices,
                              NULL);
          }
      }

      if(!device_found) {
          printf("failed to find a OpenCL device\n");
          exit(1);
      }

      for (int i = 0; i < numDevices; i++) {
          clGetDeviceInfo(devices[i],
                          CL_DEVICE_NAME,
                          4096,
                          buffer,
                          NULL);
          fprintf(f_out, "\nDevice Name: %s\n", buffer);

          clGetDeviceInfo(devices[i],
                          CL_DEVICE_VENDOR,
                          4096,
                          buffer,
                          NULL);
          fprintf(f_out, "Device Vendor: %s\n", buffer);

          clGetDeviceInfo(devices[i],
                          CL_DEVICE_MAX_COMPUTE_UNITS,
                          sizeof(buf_uint),
                          &buf_uint,
                          NULL);
          fprintf(f_out, "Device Computing Units: %u\n", buf_uint);

          clGetDeviceInfo(devices[i],
                          CL_DEVICE_GLOBAL_MEM_SIZE,
                          sizeof(unsigned long),
                          &buffer,
                          NULL);
          fprintf(f_out, "Global Memory Size: %lu\n", *((unsigned long*)buffer));

          clGetDeviceInfo(devices[i],
                          CL_DEVICE_MAX_MEM_ALLOC_SIZE,
                          sizeof(unsigned long),
                          &buffer,
                          NULL);
          fprintf(f_out, "Global Memory Allocation Size: %lu\n\n", *((unsigned long*)buffer));
      }



      //----------------------------------------------
      // Create a context
      //----------------------------------------------

      printf("\n===== Host-CPU setting up the OpenCL command queues ======\n\n");

      // Create a context using clCreateContext() and associate it with the device

      context = clCreateContext(
                      NULL,
                      1,
                      devices,
                      NULL,
                      NULL,
                      &status); CHECK(status);

      //----------------------------------------------
      // Create the program from binaries
      //----------------------------------------------
      printf("\n===== Host-CPU setting up OpenCL program and kernels ======\n\n");

      size_t binary_length;
      const unsigned char *binary;

      printf("\nAOCX file: %s\n\n", AOCX_FILE);
      // create the program using binary already compiled offline using aoc (i.e. the .aocx file)
      FILE *fp = fopen(AOCX_FILE, "rb");

      if (fp == NULL) {
          printf("Failed to open the AOCX file (fopen).\n");
          cleanup_host_side_resources();
          exit(1);
      }

      fseek(fp, 0, SEEK_END);
      long ftell_sz = ftell(fp);
      if (ftell_sz < 0) {
          printf("ftell returns a negative value.\n");
          fclose(fp);
          cleanup_host_side_resources();
          exit(1);
      }
      else {
          binary_length = ftell_sz;
      }
      binary = (unsigned char*) malloc(sizeof(unsigned char) * binary_length);
      assert(binary && "Malloc failed");
      rewind(fp);

      size_t fread_sz = fread((void*)binary, binary_length, 1, fp);
      if (fread_sz == 0) {
          printf("Failed to read from the AOCX file (fread).\n");
          fclose(fp);
          free(const_cast<unsigned char*>(binary));
          cleanup_host_side_resources();
          exit(1);
      }
      fclose(fp);

      // Create a program using clCreateProgramWithBinary()
      program = clCreateProgramWithBinary(
                      context,
                      1,
                      devices,
                      &binary_length,
                      (const unsigned char **)&binary,
                      &status,
                      NULL); CHECK(status);


      //----------------------------------------------
      // Create the kernel
      //----------------------------------------------

      status = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
      if(status != CL_SUCCESS) {
          char log[10000] = {0};
          clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 10000, log, NULL);
          printf("%s\n", log);
          CHECK(status);
      }


      for(int j=0; j<NUM_KERNELS_TO_CREATE; j++) {
          printf("Creating kernel[%d]: %s\n", j,kernel_name[j]);
          kernel[j] = clCreateKernel(program, (const char*)kernel_name[j], &status);
          CHECK(status);
      }

      unsigned int mat_a_num_vectors_in_row_of_blocks = MAT_A_NUM_VECTORS_IN_ROW_OF_BLOCKS;
      unsigned char mat_a_num_blocks_in_col = MAT_A_NUM_BLOCKS_IN_COL;
      unsigned char mat_b_num_blocks_in_row = MAT_B_NUM_BLOCKS_IN_ROW;
      
      unsigned int mat_b_num_vectors_in_col_of_blocks = MAT_B_NUM_VECTORS_IN_COL_OF_BLOCKS;
      unsigned int mat_b_num_vectors_in_matrix = MAT_B_NUM_VECTORS_IN_MATRIX;
      int mat_c_num_coalesced_words = WC * HC / PE_COLS;
      
      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A1],
          0,
          sizeof(cl_mem),
          (void*)&d_matrix_mul_inputA[0]); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A1],
          1,
          sizeof(unsigned int),
          (void*)&mat_a_num_vectors_in_row_of_blocks); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A1],
          2,
          sizeof(unsigned char),
          (void*)&mat_a_num_blocks_in_col); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A1],
          3,
          sizeof(unsigned char),
          (void*)&mat_b_num_blocks_in_row); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A1],
          4,
          sizeof(unsigned char),
          (void *)&disableA); CHECK(status);


      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B1],
          0,
          sizeof(cl_mem),
          (void*)&d_matrix_mul_inputB[0]); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B1],
          1,
          sizeof(unsigned int),
          (void*)&mat_b_num_vectors_in_col_of_blocks); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B1],
          2,
          sizeof(unsigned int),
          (void*)&mat_b_num_vectors_in_matrix); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B1],
          3,
          sizeof(unsigned char),
          (void*)&mat_a_num_blocks_in_col); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B1],
          4,
          sizeof(unsigned char),
          (void *)&disableB); CHECK(status);


      status = clSetKernelArg(
          kernel[KID_DRAIN_MAT_C1],
          0,
          sizeof(cl_mem),
          (void*)&d_matrix_mul_outputC[0]); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_DRAIN_MAT_C1],
          1,
          sizeof(int),
          (void*)&mat_c_num_coalesced_words); CHECK(status);


      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A2],
          0,
          sizeof(cl_mem),
          (void*)&d_matrix_mul_inputA[0]); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A2],
          1,
          sizeof(unsigned int),
          (void*)&mat_a_num_vectors_in_row_of_blocks); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A2],
          2,
          sizeof(unsigned char),
          (void*)&mat_a_num_blocks_in_col); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A2],
          3,
          sizeof(unsigned char),
          (void*)&mat_b_num_blocks_in_row); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A2],
          4,
          sizeof(unsigned char),
          (void *)&disableA); CHECK(status);


      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B2],
          0,
          sizeof(cl_mem),
          (void*)&d_matrix_mul_inputB[0]); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B2],
          1,
          sizeof(unsigned int),
          (void*)&mat_b_num_vectors_in_col_of_blocks); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B2],
          2,
          sizeof(unsigned int),
          (void*)&mat_b_num_vectors_in_matrix); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B2],
          3,
          sizeof(unsigned char),
          (void*)&mat_a_num_blocks_in_col); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B2],
          4,
          sizeof(unsigned char),
          (void *)&disableB); CHECK(status);


      status = clSetKernelArg(
          kernel[KID_DRAIN_MAT_C2],
          0,
          sizeof(cl_mem),
          (void*)&d_matrix_mul_outputC[0]); CHECK(status);

      status = clSetKernelArg(
          kernel[KID_DRAIN_MAT_C2],
          1,
          sizeof(int),
          (void*)&mat_c_num_coalesced_words); CHECK(status);

      //----------------------------------------------
      // Configure the work-item structure (using only tasks atm)
      //----------------------------------------------

      // Define the number of threads that will be created
      // as well as the number of work groups

      // all kernels are always tasks
      globalWorkSize[0] = 1;
      localWorkSize[0]  = 1;

      
      //----------------------------------------------
      // Create command queues
      //---------------------------------------------

      // Create a command queue using clCreateCommandQueue(),
      // and associate it with the device you want to execute on
      int i;
      for(i=0; i<NUM_QUEUES_TO_CREATE; i++) {
                      fprintf(stdout,"cmdQueue i = %d, kernel name = %s\n", i, kernel_name[i]);
                      cmdQueue[i] = clCreateCommandQueue(
                              context,
                              devices[0],
                              CL_QUEUE_PROFILING_ENABLE,
                              &status); CHECK(status);
      }
      fprintf(stdout,"cmdQueue i = %d, a queue for reading the C buffer\n", i);
      cmdQueue[i] = clCreateCommandQueue(context,
                                          devices[0],
                                          CL_QUEUE_PROFILING_ENABLE,
                                          &status); CHECK(status);
      fpga_init_done=1;
  }

  int mat_init_done=0;
  void initmat()
  {
      if(mat_init_done==1)return;
      num_elem_A = HA*WA;
      num_elem_B = HB*WB;
      num_elem_C = HC*WC;

      if((matrix_mul_inputA1 = (float*)acl_aligned_malloc(num_elem_A*sizeof(float))) == NULL) {
          perror("Failed malloc of matrix_mul_inputA1");
          exit(1);
      }
      if((matrix_mul_inputA2 = (float*)acl_aligned_malloc(num_elem_A*sizeof(float))) == NULL) {
          perror("Failed malloc of matrix_mul_inputA2");
          exit(1);
      }
      if((matrix_mul_inputB1 = (float*)acl_aligned_malloc(num_elem_B*sizeof(float))) == NULL) {
          perror("Failed malloc of matrix_mul_inputB1");
          exit(1);
      }
      if((matrix_mul_inputB2 = (float*)acl_aligned_malloc(num_elem_B*sizeof(float))) == NULL) {
          perror("Failed malloc of matrix_mul_inputB2");
          exit(1);
      }
      if((matrix_mul_outputC1 = (float*)acl_aligned_malloc(num_elem_C*sizeof(float))) == NULL) {
          perror("Failed malloc of matrix_mul_outputC1");
          exit(1);
      }
      if((matrix_mul_outputC2 = (float*)acl_aligned_malloc(num_elem_C*sizeof(float))) == NULL) {
          perror("Failed malloc of matrix_mul_outputC2");
          exit(1);
      }
      if((matrix_mul_inputB_transposed = (float*)acl_aligned_malloc(num_elem_B*sizeof(float))) == NULL) {
          perror("Failed malloc of matrix_mul_inputB_transposed");
          exit(1);
      }

      printf("Allocated memory for host-side matrices!\n");
      mat_init_done=1;

  }
extern "C" 
{
void multiply(int dimm_np,int dimn_np,int dimk_np,float* At,float* Bt,float * Ct,int test)
{
      bool tprint=false;
      if(tprint)printf("dimm_np= %d\n",dimm_np);
      if(tprint)printf("dimn_np= %d\n",dimn_np);
      if(tprint)printf("dimk_np= %d\n",dimk_np);
      if(tprint)printf("A(1,1)=%f %f\n",At[0],At[dimm_np*dimk_np-1]);
      if(tprint)printf("B(1,1)=%f %f\n",Bt[0],Bt[dimk_np*dimn_np-1]);
      if(tprint)printf("C(1,1)=%f %f\n",Ct[0],Ct[dimm_np*dimn_np-1]);
      initmat();
      initfpga();
      disableA = 0;
      disableB = 0;
      double t0=omp_get_wtime();
      cl_int status;
      
      int HA2=dimm_np/HA;
      if(dimm_np%HA!=0)HA2++;
      int WA2=dimk_np/WA;
      if(dimk_np%WA!=0)WA2++;
      int HB2=dimk_np/HB;
      if(dimk_np%HB!=0)HB2++;
      int WB2=dimn_np/WB;
      if(dimn_np%WB!=0)WB2++;
      int HC2=dimm_np/HC;
      if(dimm_np%HC!=0)HC2++;
      int WC2=dimn_np/WC;
      if(dimn_np%WC!=0)WC2++;

      int dimm=HA2*HA;
      int dimn=WB2*WB;
      int dimk=WA2*WA;
      if(tprint)printf("dimm= %d\n",dimm);
      if(tprint)printf("dimn= %d\n",dimn);
      if(tprint)printf("dimk= %d\n",dimk);

      int elema=dimm_np*dimk_np;
      int elemb=dimk_np*dimn_np;
      int elemc=dimm_np*dimn_np;

      //decompose multiplication into kernel-compatible chunks
      printf("chunkcounts of A: %d x %d\n",HA2,WA2);
      printf("chunkcounts of B: %d x %d\n",HB2,WB2);
      printf("chunkcounts of C: %d x %d\n",HC2,WC2);
      
      printf("padded size of A: %d x %d\n",HA2*HA,WA*WA2);
      printf("padded size of B: %d x %d\n",HB2*HB,WB*WB2);
      printf("padded size of C: %d x %d\n",HC2*HC,WC*WC2);

      //list multiplication we have to do
      std::vector<std::vector<int> > mults0;
      std::vector<std::vector<int> > mults;
     
      int tadd;

      std::vector<std::vector<int> > Ablocks;
      std::vector<std::vector<int> > Bblocks;
      std::vector<std::vector<int> > Cblocks;
      std::vector<std::vector<int> > Cresult;
      for(int i=0;i<HA2;i++)
      {
        for(int j=0;j<WA2;j++)
        {
          std::vector<int> block;
          block.push_back(i);
          block.push_back(j);
          block.push_back(i*HA);
          block.push_back(min((i+1)*HA-1,dimm_np-1));
          block.push_back(j*WA);
          block.push_back(min((j+1)*WA-1,dimk_np-1));
          block.push_back(0); //0 means not on device
          block.push_back(0); //0 means transfer has been planned
          Ablocks.push_back(block);
          int ind=Ablocks.size()-1;
          printf("block %d (%d,%d) of A is the chunk A(%d:%d,%d:%d)\n",Ablocks.size(),Ablocks[ind][0],Ablocks[ind][1],Ablocks[ind][2],Ablocks[ind][3],Ablocks[ind][4],Ablocks[ind][5]);
        }
      }
      printf("Blocks for A=%d\n",Ablocks.size());
      
      for(int i=0;i<HB2;i++)
      {
        for(int j=0;j<WB2;j++)
        {
          std::vector<int> block;
          block.push_back(i);
          block.push_back(j);
          block.push_back(i*HB);
          block.push_back(min((i+1)*HB-1,dimk_np-1));
          block.push_back(j*WB);
          block.push_back(min((j+1)*WB-1,dimn_np-1));
          block.push_back(0); //0 means not on device
          block.push_back(0); //0 means transfer has been planned
          Bblocks.push_back(block);
          int ind=Bblocks.size()-1;
          printf("block %d (%d,%d) of B is the chunk B(%d:%d,%d:%d)\n",Bblocks.size(),Bblocks[ind][0],Bblocks[ind][1],Bblocks[ind][2],Bblocks[ind][3],Bblocks[ind][4],Bblocks[ind][5]);
        }
      }
      printf("Blocks for B=%d\n",Bblocks.size());


      for(int crow=0;crow<HC2;crow++){
        std::vector<int> col(WC2);
        Cresult.push_back(col);
      }

      for(int crow=0;crow<HC2;crow++){
        for(int ccol=0;ccol<WC2;ccol++){
          printf("for block %d,%d of C we need:\n",crow,ccol);

          for(int acol=0;acol<WA2;acol++){
            std::vector<int> block;

            block.push_back(crow);
            block.push_back(ccol);
            block.push_back(crow*HC);
            block.push_back(min((crow+1)*HC-1,dimm_np-1));
            block.push_back(ccol*WC);
            block.push_back(min((ccol+1)*WC-1,dimn_np-1));
            block.push_back(0); //0 means not on device
            Cblocks.push_back(block);

            int brow=acol;
            int aid=get_block_id(Ablocks,crow,acol);
            int bid=get_block_id(Bblocks,brow,ccol);
            std::vector<int> mtmp;
            mtmp.push_back(aid);
            mtmp.push_back(bid);
            mtmp.push_back(Cblocks.size()-1);
            mtmp.push_back(0); //not planned yet
            
            printf("A %d,%d (%d) times B %d,%d (%d) --> C %d \n",crow,acol,aid,brow,ccol,bid,Cblocks.size()-1);

            mults0.push_back(mtmp);
          }
        }
      }
      printf("Blocks for C=%d\n",Cblocks.size());
      int nmults=mults0.size();
      printf("%i mults in total.\n",nmults);

    //replan multiplication so that at most one input matrix is copied per stop

    std::vector<int> order(nmults);
    std::vector<std::vector<int> > transferlistA(nmults/2+2);
    std::vector<std::vector<int> > transferlistB(nmults/2+2);
    std::vector<std::vector<int> > transferlistC(nmults/2+2);

/*
for block 0,0 of C we need:
A 0,0 (0) times B 0,0 (0) --> C 0 
A 0,1 (1) times B 1,0 (1) --> C 1 
for block 1,0 of C we need:
A 1,0 (2) times B 0,0 (0) --> C 2 
A 1,1 (3) times B 1,0 (1) --> C 3 
for block 2,0 of C we need:
A 2,0 (4) times B 0,0 (0) --> C 4 
A 2,1 (5) times B 1,0 (1) --> C 5 
for block 3,0 of C we need:
A 3,0 (6) times B 0,0 (0) --> C 6 
A 3,1 (7) times B 1,0 (1) --> C 7 
for block 4,0 of C we need:
A 4,0 (8) times B 0,0 (0) --> C 8 
A 4,1 (9) times B 1,0 (1) --> C 9 
for block 5,0 of C we need:
A 5,0 (10) times B 0,0 (0) --> C 10 
A 5,1 (11) times B 1,0 (1) --> C 11 
*/
/*
padded size of A: 4096 x 4096
padded size of B: 4096 x 4096
padded size of C: 4096 x 4096
block 1 (0,0) of A is the chunk A(0:2047,0:4095)
block 2 (1,0) of A is the chunk A(2048:4095,0:4095)
Blocks for A=2
block 1 (0,0) of B is the chunk B(0:4095,0:2047)
block 2 (0,1) of B is the chunk B(0:4095,2048:4095)
Blocks for B=2
for block 0,0 of C we need:
A 0,0 (0) times B 0,0 (0) --> C 0 
for block 0,1 of C we need:
A 0,0 (0) times B 0,1 (1) --> C 1 
for block 1,0 of C we need:
A 1,0 (1) times B 0,0 (0) --> C 2 
for block 1,1 of C we need:
A 1,0 (1) times B 0,1 (1) --> C 3 
Blocks for C=4
4 mults in total.
*/
    if(0==1){
      order[0]=0;
      order[1]=2;

      order[2]=1;
      order[3]=3;
/*
      order[0]=0;
      order[1]=2;

      order[2]=4;
      order[3]=6;
      
      order[4]=8;
      order[5]=10;
      
      order[6]=1;
      order[7]=3;
      
      order[8]=5;
      order[9]=7;
      
      order[10]=9;
      order[11]=11;*/
    }else if(1==1){
      order[0]=0;
      order[1]=2;

      order[2]=4;
      order[3]=6;
      
      order[4]=1;
      order[5]=3;
      
      order[6]=5;
      order[7]=7;
    }else{
      for(int i=0;i<nmults;i++)
      {
        order[i]=i;
      }
    }

    for(int i=0;i<nmults;i++)
    {
      mults.push_back(mults0[order[i]]);
    }

    

    //create chunks
    unsigned int i;
    std::streampos filesize;
    FILE *f_out = stdout;


    //perform chunked mults to test mapping
    if(test>=1){
      Ct_ref = (float*)malloc(elemc*sizeof(float));
      float alpha=1.0;
      float beta=0.0;
      //compute reference value
      double cpu0=omp_get_wtime();
      cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dimm_np,dimn_np,dimk_np,alpha,At,dimm_np,Bt,dimk_np,beta,Ct_ref,dimm_np);
      double cpu1=omp_get_wtime();
      printf("duration on cpu %f seconds, %f GFLOPS\n",cpu1-cpu0,(double)2.0*dimm_np*dimn_np*dimk_np*pow(10,-9)/(cpu1-cpu0));
    }
    if(test>=2){
      Ct_chunked = (float*)malloc(elemc*sizeof(float));
      memset(Ct_chunked, 0, elemc*sizeof(float));
      for(int i=0;i<nmults;i++){
        get_chunk_from_matrix(dimm_np,dimk_np,HA,WA,Ablocks[mults[i][0]],At,matrix_mul_inputA1);
        get_chunk_from_matrix(dimk_np,dimn_np,HB,WB,Bblocks[mults[i][1]],Bt,matrix_mul_inputB1);
        float alpha=1.0;
        float beta=0.0;
        cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,HA,WB,WA,alpha,matrix_mul_inputA1,WA,matrix_mul_inputB1,WB,beta,matrix_mul_outputC1,WC);
        set_matrix_from_chunk_add(dimm_np,dimn_np,HC,WC,Cblocks[mults[i][2]],Ct_chunked,matrix_mul_outputC1,true);
      }
      double diff=0;
      for(int j=0;j<elemc;j++)
      {
        diff+=abs(Ct_ref[j]-Ct_chunked[j]);
      }
      printf("diff_chunked=%f\n",diff/(elemc));
    }
/*    for(int i=0;i<elemc;i++)
    {
      Ct[i]=Ct_chunked[i];
    }
    return;
    */


    //----------------------------------------------
    // Enqueue the kernel for execution
    //----------------------------------------------

    double t1=omp_get_wtime();
    printf("init took %f seconds \n",t1-t0);
    if(tprint)printf("\n===== Host-CPU enqeuing the OpenCL kernels to the FPGA device ======\n\n");
    // blocking writes
    if(Ablocks[mults[0][0]][6]==0){
      double ct1=omp_get_wtime();
      if(tprint)printf("writing A %d %d\n",0,mults[0][0]);
      block_wise_reformat_chunk(dimm_np,dimk_np,At,HA,WA,Ablocks[mults[0][0]],matrix_mul_inputA1, MAT_A_BLOCK_HEIGHT, MAT_A_BLOCK_WIDTH);
      
      for(int i=d_matrix_mul_inputA.size();i<=mults[0][0];i++){
        if(tprint)printf("creating device buffer B %d\n",i);
        cl_mem tmp;
        tmp = clCreateBuffer(
              context,
              CL_MEM_READ_ONLY,
              num_elem_A*sizeof(cl_float),
              NULL,
              &status); CHECK(status);
        d_matrix_mul_inputA.push_back(tmp);
      }
      status = clEnqueueWriteBuffer(
              cmdQueue[KID_HOST_COMM],
              d_matrix_mul_inputA[mults[0][0]],
              CL_FALSE,
              0,
              num_elem_A*sizeof(cl_float),
              matrix_mul_inputA1,
              0,
              NULL,
              NULL); CHECK(status);
      double ct2=omp_get_wtime();
      printf("transfer rate %f\n",num_elem_A*sizeof(cl_float)*pow(10,-6)/(ct2-ct1));
       Ablocks[mults[0][0]][6]=1;
    }

    if(Bblocks[mults[0][1]][6]==0){
      if(tprint)printf("writing B %d %d\n",0,mults[0][1]);
      get_chunk_from_matrix(dimk_np,dimn_np,HB,WB,Bblocks[mults[0][1]],Bt,matrix_mul_inputB1);
      transpose_matrix(matrix_mul_inputB1, matrix_mul_inputB_transposed, HB, WB);
      block_wise_reformat_matrix(matrix_mul_inputB_transposed, matrix_mul_inputB1, WB, HB, MAT_B_BLOCK_WIDTH, MAT_B_BLOCK_HEIGHT);
      
      for(int i=d_matrix_mul_inputB.size();i<=mults[0][1];i++){
        if(tprint)if(tprint)printf("creating device buffer B %d\n",i);
        cl_mem tmp;
        tmp = clCreateBuffer(
              context,
              CL_MEM_READ_ONLY,
              num_elem_B*sizeof(cl_float),
              NULL,
              &status); CHECK(status);
        d_matrix_mul_inputB.push_back(tmp);
      }
      status = clEnqueueWriteBuffer(
              cmdQueue[KID_HOST_COMM],
              d_matrix_mul_inputB[mults[0][1]],
              CL_FALSE,
              0,
              num_elem_B*sizeof(cl_float),
              matrix_mul_inputB1,
              0,
              NULL,
              NULL); CHECK(status);
        Bblocks[mults[0][1]][6]=1;
    }

    if(Ablocks[mults[1][0]][6]==0){
      if(tprint)printf("writing A %d %d\n",0,mults[1][0]);
      block_wise_reformat_chunk(dimm_np,dimk_np,At,HA,WA,Ablocks[mults[1][0]],matrix_mul_inputA2, MAT_A_BLOCK_HEIGHT, MAT_A_BLOCK_WIDTH);
      for(int i=d_matrix_mul_inputA.size();i<=mults[1][0];i++){
        if(tprint)printf("creating device buffer A %d\n",i);
        cl_mem tmp;
        tmp = clCreateBuffer(
              context,
              CL_MEM_READ_ONLY,
              num_elem_A*sizeof(cl_float),
              NULL,
              &status); CHECK(status);
        d_matrix_mul_inputA.push_back(tmp);
      }
      status = clEnqueueWriteBuffer(
              cmdQueue[KID_HOST_COMM],
              d_matrix_mul_inputA[mults[1][0]],
              CL_FALSE,
              0,
              num_elem_A*sizeof(cl_float),
              matrix_mul_inputA2,
              0,
              NULL,
              NULL); CHECK(status);
      Ablocks[mults[1][0]][6]=1;
    }

    if(Bblocks[mults[1][1]][6]==0){
      if(tprint)printf("writing B %d %d\n",0,mults[1][1]);
      get_chunk_from_matrix(dimk_np,dimn_np,HB,WB,Bblocks[mults[1][1]],Bt,matrix_mul_inputB2);
      transpose_matrix(matrix_mul_inputB2, matrix_mul_inputB_transposed, HB, WB);
      block_wise_reformat_matrix(matrix_mul_inputB_transposed, matrix_mul_inputB2, WB, HB, MAT_B_BLOCK_WIDTH, MAT_B_BLOCK_HEIGHT);
      for(int i=d_matrix_mul_inputB.size();i<=mults[1][1];i++){
        if(tprint)printf("creating device buffer B %d\n",i);
        cl_mem tmp;
        tmp = clCreateBuffer(
              context,
              CL_MEM_READ_ONLY,
              num_elem_B*sizeof(cl_float),
              NULL,
              &status); CHECK(status);
        d_matrix_mul_inputB.push_back(tmp);
      }
      status = clEnqueueWriteBuffer(
              cmdQueue[KID_HOST_COMM],
              d_matrix_mul_inputB[mults[1][1]],
              CL_FALSE,
              0,
              num_elem_B*sizeof(cl_float),
              matrix_mul_inputB2,
              0,
              NULL,
              NULL); CHECK(status);
        Bblocks[mults[1][1]][6]=1;
    }
        status = clFlush(cmdQueue[KID_HOST_COMM]);
        status = clFinish(cmdQueue[KID_HOST_COMM]); CHECK(status);

    double t2=omp_get_wtime();
    printf("first pciewrite took %f seconds \n",t2-t1);
    if(tprint)printf("\n");
    if(tprint)printf("loop start\n");

    std::vector<int>cdone;

    for(int iter=0;iter<mults.size();iter+=2)
    {
      for(int i=d_matrix_mul_outputC.size();i<=mults[iter][2];i++){
        if(tprint)printf("creating device buffer C %d\n",i);
        cl_mem tmp;
        tmp = clCreateBuffer(
              context,
              CL_MEM_WRITE_ONLY,
              num_elem_C*sizeof(cl_float),
              NULL,
              &status); CHECK(status);
        d_matrix_mul_outputC.push_back(tmp);
      }
      int transferredA=0;
      int transferredB=0;
      double t1=omp_get_wtime();
      if(tprint)printf("setting A %d %d as argument\n",iter,mults[iter][0]);
      status = clSetKernelArg(
          kernel[KID_FEED_MAT_A1],
          0,
          sizeof(cl_mem),
          (void*)&(d_matrix_mul_inputA[mults[iter][0]])); CHECK(status);
      
      if(tprint)printf("setting B %d %d as argument\n",iter,mults[iter][1]);
      status = clSetKernelArg(
          kernel[KID_FEED_MAT_B1],
          0,
          sizeof(cl_mem),
          (void*)&(d_matrix_mul_inputB[mults[iter][1]])); CHECK(status);

      if(tprint)printf("setting C %d %d as argument\n",iter,mults[iter][2]);
      status = clSetKernelArg(
          kernel[KID_DRAIN_MAT_C1],
          0,
          sizeof(cl_mem),
          (void*)&(d_matrix_mul_outputC[mults[iter][2]])); CHECK(status);
     
     if(iter+1<nmults){ 
        for(int i=d_matrix_mul_outputC.size();i<=mults[iter+1][2];i++){
          printf("creating device buffer C %d\n",i);
          cl_mem tmp;
          tmp = clCreateBuffer(
                context,
                CL_MEM_WRITE_ONLY,
                num_elem_C*sizeof(cl_float),
                NULL,
                &status); CHECK(status);
          d_matrix_mul_outputC.push_back(tmp);
        }
        if(tprint)printf("setting A %d %d as argument\n",iter+1,mults[iter+1][0]);
        status = clSetKernelArg(
            kernel[KID_FEED_MAT_A2],
            0,
            sizeof(cl_mem),
            (void*)&(d_matrix_mul_inputA[mults[iter+1][0]])); CHECK(status);
        
        if(tprint)printf("setting B %d %d as argument\n",iter+1,mults[iter+1][1]);
        status = clSetKernelArg(
            kernel[KID_FEED_MAT_B2],
            0,
            sizeof(cl_mem),
            (void*)&(d_matrix_mul_inputB[mults[iter+1][1]])); CHECK(status);

        if(tprint)printf("setting C %d %d as argument\n",iter+1,mults[iter+1][2]);
        status = clSetKernelArg(
            kernel[KID_DRAIN_MAT_C2],
            0,
            sizeof(cl_mem),
            (void*)&(d_matrix_mul_outputC[mults[iter+1][2]])); CHECK(status);
      }
        
      for(i=0; i<NUM_QUEUES_TO_CREATE; i++) {
        clReleaseCommandQueue(cmdQueue[i]);
      }
      clReleaseCommandQueue(cmdQueue[KID_HOST_COMM]);

      for(i=0; i<NUM_QUEUES_TO_CREATE; i++) {
                      cmdQueue[i] = clCreateCommandQueue(
                              context,
                              devices[0],
                              CL_QUEUE_PROFILING_ENABLE,
                              &status); CHECK(status);
      }

      cmdQueue[KID_HOST_COMM] = clCreateCommandQueue(context,
                                          devices[0],
                                          CL_QUEUE_PROFILING_ENABLE,
                                          &status); CHECK(status);

      double t2=omp_get_wtime();
      if(iter%100==0 && iter>0){
        if(tprint)printf("resetting queues to prevent slowdown iter=%d\n",iter);

        for(i=0; i<NUM_QUEUES_TO_CREATE; i++) {
          clReleaseCommandQueue(cmdQueue[i]);
        }
        clReleaseCommandQueue(cmdQueue[KID_HOST_COMM]);

        for(i=0; i<NUM_QUEUES_TO_CREATE; i++) {
                        cmdQueue[i] = clCreateCommandQueue(
                                context,
                                devices[0],
                                CL_QUEUE_PROFILING_ENABLE,
                                &status); CHECK(status);
        }

        cmdQueue[KID_HOST_COMM] = clCreateCommandQueue(context,
                                            devices[0],
                                            CL_QUEUE_PROFILING_ENABLE,
                                            &status); CHECK(status);
      }
      double t3=omp_get_wtime();

      if(tprint)printf("running multiplication %d\n",iter);
      for(i=0; i<NUM_KERNELS_TO_CREATE; i++) 
      {

          // Alternatively, can use clEnqueueTaskKernel
          //printf("clEnqueueNDRangeKernel[%d]: %s!\n", i,kernel_name[i]);
          status = clEnqueueNDRangeKernel(
                          cmdQueue[i],
                          kernel[i],
                          1,
                          NULL,
                          globalWorkSize,
                          localWorkSize,
                          0,
                          NULL,
                          &kernel_exec_event[i]
                          );
          CHECK(status);
      }

      double t4=omp_get_wtime();
      if(iter+2<mults.size()){
        //copy new A and B to device
        int ind=iter+2;
        if(Ablocks[mults[ind][0]][6]==0){
          if(tprint)printf("writing A %d %d\n",ind,mults[ind][0]);
          block_wise_reformat_chunk(dimm_np,dimk_np,At,HA,WA,Ablocks[mults[ind][0]],matrix_mul_inputA1, MAT_A_BLOCK_HEIGHT, MAT_A_BLOCK_WIDTH);
          double ct1=omp_get_wtime();

          for(int i=d_matrix_mul_inputA.size();i<=mults[ind][0];i++){
            if(tprint)printf("creating device buffer A %d\n",i);
            cl_mem tmp;
            tmp = clCreateBuffer(
                  context,
                  CL_MEM_READ_ONLY,
                  num_elem_A*sizeof(cl_float),
                  NULL,
                  &status); CHECK(status);
            d_matrix_mul_inputA.push_back(tmp);
          }
          status = clEnqueueWriteBuffer(
                  cmdQueue[KID_HOST_COMM],
                  d_matrix_mul_inputA[mults[ind][0]],
                  CL_TRUE,
                  0,
                  num_elem_A*sizeof(cl_float),
                  matrix_mul_inputA1,
                  0,
                  NULL,
                  NULL); CHECK(status);
          double ct2=omp_get_wtime();
          printf("transfer rate %f\n",num_elem_A*sizeof(cl_float)*pow(10,-6)/(ct2-ct1));
          Ablocks[mults[ind][0]][6]=1;
          transferredA++;
        }

        if(Bblocks[mults[ind][1]][6]==0){
          if(tprint)printf("writing B %d %d\n",ind,mults[mults[ind][1]][1]);
          get_chunk_from_matrix(dimk_np,dimn_np,HB,WB,Bblocks[mults[ind][1]],Bt,matrix_mul_inputB1);
          transpose_matrix(matrix_mul_inputB1, matrix_mul_inputB_transposed, HB, WB);
          block_wise_reformat_matrix(matrix_mul_inputB_transposed, matrix_mul_inputB1, WB, HB, MAT_B_BLOCK_WIDTH, MAT_B_BLOCK_HEIGHT);
          double ct1=omp_get_wtime();
          
          for(int i=d_matrix_mul_inputB.size();i<=mults[ind][1];i++){
            if(tprint)printf("creating device buffer B %d\n",i);
            cl_mem tmp;
            tmp = clCreateBuffer(
                  context,
                  CL_MEM_READ_ONLY,
                  num_elem_B*sizeof(cl_float),
                  NULL,
                  &status); CHECK(status);
            d_matrix_mul_inputB.push_back(tmp);
          }
          status = clEnqueueWriteBuffer(
                  cmdQueue[KID_HOST_COMM],
                  d_matrix_mul_inputB[mults[ind][1]],
                  CL_TRUE,
                  0,
                  num_elem_B*sizeof(cl_float),
                  matrix_mul_inputB1,
                  0,
                  NULL,
                  NULL); CHECK(status);
          double ct2=omp_get_wtime();
          printf("transfer rate %f\n",num_elem_B*sizeof(cl_float)*pow(10,-6)/(ct2-ct1));
          Bblocks[mults[ind][1]][6]=1;
          transferredB++;
        }
      }
      if(iter+3<mults.size()){
        int ind=iter+3;
        if(Ablocks[mults[ind][0]][6]==0){
          if(tprint)printf("writing A %d %d\n",ind,mults[ind][0]);
          block_wise_reformat_chunk(dimm_np,dimk_np,At,HA,WA,Ablocks[mults[ind][0]],matrix_mul_inputA2, MAT_A_BLOCK_HEIGHT, MAT_A_BLOCK_WIDTH);
          double ct1=omp_get_wtime();
          
          for(int i=d_matrix_mul_inputA.size();i<=mults[ind][0];i++){
            if(tprint)printf("creating device buffer A %d\n",i);
            cl_mem tmp;
            tmp = clCreateBuffer(
                  context,
                  CL_MEM_READ_ONLY,
                  num_elem_A*sizeof(cl_float),
                  NULL,
                  &status); CHECK(status);
            d_matrix_mul_inputA.push_back(tmp);
          }
          status = clEnqueueWriteBuffer(
                  cmdQueue[KID_HOST_COMM],
                  d_matrix_mul_inputA[mults[ind][0]],
                  CL_TRUE,
                  0,
                  num_elem_A*sizeof(cl_float),
                  matrix_mul_inputA2,
                  0,
                  NULL,
                  NULL); CHECK(status);
          double ct2=omp_get_wtime();
          printf("transfer rate %f\n",num_elem_A*sizeof(cl_float)*pow(10,-6)/(ct2-ct1));
          Ablocks[mults[ind][0]][6]=1;
          transferredA++;
        }

        if(Bblocks[mults[ind][1]][6]==0){
          if(tprint)printf("writing B %d %d\n",ind,mults[ind][1]);
          get_chunk_from_matrix(dimk_np,dimn_np,HB,WB,Bblocks[mults[ind][1]],Bt,matrix_mul_inputB2);
          transpose_matrix(matrix_mul_inputB2, matrix_mul_inputB_transposed, HB, WB);
          block_wise_reformat_matrix(matrix_mul_inputB_transposed, matrix_mul_inputB2, WB, HB, MAT_B_BLOCK_WIDTH, MAT_B_BLOCK_HEIGHT);
          double ct1=omp_get_wtime();
          
          for(int i=d_matrix_mul_inputB.size();i<=mults[ind][1];i++){
            if(tprint)printf("creating device buffer B %d\n",i);
            cl_mem tmp;
            tmp = clCreateBuffer(
                  context,
                  CL_MEM_READ_ONLY,
                  num_elem_B*sizeof(cl_float),
                  NULL,
                  &status); CHECK(status);
            d_matrix_mul_inputB.push_back(tmp);
          }
          status = clEnqueueWriteBuffer(
                  cmdQueue[KID_HOST_COMM],
                  d_matrix_mul_inputB[mults[ind][1]],
                  CL_TRUE,
                  0,
                  num_elem_B*sizeof(cl_float),
                  matrix_mul_inputB2,
                  0,
                  NULL,
                  NULL); CHECK(status);
          double ct2=omp_get_wtime();
          printf("transfer rate %f\n",num_elem_B*sizeof(cl_float)*pow(10,-6)/(ct2-ct1));
          Bblocks[mults[ind][1]][6]=1;
          transferredB++;
        }
      }
      
      status = clFlush(cmdQueue[KID_HOST_COMM]);
      status = clFinish(cmdQueue[KID_HOST_COMM]); CHECK(status);
      double t5=omp_get_wtime();
      double t6=omp_get_wtime();
      //copy new C from device
      if(transferredB==0){
        int kmax=6;
        if(transferredA>0)kmax=1;
        for(int k=0;k<kmax;k++)
        {
          if(cdone.size()>0){
            int ind1=cdone[0];
            int ind2=cdone[1];
            if(tprint)printf("reading C %d %d\n",ind1,mults[ind1][2]);
            clEnqueueReadBuffer(
                        cmdQueue[KID_HOST_COMM], // using a special queue for reading buffer C
                        d_matrix_mul_outputC[mults[ind1][2]],
                        CL_TRUE,
                        0,
                        num_elem_C*sizeof(cl_float),
                        matrix_mul_outputC1,
                        0,
                        NULL,
                        NULL); CHECK(status);
            
            if(tprint)printf("reading C %d %d\n",ind2,mults[ind2][2]);
            clEnqueueReadBuffer(
                        cmdQueue[KID_HOST_COMM], // using a special queue for reading buffer C
                        d_matrix_mul_outputC[mults[ind2][2]],
                        CL_FALSE,
                        0,
                        num_elem_C*sizeof(cl_float),
                        matrix_mul_outputC2,
                        0,
                        NULL,
                        NULL); CHECK(status);

            t6=omp_get_wtime();
            
            tadd=Cresult[Cblocks[mults[ind1][2]][0]][Cblocks[mults[ind1][2]][1]];
            set_matrix_from_chunk_reorder_add(dimm_np,dimn_np,Cblocks[mults[ind1][2]],Ct,matrix_mul_outputC1,HC,WC,MAT_C_BLOCK_HEIGHT,MAT_C_BLOCK_WIDTH,PE_COLS,tadd);
            Cresult[Cblocks[mults[ind1][2]][0]][Cblocks[mults[ind1][2]][1]]=1;

            status = clFlush(cmdQueue[KID_HOST_COMM]);
            status = clFinish(cmdQueue[KID_HOST_COMM]); CHECK(status);

            tadd=Cresult[Cblocks[mults[ind2][2]][0]][Cblocks[mults[ind2][2]][1]];
            set_matrix_from_chunk_reorder_add(dimm_np,dimn_np,Cblocks[mults[ind2][2]],Ct,matrix_mul_outputC2,HC,WC,MAT_C_BLOCK_HEIGHT,MAT_C_BLOCK_WIDTH,PE_COLS,tadd);
            Cresult[Cblocks[mults[ind2][2]][0]][Cblocks[mults[ind2][2]][1]]=1;
            
            cdone.erase(cdone.begin()+1);
            cdone.erase(cdone.begin()+0);
          }
        }
      }
      cdone.push_back(iter);
      cdone.push_back(iter+1);

      double t7=omp_get_wtime();

      for(i=0; i < NUM_KERNELS_TO_CREATE ; i++) {
          status = clFlush(cmdQueue[i]);
          CHECK(status);
      }
      for(i=0; i < NUM_QUEUES_TO_FINISH; i++) {
           status = clFinish(cmdQueue[i]); CHECK(status);
      }
      status = clFlush(cmdQueue[KID_HOST_COMM]);
      status = clFinish(cmdQueue[KID_HOST_COMM]); CHECK(status);
      double t8=omp_get_wtime();
      //printf("T %d %f %f %f %f %f %f %f %f\n",iter,t2-t1,t3-t2,t4-t3,t5-t4,t6-t5,t7-t6,t8-t7);
      printf("T %d %f %f %f %f %f %f %f %f\n",iter,t2-t1,t3-t1,t4-t1,t5-t1,t6-t1,t7-t1,t8-t1);
    
    }
    for(i=0; i < NUM_KERNELS_TO_CREATE ; i++) {
        status = clFlush(cmdQueue[i]);
        CHECK(status);
    }
    for(i=0; i < NUM_QUEUES_TO_FINISH; i++) {
         status = clFinish(cmdQueue[i]); CHECK(status);
    }
    double t3=omp_get_wtime();

    if(tprint)printf("loop took %f seconds\n",t3-t2);


    double k_start_time[NUM_QUEUES_TO_FINISH];
    double k_end_time[NUM_QUEUES_TO_FINISH];
    double k_exec_time[NUM_QUEUES_TO_FINISH];

    for (i=0; i<NUM_QUEUES_TO_FINISH; i++) {
        k_exec_time[i] = compute_kernel_execution_time(kernel_exec_event[i], k_start_time[i], k_end_time[i]);
    }
    if(tprint)printf("\n\n");

    if(tprint)printf("\n===== Host-CPU transferring result matrix C from the FPGA device global memory (DDR4) via PCIe ======\n\n");

    // Read the results back from the device, blocking read
    if(tprint)printf("remaining to be transferred %d\n",cdone.size());
    while(cdone.size()>0)
    {
      int ind1=cdone[0];
      int ind2=cdone[1];
      if(tprint)printf("reading C %d %d\n",ind1,mults[ind1][2]);
      clEnqueueReadBuffer(
                  cmdQueue[KID_HOST_COMM], // using a special queue for reading buffer C
                  d_matrix_mul_outputC[mults[ind1][2]],
                  CL_TRUE,
                  0,
                  num_elem_C*sizeof(cl_float),
                  matrix_mul_outputC1,
                  0,
                  NULL,
                  NULL); CHECK(status);
      
      if(tprint)printf("reading C %d %d\n",ind2,mults[ind2][2]);
      clEnqueueReadBuffer(
                  cmdQueue[KID_HOST_COMM], // using a special queue for reading buffer C
                  d_matrix_mul_outputC[mults[ind2][2]],
                  CL_FALSE,
                  0,
                  num_elem_C*sizeof(cl_float),
                  matrix_mul_outputC2,
                  0,
                  NULL,
                  NULL); CHECK(status);

      tadd=Cresult[Cblocks[mults[ind1][2]][0]][Cblocks[mults[ind1][2]][1]];
      set_matrix_from_chunk_reorder_add(dimm_np,dimn_np,Cblocks[mults[ind1][2]],Ct,matrix_mul_outputC1,HC,WC,MAT_C_BLOCK_HEIGHT,MAT_C_BLOCK_WIDTH,PE_COLS,tadd);
      Cresult[Cblocks[mults[ind1][2]][0]][Cblocks[mults[ind1][2]][1]]=1;
      
      status = clFlush(cmdQueue[KID_HOST_COMM]);
      status = clFinish(cmdQueue[KID_HOST_COMM]); CHECK(status);
        
      tadd=Cresult[Cblocks[mults[ind2][2]][0]][Cblocks[mults[ind2][2]][1]];
      set_matrix_from_chunk_reorder_add(dimm_np,dimn_np,Cblocks[mults[ind2][2]],Ct,matrix_mul_outputC2,HC,WC,MAT_C_BLOCK_HEIGHT,MAT_C_BLOCK_WIDTH,PE_COLS,tadd);
      Cresult[Cblocks[mults[ind2][2]][0]][Cblocks[mults[ind2][2]][1]]=1;
      
      cdone.erase(cdone.begin()+1);
      cdone.erase(cdone.begin()+0);
    }

    double t4=omp_get_wtime();
    printf("final pcieread took %f seconds \n",t4-t3);
    printf("total took %f seconds \n",t4-t1);
    printf("total GFLOPS without init, windup and winddown=%f %f \n",(double)1.0e-9*((double)2.0 * WA * HC * WC*mults.size())/(t3-t2),t3-t2);


    if(test>=1){ 
      double diff2=0;
      for(int j=0;j<elemc;j++)
      {
        diff2+=abs(Ct_ref[j]-Ct[j]);
      }
      printf("diff_sgemmcpu_sgemmfpga=%f\n",diff2/(elemc));
      free(Ct_ref);
      
      double diff3=0;
      if(test>=2){ 
        for(int j=0;j<elemc;j++)
        {
          diff3+=abs(Ct_chunked[j]-Ct[j]);
        }
        printf("diff_sgemmcpuchunked_sgemmfpga=%f\n",diff3/(elemc));
        free(Ct_chunked);
      }
    }
    fflush(stdout); 
    //cleanup();
}
}

void cleanup_host_side_resources() {
    
//    if (matrix_mul_inputB_transposed) acl_aligned_free(matrix_mul_inputB_transposed);
//    if (matrix_mul_inputB_block_wise) acl_aligned_free(matrix_mul_inputB_block_wise);
//    matrix_mul_inputA_block_wise = NULL;
    matrix_mul_inputB_transposed = NULL;
//    matrix_mul_inputB_block_wise = NULL;
}

// Free the resources allocated during initialization
void cleanup() {
    //----------------------------------------------
    // Release the OpenCL resources
    //----------------------------------------------
    int i;
    // Free resources
    for(i=0; i<NUM_KERNELS_TO_CREATE; i++) {
        clReleaseKernel(kernel[i]);
    }

    for(i=0; i<NUM_QUEUES_TO_CREATE; i++) {
        clReleaseCommandQueue(cmdQueue[i]);
    }

    for(i=0; i<NUM_QUEUES_TO_FINISH; i++) {
        clReleaseEvent(kernel_exec_event[i]);
    }

    for(i=0;i<d_matrix_mul_inputA.size();i++){
      clReleaseMemObject(d_matrix_mul_inputA[i]);
    }
    for(i=0;i<d_matrix_mul_inputB.size();i++){
      clReleaseMemObject(d_matrix_mul_inputB[i]);
    }
    for(i=0;i<d_matrix_mul_outputC.size();i++){
      clReleaseMemObject(d_matrix_mul_outputC[i]);
    }
    /*clReleaseMemObject(d_matrix_mul_inputA1);
    clReleaseMemObject(d_matrix_mul_inputB1);
    clReleaseMemObject(d_matrix_mul_outputC1);
    clReleaseMemObject(d_matrix_mul_inputA2);
    clReleaseMemObject(d_matrix_mul_inputB2);
    clReleaseMemObject(d_matrix_mul_outputC2);*/


    /*acl_aligned_free(matrix_mul_inputA);
    acl_aligned_free(matrix_mul_inputB);
    acl_aligned_free(matrix_mul_outputC);*/

    clReleaseProgram(program);
    clReleaseContext(context);

    acl_aligned_free(devices);

    cleanup_host_side_resources();
}


