#include<time.h>
#include<omp.h>
#include<string.h>
#include<stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <cuda_fp16.h>
#include <omp.h>
#include <stdint.h>
#include <immintrin.h>

#include "WHATTYPE.h"
//#define __type1 __half 
//#define __type1 float
//#define __type1 double 

//#define __type2 __half 
//#define __type2 float
//#define __type2 double 




template<typename AB_type, typename C_type>
struct CublasContext
{
    cublasHandle_t handle;
    AB_type* devPtrX;
    AB_type* devPtrY1;
    AB_type* devPtrY2;
    AB_type* devPtrY3;
    AB_type* devPtrones;
    AB_type* x;
};

bool inited=false;
void init_device()
{
    if(inited)return;
    cudaError_t cudaStat;
    //Query Device Properties
    int device;
    cudaStat = cudaGetDevice(&device);
    if(cudaStat == cudaSuccess)
    {
        struct cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device);
        printf("%s\n", prop.name);
    }
    else
    {
        printf("Failed to query device\n");
        exit(1);
    }
    inited=true;
}

template<typename AB_type, typename C_type>
int create_context(CublasContext<AB_type, C_type>& cont, int n)
{   
    cudaError_t cudaStat;
    cublasStatus_t stat;
    cont.handle = 0;
    cont.devPtrX = nullptr;
    cont.devPtrY1 = nullptr;
    cont.devPtrY2 = nullptr;
    cont.devPtrY3 = nullptr;
    cont.devPtrones = nullptr;
    cont.x = nullptr;

    //Allocate Host Memory
    cont.x = (AB_type*)malloc(n*n*sizeof(AB_type));
    
    if(!cont.x)
    {
        printf("Host memory allocation failed!\n");
        destroy_context(cont);
        return EXIT_FAILURE;
    }
    

    //Create Cublas Handle
    stat = cublasCreate(&cont.handle);
    if (stat != CUBLAS_STATUS_SUCCESS) 
    {
        printf ("CUBLAS initialization failed\n");
        destroy_context(cont);
        return EXIT_FAILURE;
    }

    //Allocate Device Memory
    cudaStat = cudaMalloc((void**)&cont.devPtrX, n*n*sizeof(AB_type));
    if(cudaStat != cudaSuccess)exit(1);
    cudaStat = cudaMalloc((void**)&cont.devPtrY1, n*n*sizeof(AB_type));
    if(cudaStat != cudaSuccess)exit(1);
    cudaStat = cudaMalloc((void**)&cont.devPtrY2, n*n*sizeof(AB_type));
    if(cudaStat != cudaSuccess)exit(1);
    cudaStat = cudaMalloc((void**)&cont.devPtrY3, n*n*sizeof(AB_type));
    if(cudaStat != cudaSuccess)exit(1);
    cudaStat = cudaMalloc((void**)&cont.devPtrones, sizeof(AB_type));
    if(cudaStat != cudaSuccess)exit(1);
    
    if(cudaStat != cudaSuccess)
    {
        printf("Device memory allocation failed!\n");
        destroy_context(cont);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

template<typename AB_type, typename C_type>
void destroy_context(CublasContext<AB_type, C_type>& cont)
{
  double t=omp_get_wtime();
    if(cont.devPtrX) cudaFree(cont.devPtrX);
    if(cont.devPtrY1) cudaFree(cont.devPtrY1);
    if(cont.devPtrY2) cudaFree(cont.devPtrY2);
    if(cont.devPtrY3) cudaFree(cont.devPtrY3);
    if(cont.devPtrones) cudaFree(cont.devPtrones);
    if(cont.x) free(cont.x);
    printf("t dealloc=%f\n",omp_get_wtime()-t);
    
    if(cont.handle) cublasDestroy(cont.handle);
    printf("t destroy handle=%f\n",omp_get_wtime()-t);
}

void scale_matrix(int size, double* matrix)
{
    double frob=0;
    for(int i=0;i<size*size;i++){
      frob+=matrix[i]*matrix[i];
    }
    frob=sqrt(frob);
    
    double gersh=0;
    for(int j=0;j<size;j++)
    {
      double sum=0;
      for(int i=0;i<size;i++)
      {
        sum+=fabs(matrix[i*size+j]);
      }
      if(gersh<sum)gersh=sum;
    }

    double sc=1.0/fmin(frob,gersh);
    printf("scale=%f frob=%f gersh=%f\n",sc,frob,gersh);
    for(int i=0;i<size*size;i++){
      matrix[i]=matrix[i]*sc;
    }


}

// S = sgn(A), A = X_0
extern "C" int sign_gpu_(int* n_ptr, double* A, double* S, double* conv_ptr, int* order_ptr)
{
    int n=*n_ptr;
    int order=*order_ptr;
    double conv=*conv_ptr;
    int totmult=0;
    
    //scale_matrix(n,A);

    double tstart=omp_get_wtime();
    double ti=tstart;

    
    CublasContext<__type1, __type1> cont;
    cudaError_t cudaStat;
    cublasStatus_t stat;

    create_context(cont, n);
    
    printf("dim=%d\n",n);
    printf("Tinit context=%f\n",omp_get_wtime()-ti);
    ti=omp_get_wtime();

    //type casting input A to half
    //float* Af = (float*)malloc(n*n*sizeof(float));
	
#pragma omp parallel for schedule(static)
    for(int i = 0;i < n*n;i++)
    {
        cont.x[i] = A[i];
    }
    ti=omp_get_wtime();
    printf("Tinit convert1=%f\n",omp_get_wtime()-ti);
    
    
    printf("Tinit convert=%f\n",omp_get_wtime()-tstart);

    stat = cublasSetMatrix (n, n, sizeof(*cont.x), cont.x, n, cont.devPtrX, n); //column major format required
    if (stat != CUBLAS_STATUS_SUCCESS) 
    {
        printf ("Data download failed for A");
        destroy_context(cont);
        return EXIT_FAILURE;
    }

    __type1* ones = (__type1*)malloc(1*sizeof(__type1));
    ones[0] = (__type1)-1.0;

    stat = cublasSetVector (1, sizeof(*ones), ones, 1, cont.devPtrones, 1); //column major format required
    if (stat != CUBLAS_STATUS_SUCCESS) 
    {
        printf ("Data download failed for A");
        destroy_context(cont);
        return EXIT_FAILURE;
    }
    cudaDeviceSynchronize();
    printf("Tinit copy=%f\n",omp_get_wtime()-tstart);


    // calculate!
    float alpha0f = 1.0;
    double alpha0d = 1.0;

    __type2 alpha1 = 1.0;
    __type2 beta1 = 0.0;
    __type2 alpha2 = 3.0;
    __type2 beta2 = -10.0;
    __type2 alpha3 = 1.0/8.0;
    __type2 beta3 = 15.0/8.0;
    
    __type2 alpha4 = -0.5;
    __type2 beta4 = 1.5;

    double conv_res = 0;
    __type1 nrm=0.0;
    __type1 nrm_base=0.0;


    printf("Tinit=%f\n",omp_get_wtime()-tstart);

    int mixed=0;
    cublasGemmAlgo_t algo=CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    if(sizeof(nrm)!=sizeof(alpha1))mixed=1;
    cudaDataType_t cutype=CUDA_R_16F;
    cudaDataType_t cutypecalc=CUDA_R_16F;
    cudaDataType_t cutypeaxpy=CUDA_R_32F;
    if(sizeof(cont.x[0])==2){
      //Enable Math Mode
      stat = cublasSetMathMode(cont.handle, CUBLAS_TENSOR_OP_MATH);
      if(stat != CUBLAS_STATUS_SUCCESS)
      {
        printf("Failed setting Math Mode!\n");
        destroy_context(cont);
        return EXIT_FAILURE;
      }
      if(mixed==0){
        cutype=CUDA_R_16F;
        cutypecalc=CUDA_R_16F;
        cutypeaxpy=CUDA_R_32F;
        algo=CUBLAS_GEMM_DEFAULT_TENSOR_OP;
      }else{
        cutype=CUDA_R_16F;
        cutypecalc=CUDA_R_32F;
        cutypeaxpy=CUDA_R_32F;
        algo=CUBLAS_GEMM_DEFAULT_TENSOR_OP;
      }
    }
    else if(sizeof(cont.x[0])==4){
      cutype=CUDA_R_32F;
      cutypecalc=CUDA_R_32F;
      cutypeaxpy=CUDA_R_32F;
        algo=CUBLAS_GEMM_DEFAULT;
    }
    else if(sizeof(cont.x[0])==8){
      cutype=CUDA_R_64F;
      cutypecalc=CUDA_R_64F;
      cutypeaxpy=CUDA_R_64F;
      algo=CUBLAS_GEMM_DEFAULT;
    }

    printf("Start Iteration!\n");
    for(int iter=0;iter<order;iter++)
    {
        double t1=omp_get_wtime();
        int nmult=0;
        bool ischeckiter=(iter+1)%1==0;

        //Y=X^2
        nmult++;
        stat= cublasGemmEx(cont.handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha1, cont.devPtrX,cutype,n,cont.devPtrX,cutype,n,&beta1, cont.devPtrY1,cutype,n,cutypecalc,algo);
        //stat=cublasHgemm(cont.handle,CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha1, cont.devPtrX,n,cont.devPtrX,n,&beta1, cont.devPtrY1,n);
        if(stat!=CUBLAS_STATUS_SUCCESS)exit(1);
        cudaMemcpy(cont.devPtrY2,cont.devPtrX,n*n*sizeof(*cont.x),cudaMemcpyDeviceToDevice);
       
        if(1==1){
          nmult++;
          stat= cublasGemmEx(cont.handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha2, cont.devPtrX,cutype,n,cont.devPtrY1,cutype,n,&beta2, cont.devPtrY2,cutype,n,cutypecalc,algo);
          if(stat!=CUBLAS_STATUS_SUCCESS)exit(1);
          
          nmult++;
          stat= cublasGemmEx(cont.handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha3, cont.devPtrY2,cutype,n,cont.devPtrY1,cutype,n,&beta3, cont.devPtrX,cutype,n,cutypecalc,algo);
          if(stat!=CUBLAS_STATUS_SUCCESS)exit(1);
        }else{
          nmult++;
          stat= cublasGemmEx(cont.handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha4, cont.devPtrX,cutype,n,cont.devPtrY1,cutype,n,&beta4, cont.devPtrX,cutype,n,cutypecalc,algo);
          if(stat!=CUBLAS_STATUS_SUCCESS)exit(1);
        }

        if(ischeckiter){
          cudaMemcpy(cont.devPtrY3,cont.devPtrY1,n*n*sizeof(*cont.x),cudaMemcpyDeviceToDevice);

          stat=cublasNrm2Ex(cont.handle,n*n,cont.devPtrY3,cutype,1,&nrm_base,cutype,cutypeaxpy);
          if(stat!=CUBLAS_STATUS_SUCCESS){printf("cublasNrm2Ex_2 %d",stat);exit(1);}

          if(cutypeaxpy==CUDA_R_32F){
            stat=cublasAxpyEx(cont.handle,n,&alpha0f,cutypeaxpy,cont.devPtrones,cutype,0,cont.devPtrY3,cutype,n+1,cutypeaxpy);
          }else{
            stat=cublasAxpyEx(cont.handle,n,&alpha0d,cutypeaxpy,cont.devPtrones,cutype,0,cont.devPtrY3,cutype,n+1,cutypeaxpy);
          }
          if(stat!=CUBLAS_STATUS_SUCCESS){printf("cublasAxpyEx %d",stat);exit(1);}
          
          stat=cublasNrm2Ex(cont.handle,n*n,cont.devPtrY3,cutype,1,&nrm,cutype,cutypeaxpy);
          if(stat!=CUBLAS_STATUS_SUCCESS){printf("cublasNrm2Ex_1 %d",stat);exit(1);}
        }


        if(ischeckiter){
          printf("Iteration = %i, convergence = %f, %f\n", iter, (double)nrm/(double)nrm_base,pow(((double)nrm)/((double)nrm_base),2));
          if((double)nrm*(double)nrm < conv*(double)nrm_base*(double)nrm_base){
            break;
          }
        }

        if(stat != CUBLAS_STATUS_SUCCESS && cudaStat != cudaSuccess)
        {
            printf("iteration failed!\n");
            destroy_context(cont);
            return EXIT_FAILURE;
        }   
      
        cudaDeviceSynchronize();
        double t2=omp_get_wtime();
        double flops=nmult*2.0*pow(n,3)*1e-9;
        printf("titer=%f GFLOPS=%f GFLOPS/s=%f \n",t2-t1,flops,flops/(t2-t1));
        totmult+=nmult;

  }
  double tend=omp_get_wtime();
  //get current sign approximation         
  stat= cublasGetMatrix(n, n, sizeof(*cont.x), cont.devPtrX, n, cont.x, n);
  if(stat!=CUBLAS_STATUS_SUCCESS)exit(1);
  printf("Tcopy=%f\n",omp_get_wtime()-tend);

  ti=omp_get_wtime();
#pragma omp parallel for schedule(static)
  for(int i = 0;i < n*n; i++)S[i] = (double)cont.x[i];
  printf("Tend convert naiv=%f\n",omp_get_wtime()-ti);
  ti=omp_get_wtime();

  destroy_context(cont);

  printf("Tend=%f\n",omp_get_wtime()-tend);
  tend=omp_get_wtime();
  double flops=totmult*2.0*pow(n,3)*1e-9;
  printf("total T=%f GFLOPS=%f GFLOPS/s=%f \n",tend-tstart,flops,flops/(tend-tstart));
  
  return EXIT_SUCCESS;
}
