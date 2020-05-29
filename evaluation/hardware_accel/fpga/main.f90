      MODULE kinds
       IMPLICIT NONE

       PUBLIC :: real_4, real_8, real_4_size, real_8_size

       INTEGER, PARAMETER :: real_4 = SELECTED_REAL_KIND(6, 30)
       INTEGER, PARAMETER :: real_8 = SELECTED_REAL_KIND(14, 200)
       INTEGER, PARAMETER :: real_4_size = 4
       INTEGER, PARAMETER :: real_8_size = 8
      END MODULE
      
      module fpga_sgemm_mod
        use, intrinsic :: iso_c_binding
       interface
         subroutine multiply_c(m,n,k,A,B,C,test) bind(c,name="multiply")
            use, intrinsic :: iso_c_binding
            integer( kind=c_int32_t), value :: m,n,k,test
            REAL(KIND=c_float) :: A(m,k)
            REAL(KIND=c_float) :: B(k,n)
            REAL(KIND=c_float) :: C(m,n)
          end subroutine multiply_c
        end interface
      contains
        subroutine sgemm_fpga(TRANSA,TRANSB,M,N,K,alpha,A,LDA,B,ldb,beta,C,ldc)
        use kinds 
        use omp_lib
        implicit none
        CHARACTER,intent(in) :: TRANSA,TRANSB
        integer,intent(in)    :: M,N,K,lda,ldb,ldc
        real(KIND=real_4),intent(in),target      :: A(m,k)     
        real(KIND=real_4),intent(in),target      :: B(k,n)     
        real(KIND=real_8),intent(in)      :: alpha,beta
        real(KIND=real_4),intent(inout),target   :: C(m,n)
        real(KIND=real_4),target   :: D(m,n)
        real(kind=8)    :: t1,t2,t3,t4

        integer :: test=0
        t1=omp_get_wtime()

        if(TRANSA.eq.'t'.or.transa.eq.'T')then
          print*,"TransA=t/T not implemented"
          stop
        endif
        
        if(TRANSB.eq.'t'.or.transB.eq.'T')then
          print*,"TransB=t/T not implemented"
          stop
        endif


        if(LDA.ne.m)then
          print*,"LDA!=M not implemented"
          stop
        endif
        
        if(LDB.ne.k)then
          print*,"LDB!=k not implemented"
          stop
        endif
        
        if(LDC.ne.m)then
          print*,"LDC!=m not implemented"
          stop
        endif
        t2=omp_get_wtime()
        call multiply_c(m,n,k,A,B,D,test)
        t3=omp_get_wtime()

        call saxpby(m*k,real(alpha,kind=4),D,1,real(beta,kind=4),C,1)
!        C=real(beta,kind=4)*C+real(alpha,kind=4)*D
        t4=omp_get_wtime()
        
        print*,"TFlops_sgemm_fpga=",2.0D0*real(m,kind=8)**3*1d-12/(t3-t2),2.0D0*real(m,kind=8)**3*1d-12/(t4-t1)
        end subroutine sgemm_fpga
      end module

      program test
      use kinds 
      use omp_lib
      implicit none
      integer(4)            :: n,maxiter,i,bs,col=1,n2
      real(8),allocatable   :: A(:,:),S(:,:),Ascaled(:,:),S2(:,:)
      real(8)               :: conv,Eref,INV,work,frob_matrix,gersh_matrix
      CHARACTER(len=32) :: arg
      REAL(KIND=8), EXTERNAL                            :: dlange
    
      CALL get_command_argument(1, arg)
      read(arg,*) n
      CALL get_command_argument(2, arg)
      read(arg,*) bs
      CALL get_command_argument(3, arg)
      read(arg,*) conv
      CALL get_command_argument(4, arg)
      read(arg,*)maxiter

      print*,"Matrix size=",n
      print*,"Block size=",bs
      print*,"convergence crit=",conv

      allocate(A(n,n))

!      call makerandommat(n,A)
      OPEN(10,file="ksmat.bin",access='stream')
      READ(10)A
      CLOSE(10)

      if(modulo(n,8).eq.0)then
        n2=n
      else
        n2=n+8-modulo(n,8)
      endif
      print*,"padded size",n2
      
      allocate(Ascaled(n2,n2))
      allocate(S(n2,n2))
      Ascaled=0
      do i=1,n2
        Ascaled(i,i)=Ascaled(i,i)+1.0D0
      enddo
      Ascaled(1:n,1:n)=A(1:n,1:n)
      call scalemat(n2,Ascaled)

      !one iteration for gpu init
      call sign_fpga(n2,Ascaled,S,conv,1)

      call sign_fpga(n2,Ascaled,S,conv,maxiter)
      
      allocate(S2(n2,n2))
      CALL dgemm('N', 'N', n2, n2, n2, 1.0D0, S, n2, S, n2, 0.0D0,S2, n2)
      do i=1,N2
        S2(i,i)=S2(i,i)-1.0D0
      enddo
      frob_matrix  = dlange('F', n2, n2, S2, n2, work) !dbcsr_frobenius_norm(matrix_sign)
      gersh_matrix = dlange('1', n2, n2, S2, n2, work) !dbcsr_gershgorin_norm(matrix_sign)
      
      !from sign function to density matrix 
      S=-0.5D0*S
      do i=1,n
        S(i,i)=S(i,i)+0.5D0
      enddo
  
      Eref=sum(S(1:n,(col-1)*bs+1:col*bs)*A(1:n,(col-1)*bs+1:col*bs))
      print*,"E=",col,Eref
      print*,"INV=",col,frob_matrix,gersh_matrix

      end program

      subroutine makerandommat(n,A)
      implicit none
      integer(4),intent(in)        :: n
      real(8),intent(inout)        :: A(n,n)
      integer(4)                   :: i,j
      real(8)                      :: svar

      do i=1,n
        call random_number(svar)
        A(i,i)=2.0*svar-1.0
        do j=i+1,n
          call random_number(svar)
          A(i,j)=2.0*svar-1.0
          A(j,i)=A(i,j)

        enddo
      enddo
      call scalemat(n,A)
      end subroutine

      subroutine scalemat(n,A)
      implicit none
      integer(4),intent(in)   :: N
      real(8),intent(inout)   :: A(N,N)
      real(8)                 :: work
      real(8)                      :: frob_matrix,gersh_matrix,scaling_factor
      REAL(KIND=8), EXTERNAL                            :: dlange

      frob_matrix  = dlange('F', n, n, A, n, work) !dbcsr_frobenius_norm(matrix_sign)
      gersh_matrix = dlange('1', n, n, A, n, work) !dbcsr_gershgorin_norm(matrix_sign)
      scaling_factor = 1/MIN(frob_matrix, gersh_matrix)
      !print*,"norms",sz,frob_matrix,gersh_matrix,scaling_factor
      A = A * scaling_factor

      end subroutine

      subroutine sign_fpga(n,Xin,S,conv,maxiter)
      use kinds 
      use fpga_sgemm_mod
      use omp_lib
      implicit none
      integer(4),intent(in)   :: n,maxiter
      real(kind=real_8),intent(in) :: Xin(n,n)
      real(kind=real_8),intent(out) :: S(n,n)
      real(kind=8),intent(in)   :: conv
      integer(4)    :: i,j,nmult,iter
      real(4)   :: alpha=1.0D0,beta=0.0D0
      real(kind=real_8) :: diff,t1,t2,e0,e1
      real(kind=4),allocatable    :: X(:,:),Y1(:,:),Y2(:,:)
      REAL(KIND=4), EXTERNAL                            :: snrm2

      allocate(X(n,n))
      allocate(Y1(n,n))
      allocate(Y2(n,n))
      X=Xin

      t1=omp_get_wtime()
      nmult=0

      !third-order Pade sign iteration
      do iter=1,maxiter
        !call multiplication on fpga

        call sgemm_fpga('N','N',n,n,n,1.0D0,X,n,X,n,0.0D0,Y1,n)
        nmult=nmult+1

        call scopy(N*N,X,1,Y2,1)
        call sgemm_fpga('N','N',n,n,n,3.0D0,X,n,Y1,n,-10.0D0,Y2,n)
        nmult=nmult+1
        
        call sgemm_fpga('N','N',n,n,n,1.0D0/8.0D0,Y2,n,Y1,n,15.0D0/8.0D0,X,n)
        nmult=nmult+1

        e0=snrm2(n*n,Y1,1)
        do i=1,n
          Y1(i,i)=Y1(i,i)-1.0D0
        enddo
        e1=snrm2(n*n,Y1,1)
        print*,"iter",iter,e0,e1
        if(e1*e1.lt.e0*e0*conv)then
          exit
        endif

      enddo
      t2=omp_get_wtime()
      print*,"TFlops_sign=",2.0D0*real(n,kind=8)**3*nmult*1d-12/(t2-t1)
      S=X

      end subroutine

