      program test
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
      call sign_gpu(n2,Ascaled,S,conv,1)

      call sign_gpu(n2,Ascaled,S,conv,maxiter)
      
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
