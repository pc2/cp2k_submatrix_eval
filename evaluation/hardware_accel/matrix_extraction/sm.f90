        module sm
        integer(4)    :: nb,nbmax
        real(8),allocatable :: blk(:,:,:)
        integer(4),allocatable    :: indexi(:)
        integer(4),allocatable    :: indexj(:)
        end module

        program test
        call readmat(32*5**3,6)
        end program

        subroutine readmat(bdim,bs)
        use omp_lib
        use sm
        use RCIMAGEIO
        implicit none
        integer(4),intent(in)   :: bdim,bs
        integer(4)    :: i,j,k,l,iatom,jatom,cols
        logical(4)  :: file_exists
        character(len=100) :: str
        CHARACTER(len=100) :: arg
        integer(4),allocatable  :: meta(:,:)
        integer(4)    :: col,smdim,is,js,bdim2,b2,smdim2
        real(8),allocatable   :: smat(:,:)

        allocate(meta(bdim,bdim))
        meta=0

        CALL get_command_argument(1, arg)
        
        print*,"reading strcture of matrix"
        !$omp parallel do schedule(dynamic) private(jatom,str,file_exists)
        do iatom=1,bdim
          print*,iatom
          do jatom=1,iatom
            write(str,'(A,A,I0.8,A,I0.8)')trim(arg),"/mat/",iatom,"/",jatom
            INQUIRE(FILE=trim(str), EXIST=file_exists)
            if(file_exists)then
              meta(iatom,jatom)=1
            endif
          enddo
        enddo
        
        !build submatrix for meta
        call submatrix_meta(bdim,meta)

        col=1
        cols=32
        call submatrix_smdim(bdim,meta,col,cols,smdim)

        allocate(smat(smdim*bs,smdim*bs))
        call submatrix_create(bdim,meta,col,cols,bs,smdim,smat)
        
        !write as binary
        OPEN(10,file="ksmat.bin",access='stream')
        WRITE(10)smat
        CLOSE(10)

        end subroutine

        subroutine submatrix_create(metadim,meta,col,cols,bs,smdim,sm)
        implicit none
        integer(4),intent(in)   :: metadim,col,cols,bs,smdim
        integer(4),intent(in)   :: meta(metadim,metadim)
        real(8),intent(inout)      :: sm(smdim*bs,smdim*bs)
        integer(4)    :: map(metadim),smdim2
        integer(4)    :: imap(metadim),is,js,ib,jb,a,b,i,iatom,jatom,j
        real(8)       :: svar
        logical(4)  :: file_exists
        character(len=100) :: str
        CHARACTER(len=100) :: arg

        map=0
        imap=0

        smdim2=0
        do i=1,metadim
          do j=col,col+cols-1
            if(meta(i,j).eq.1)then
              smdim2=smdim2+1
              map(smdim2)=i
              imap(i)=smdim2
              exit
            endif
          enddo
        enddo
        CALL get_command_argument(1, arg)

        sm=0
        print*,"reading blocks"
        do is=1,smdim
          print*,"progress",is,"of",smdim
          iatom=map(is)
          do js=1,is !smdim
            jatom=map(js)
            if(meta(iatom,jatom).eq.0)cycle

            write(str,'(A,A,I0.8,A,I0.8)')trim(arg),"/matrix_ssqrtinv_ks_ssqrtinv/",iatom,"/",jatom
            INQUIRE(FILE=trim(str), EXIST=file_exists)
            if(file_exists)then
              !read file
              OPEN(42,file=trim(str))
                  
              do ib=1,bs
                do jb=1,bs
                  read(42,*)a,b,svar
                  sm((is-1)*bs+ib,(js-1)*bs+jb)=svar
                  sm((js-1)*bs+jb,(is-1)*bs+ib)=svar
                enddo
              enddo
              close(42)
            endif
          enddo
        enddo
        end subroutine
        
        subroutine submatrix_smdim(metadim,meta,col,cols,smdim)
        implicit none
        integer(4),intent(in)   :: metadim,col,cols
        integer(4),intent(in)   :: meta(metadim,metadim)
        integer(4),intent(out)   :: smdim
        integer(4)    :: map(metadim),j
        integer(4)    :: imap(metadim),is,js,ib,jb,a,b,i,iatom,jatom
        real(8)       :: svar
        logical(4)  :: file_exists
        character(len=100) :: str
        CHARACTER(len=100) :: arg

        smdim=0
        do i=1,metadim
          do j=col,col+cols-1
            if(meta(i,j).eq.1)then
              smdim=smdim+1
              exit
            endif
          enddo
        enddo
        end subroutine

        subroutine submatrix_create_block(metadim,meta,col,cols,bs,smdim,sm)
        implicit none
        integer(4),intent(in)   :: metadim,col,bs,smdim,cols
        integer(4),intent(in)   :: meta(metadim,metadim)
        integer(4),intent(inout)      :: sm(smdim,smdim)
        integer(4)    :: map(metadim),smdim2
        integer(4)    :: imap(metadim),is,js,ib,jb,a,b,i,iatom,jatom,j
        real(8)       :: svar
        logical(4)  :: file_exists
        character(len=100) :: str
        CHARACTER(len=100) :: arg

        map=0
        imap=0

        smdim2=0
        do i=1,metadim
          do j=col,col+cols-1
            if(meta(i,j).eq.1)then
              smdim2=smdim2+1
              map(smdim2)=i
              imap(i)=smdim2
              exit
            endif
          enddo
        enddo
        CALL get_command_argument(1, arg)

        sm=0
        do is=1,smdim
          iatom=map(is)
          do js=1,smdim
            jatom=map(js)
            if(meta(iatom,jatom).eq.1)then
              sm(is,js)=1
            endif
          enddo
        enddo
        end subroutine


        subroutine submatrix_meta(metadim,meta)
        implicit none
        integer(4),intent(in)   :: metadim
        integer(4),intent(inout)   :: meta(metadim,metadim)
        integer(4)  :: i,j,k,l,col,smdim,ov
        integer(4)    :: map(metadim)
        integer(4)    :: imap(metadim)
        integer(4)    :: refmap(metadim)
        integer(4)    :: refimap(metadim),refsmdim=0

        integer(4)    :: refcol=1
        
        smdim=0
        refmap=0
        refimap=0
        do i=1,metadim
          if(meta(i,refcol).eq.1)then
            refsmdim=refsmdim+1
            refmap(refsmdim)=i
            refimap(i)=refsmdim
          endif
        enddo

        do col=1,metadim
          smdim=0
          do i=1,metadim
            if(meta(i,col).eq.1)then
              smdim=smdim+1
              map(smdim)=i
              imap(i)=smdim
            endif
          enddo
          !overlap to refmap
          ov=0
          do i=1,smdim
            do j=1,refsmdim
              if(map(i).eq.refmap(j))then
                ov=ov+1
                exit
              endif
            enddo
          enddo
          print*,"sm",col,smdim,ov,int(ov/real(smdim,kind=8)*100.0D0)
        enddo
        end subroutine
