        module sm_module
        integer(4)    :: nb,nbmax
        real(8),allocatable :: blk(:,:,:)
        integer(4),allocatable    :: indexi(:)
        integer(4),allocatable    :: indexj(:)
        real(8)        :: basis(3)
        end module

        subroutine initseed()
        INTEGER :: i, n
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        real(8) :: svar
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
                               
        seed = 42 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)
        call random_number(svar)
        end subroutine

        program test
        call readmat()
        end program

        subroutine readmat()
        use omp_lib
        use sm_module
        !use RCIMAGEIO
        implicit none
        integer(4)  :: nrep,bs
        integer(4)    :: i,j,k,l,iatom,jatom,cols
        logical(4)  :: file_exists
        character(len=100) :: str,str2
        CHARACTER(len=100) :: arg
        CHARACTER(len=100) :: path
        integer(4),allocatable  :: meta(:,:),sm(:,:)
        integer(4)    :: col,smdim,is,js,bdim2,b2,smdim2,bdim,reason,mindim,nsuper
        real(8),allocatable   :: smat(:,:)
        real(8),allocatable   :: xyz(:,:),xyzd(:,:)
        real(8)               :: x,y,z,sx,sy,sz,svar,minops,mind
        real(8)               :: svar1,svar2,svar3
        integer(4),allocatable :: insm(:),insmall(:),insmmin(:),insmtmp(:)
        real(8)        :: dmax,tsingle,sglob,tmulti

        integer(4)    :: nsm,ierr,nsm2
        real(8),allocatable :: kmeans_work(:),kmeans_C(:,:)
        integer(4),allocatable :: kmeans_z(:),smdimtmp(:)


        CALL get_command_argument(1, path)
        CALL get_command_argument(2, arg)
        read(arg,*)nrep
        CALL get_command_argument(3, arg)
        read(arg,*)bs
        
        bdim=32*nrep**3
        
        allocate(meta(bdim,bdim))
        meta=0

        print*,"running for",path,"with NREP=",nrep," and bs=",bs
        
        print*,"reading strcture of matrix"
        write(str,'(A,A)')trim(path),"/blocks.dat"
        OPEN(42,file=trim(str))
        do i=1,bdim*bdim
          read(42,*,IOSTAT=reason)str,str2,iatom,jatom
          if(reason>0)then
            stop
          else if(reason<0)then
            exit
          else
            meta(iatom,jatom)=1
            meta(jatom,iatom)=1
          endif
        enddo
        close(42)
        !print*,"sum(meta)=",sum(meta)
        print*,"sparsity of full matrix (blocks)",sum(meta)
       
        allocate(insmall(bdim))
        allocate(insm(bdim))
        if(0.eq.1)then
          !write graph for metis
          call EXECUTE_COMMAND_LINE("rm -rf graph")
          OPEN(42,file='graph',status='new',action='write')
          write(42,*)bdim,(sum(meta)-bdim)/2
          do i=1,bdim
            do j=1,bdim
              if(meta(j,i).eq.1.and.i.ne.j)then
                write(42,'(I8)', advance="no")j
              endif
            enddo
            write(42,*)" "
          enddo
          close(42)
        
          insm(:)=0
          col=1!int(32*nrep**3*0.5)
          insm(col)=1
          call submatrix_smdim_multi(bdim,meta,insm,smdim)
          tsingle=(real(smdim,kind=8)**3)

          !do nsm=2,600,50
          !do nsm=2,int(bdim*0.3),int(real(bdim,kind=8)/100.0)
          do nsm=2,100
            write(str,'(A,I8)')'bash metis.sh ',nsm
            call EXECUTE_COMMAND_LINE(str)
            OPEN(42,file='graph.out')
            do i=1,bdim
              read(42,*)insmall(i)
            enddo
            close(42)
            
            sglob=0
            do j=1,nsm
              insm=0
              do i=1,bdim
                if(insmall(i)+1.eq.j)then
                  insm(i)=1
                endif
              enddo
              call submatrix_smdim_multi(bdim,meta,insm,smdim)
              svar=real(smdim,kind=8)**3/real(sum(insm),kind=8)/tsingle
!              print*,"metis",nsm,j,svar,smdim,sum(insm)
              sglob=sglob+svar
            enddo
            print*,"metis performance",nsm,sglob/nsm
          enddo
          deallocate(insmall)
        endif

        col=1
        !single column
        insm(:)=0
        insm(col)=1
        call submatrix_smdim_multi(bdim,meta,insm,smdim)
        tsingle=real(smdim,kind=8)**3
        
        allocate(sm(smdim,smdim))
        sm(:,:)=0
        call submatrix_create_block_multi(bdim,meta,insm,bs,smdim,sm)
        print*,"single column smdim=",col,smdim,Tsingle
        print*,"sparsity of submatrix (blocks)",sum(sm)
        deallocate(sm)

        if(1.eq.1)then        
          !build submatrix for meta
          col=1
          cols=1
          call submatrix_smdim(bdim,meta,col,cols,smdim)
          print*,"smdim",smdim

          allocate(smat(smdim*bs,smdim*bs))
          call submatrix_create(bdim,meta,col,cols,bs,smdim,smat)
          !measure sparsity of submatrix in elements
          k=0
          do i=1,smdim*bs
            do j=1,smdim*bs
              if(abs(smat(i,j)).gt.1e-5)then
                k=k+1
              endif 
            enddo
          enddo
          print*,"non-zero elements",k,smdim*bs
        
          !write as binary
          OPEN(10,file="ksmat.bin",access='stream')
          WRITE(10)smat
          CLOSE(10)
        endif
stop

        if(0.eq.1)then 
          !linear
          insm(:)=0
          do i=1,100
            insm(i)=1
            call submatrix_smdim_multi(bdim,meta,insm,smdim)
            allocate(sm(smdim,smdim))
            sm(:,:)=0
            call submatrix_create_block_multi(bdim,meta,insm,bs,smdim,sm)
            print*,"linear",col,i,smdim,real(smdim,kind=8)**3/real(sum(insm),kind=8)/tsingle,sum(insm)
            print*,"sparsity of submatrix (blocks)",sum(sm)/real(smdim,kind=8)**2,sum(sm)/real(smdim,kind=8)
            deallocate(sm)
          enddo
        endif
        
        
        !load spatial structure from xyz
        allocate(xyz(3,bdim))
        xyz=0
        !H2O-COORD.XYZ-pos-1.xyz
        !read file
        write(str,'(A,A)')trim(path),"/H2O-COORD.XYZ-pos-1.xyz"
        OPEN(42,file=trim(str))
        read(42,*)str
        read(42,*)str
        do iatom=1,bdim
          sx=0
          sy=0
          sz=0
          read(42,*)str,x,y,z
          sx=sx+x
          sy=sy+y
          sz=sz+z
          read(42,*)str,x,y,z
          sx=sx+x
          sy=sy+y
          sz=sz+z
          read(42,*)str,x,y,z
          sx=sx+x
          sy=sy+y
          sz=sz+z
          xyz(1,iatom)=sx/3.0D0
          xyz(2,iatom)=sy/3.0D0
          xyz(3,iatom)=sz/3.0D0
        enddo
        close(42)
        print*,(xyz(:,1))

        basis(1)=9.8528D0*nrep
        basis(2)=9.8528D0*nrep
        basis(3)=9.8528D0*nrep

        if(0.eq.1)then 
          !multi-column sm
          !do until all block columns are part of some submatrix
          allocate(insmall(bdim))
          allocate(insmmin(bdim))
          insmall=0
          nsm=0

          do while (sum(insmall).lt.bdim)
            minops=huge(minops)
            mind=0.0
            mindim=0
            !find first column that is not yet part of a submatrix
            do j=1,bdim
              if(insmall(j).eq.0)then
                col=j
                exit
              endif
            enddo

            do j=0,100
              insm(:)=0
              insm(col)=1
              dmax=j*0.1D0
              !find all block columns in radius
              do i=1,bdim
                call dist(xyz(:,col),xyz(:,i),basis,svar)
                if(svar.le.dmax)then
                  insm(i)=1
                endif
              enddo
              !remove columns that we already have in another submatrix
              do i=1,bdim
                if(insmall(i).eq.1.and.insm(i).eq.1)then
                  insm(i)=0
                endif
              enddo
              call submatrix_smdim_multi(bdim,meta,insm,smdim)
             
              if (j.eq.0)then
                tsingle=(real(smdim,kind=8)**3)
              endif

              svar=(real(smdim,kind=8)**3)/real(sum(insm),kind=8)/tsingle
              if (svar.lt.minops)then
                minops=svar
                mind=dmax
                insmmin=insm
                mindim=smdim
              endif
            enddo
            
            do j=1,bdim
              if(insmall(j).eq.0.and.insmmin(j).eq.1)then
                insmall(j)=1
              endif
            enddo
            print*,"multi column smdim: refcol=",col,"ncols=",sum(insmmin)," smdim=",mindim,"minops=",minops,mind
            sglob=sglob+minops
            nsm=nsm+1
          enddo
          print*,"multi column performance",sglob/nsm
          deallocate(insmall)
          deallocate(insmmin)
        endif

        insm(:)=0
        col=1!int(32*nrep**3*0.5)
        insm(col)=1
        call submatrix_smdim_multi(bdim,meta,insm,smdim)
        tsingle=(real(smdim,kind=8)**3)

        !do nsm=1,bdim,int(real(bdim,kind=8)/200.0)
        do nsm=2,int(bdim*0.3),int(real(bdim,kind=8)/100.0)
          if(1.eq.1)then
            !unperiodic kmeans
            sglob=0
!            basis(:)=huge(svar)
            allocate(kmeans_z(bdim))
            allocate(kmeans_work(bdim))
            allocate(kmeans_C(3,nsm))
            call initseed()
            CALL KMPP (xyz, 3, bdim, nsm, kmeans_c, kmeans_Z, kmeans_WORK, I)
!          print*,"kmeans_fault",i
!          do j=1,bdim
!            print*,"kmeans",j,kmeans_Z(j)
!          enddo 
            do j=0,nsm-1
              insm(:)=0
              do i=1,bdim
                if(j.eq.kmeans_z(i))then
                  insm(i)=1
                endif
              enddo

              call submatrix_smdim_multi(bdim,meta,insm,smdim)
              if(sum(insm).gt.0)then
                svar=(real(smdim,kind=8)**3)/real(sum(insm),kind=8)/tsingle
!              print*,"kmeans column smdim: ncols=",sum(insm)," smdim=",smdim,"minops=",svar
                sglob=sglob+svar
              endif
            enddo
            
            !number of block comuns combined for submatrix of first
            !molecule
            insm(:)=0
            do i=1,bdim
              if(kmeans_z(col).eq.kmeans_z(i))then
              !if(kmeans_z(1).eq.kmeans_z(i))then
                insm(i)=1
              endif
            enddo
            call submatrix_smdim_multi(bdim,meta,insm,smdim)
            svar=(real(smdim,kind=8)**3)/real(sum(insm),kind=8)/tsingle

            print*,"unperiodic kmeans fortran performance",nsm,sglob/nsm,sum(insm),svar
            deallocate(kmeans_z,kmeans_c,kmeans_work)
          endif
          
          if(0.eq.1)then
            !1d periodic kmeans
            sglob=0
!            basis(:)=huge(svar)

            nsuper=1
            allocate(kmeans_z(bdim*(2*nsuper+1)))
            allocate(kmeans_work(bdim*(2*nsuper+1)))
            allocate(kmeans_C(3,(2*nsuper+1)*nsm))
            allocate(xyzd(3,bdim*(2*nsuper+1)))
            j=0
            do i=1,bdim
              xyzd(1,1+j*bdim:bdim*(j+1))=xyz(1,1:bdim)+j*basis(1)
              xyzd(2,1+j*bdim:bdim*(j+1))=xyz(2,1:bdim)
              xyzd(3,1+j*bdim:bdim*(j+1))=xyz(3,1:bdim)
            enddo
            do j=1,nsuper
              do i=1,bdim
                xyzd(1,1+j*bdim:bdim*(j+1))=xyz(1,1:bdim)+j*basis(1)
                xyzd(2,1+j*bdim:bdim*(j+1))=xyz(2,1:bdim)
                xyzd(3,1+j*bdim:bdim*(j+1))=xyz(3,1:bdim)
              enddo
            enddo
            do j=1,nsuper
              do i=1,bdim
                xyzd(1,1+(nsuper+j)*bdim:bdim*(j+nsuper+1))=xyz(1,1:bdim)-j*basis(1)
                xyzd(2,1+(nsuper+j)*bdim:bdim*(j+nsuper+1))=xyz(2,1:bdim)
                xyzd(3,1+(nsuper+j)*bdim:bdim*(j+nsuper+1))=xyz(3,1:bdim)
              enddo
            enddo

            call initseed()
            CALL KMPP (xyzd, 3, (2*nsuper+1)*bdim, (2*nsuper+1)*nsm, kmeans_c, kmeans_Z, kmeans_WORK, IERR)
            deallocate(xyzd)
!          print*,"kmeans_fault",i
!          do j=1,bdim
!            print*,"kmeans",j,kmeans_Z(j)
!          enddo 
            nsm2=0
            do j=1,(2*nsuper+1)*nsm
              insm(:)=0
              do i=1,bdim
                if(j.eq.kmeans_z(i))then
                  insm(i)=1
                endif
              enddo
              if(sum(insm).gt.0)then
                call submatrix_smdim_multi(bdim,meta,insm,smdim)
                svar=(real(smdim,kind=8)**3)/real(sum(insm),kind=8)/tsingle
!              print*,"kmeans column smdim: ncols=",sum(insm)," smdim=",smdim,"minops=",svar
                sglob=sglob+svar
                nsm2=nsm2+1
              endif
            enddo
            print*,"1d periodic kmeans fortran performance",nsm,sglob/nsm2,IERR
            deallocate(kmeans_z,kmeans_c,kmeans_work)
          endif
         
          if(0.eq.1)then 
            !unperiodic kmeans
            sglob=0
            basis(:)=huge(svar)

            allocate(kmeans_z(bdim))
            !write to file
            write(100,'(A,A,A,I8,A,I8,A)')"python3 clustering.py ",trim(path),"/H2O-COORD.XYZ-pos-1.xyz",&
              nrep," ",nsm," -1 | tee | grep result > result"
            call EXECUTE_COMMAND_LINE("bash fort.100")

            OPEN(42,file='result')
            do j=1,bdim
              read(42,*)str,i,svar1,svar2,svar3,kmeans_z(j)
              kmeans_z(j)=kmeans_z(j)+1
            enddo
            close(42)

            do j=1,nsm
              insm(:)=0
              do i=1,bdim
                if(j.eq.kmeans_z(i))then
                  insm(i)=1
                endif
              enddo
              call submatrix_smdim_multi(bdim,meta,insm,smdim)
              if(sum(insm).gt.0)then
                svar=(real(smdim,kind=8)**3)/real(sum(insm),kind=8)/tsingle
!              print*,"kmeans column smdim: ncols=",sum(insm)," smdim=",smdim,"minops=",svar
                sglob=sglob+svar
              endif
            enddo
            print*,"unperiodic kmeans scikit performance",nsm,sglob/nsm
            deallocate(kmeans_z)
          endif

          if(0.eq.1)then
            !unperiodic hierarchical
            sglob=0
            basis(:)=huge(svar)

            allocate(kmeans_z(bdim))
            !write to file
            write(100,'(A,A,A,I8,A,I8,A)')"python3 clustering.py ",trim(path),"/H2O-COORD.XYZ-pos-1.xyz",&
              nrep," ",nsm," 0 | tee | grep result > result"
            call EXECUTE_COMMAND_LINE("bash fort.100")

            OPEN(42,file='result')
            do j=1,bdim
              read(42,*)str,i,svar1,svar2,svar3,kmeans_z(j)
              kmeans_z(j)=kmeans_z(j)+1
            enddo
            close(42)

            do j=1,nsm
              insm(:)=0
              do i=1,bdim
                if(j.eq.kmeans_z(i))then
                  insm(i)=1
                endif
              enddo
              call submatrix_smdim_multi(bdim,meta,insm,smdim)
              if(sum(insm).gt.0)then
                svar=(real(smdim,kind=8)**3)/real(sum(insm),kind=8)/tsingle
!                print*,"kmeans column smdim: ncols=",sum(insm)," smdim=",smdim,"minops=",svar
                sglob=sglob+svar
              endif
            enddo
            print*,"unperiodic hierarchical performance",nsm,sglob/nsm
            deallocate(kmeans_z)
          endif
          
          if(0.eq.1)then
            !periodic hierarchical
            sglob=0
            basis(:)=huge(svar)

            allocate(kmeans_z(bdim))
            !write to file
            write(100,'(A,A,A,I8,A,I8,A)')"python3 clustering.py ",trim(path),"/H2O-COORD.XYZ-pos-1.xyz",&
              nrep," ",nsm," 1 | tee | grep result > result"
            call EXECUTE_COMMAND_LINE("bash fort.100")

            OPEN(42,file='result')
            do j=1,bdim
              read(42,*)str,i,svar1,svar2,svar3,kmeans_z(j)
              kmeans_z(j)=kmeans_z(j)+1
            enddo
            close(42)

            do j=1,nsm
              insm(:)=0
              do i=1,bdim
                if(j.eq.kmeans_z(i))then
                  insm(i)=1
                endif
              enddo
              call submatrix_smdim_multi(bdim,meta,insm,smdim)
              if(sum(insm).gt.0)then
                svar=(real(smdim,kind=8)**3)/real(sum(insm),kind=8)/tsingle
!              print*,"kmeans column smdim: ncols=",sum(insm)," smdim=",smdim,"minops=",svar
                sglob=sglob+svar
              endif
            enddo
            print*,"periodic hierarchical performance",nsm,sglob/nsm
            deallocate(kmeans_z)
          endif
        enddo

        if(0.eq.1)then
          col=1
          print*,"start greedy",col
          !greedy for the first block column
          allocate( insmall(bdim))
          insmall(:)=0
          allocate(smdimtmp(bdim))

          nsm=0
          do i=1,30
            smdimtmp(:)=huge(i)

            !$OMP parallel do private(insm)
            do j=1,bdim
              insm=insmall
              if(insm(j).eq.1)then
                cycle
              endif
              insm(j)=1
              call submatrix_smdim_multi(bdim,meta,insm,smdimtmp(j))
            enddo

            smdim=huge(smdim)
            b2=0
            do j=1,bdim
              if(smdimtmp(j).lt.smdim)then
                b2=j
                smdim=smdimtmp(j)
              endif
            enddo

            insmall(b2)=1
            call submatrix_smdim_multi(bdim,meta,insmall,smdim)
            print*,"greedy",col,i,b2,smdim,real(smdim,kind=8)**3/real(sum(insmall),kind=8)/tsingle,sum(insmall)
            !allocate(sm(smdim,smdim))
            !sm(:,:)=0
            !call submatrix_create_block_multi(bdim,meta,insm,bs,smdim,sm)
            !print*,"sparsity of submatrix (blocks)",sum(sm)/real(smdim,kind=8)**2,sum(sm)/real(smdim,kind=8)
            !deallocate(sm)
          enddo
          deallocate(smdimtmp,insmall)
        endif
       

        if(0.eq.1)then
          col=1
          print*,"start greedy all",col
          allocate(insmall(bdim))
          insmall(:)=0
          insmall(col)=1
          allocate(smdimtmp(bdim))
          allocate( insmtmp(bdim))
          allocate(insmmin(bdim))

          sglob=0

          nsm=0
          do while(sum(insmall).lt.bdim)
            print*,sum(insmall)
            do i=1,bdim
              if(insmall(i).eq.0) then
                col=i
                exit
              endif
            enddo

            insmtmp=0
            insmtmp(col)=1
            call submatrix_smdim_multi(bdim,meta,insmtmp,smdim)
            tsingle=(real(smdim,kind=8)**3)

            nsm=nsm+1
            insm=0
            insm(col)=1

            tmulti=huge(tmulti)
            do i=0,min(50,bdim-sum(insmall))
!            do i=0,int(bdim/40.0),int(bdim/40.0)-1
              if(i.ne.0)then
                smdimtmp(:)=huge(i)
                
                !$OMP parallel do private(insmtmp)
                do j=1,bdim
                  insmtmp=insm
                  if(insmall(j).eq.1.or.insm(j).eq.1)then
                    cycle
                  endif
                  insmtmp(j)=1
                  call submatrix_smdim_multi(bdim,meta,insmtmp,smdimtmp(j))
                enddo

                smdim=huge(smdim)
                b2=0
                do j=1,bdim
                  if(smdimtmp(j).lt.smdim)then
                    b2=j
                    smdim=smdimtmp(j)
                  endif
                enddo
                if(b2.ne.0)insm(b2)=1
              else
                call submatrix_smdim_multi(bdim,meta,insm,smdim)
              endif
              svar=real(smdim,kind=8)**3/real(sum(insm),kind=8)/tsingle

              print*,"greedy step",nsm,i,b2,svar,sum(insm)
              if(svar.lt.tmulti)then
                tmulti=svar
                insmmin=insm
              endif
            enddo

            do i=1,bdim
              insmall(i)=max(insmall(i),insmmin(i))
            enddo
            print*,"greedy steps result",nsm,tmulti,sum(insmmin)

            call submatrix_smdim_multi(bdim,meta,insmmin,smdim)
            svar=real(smdim,kind=8)**3/real(sum(insmmin),kind=8)/tsingle
            print*,"greedy cluster",nsm,svar,sum(insmmin),sum(insmall)
            sglob=sglob+svar
          enddo
          print*,"total",sglob/nsm,nsm
          deallocate(smdimtmp,insmall)
        endif
        
        if(1.eq.1)then
          col=1
          print*,"start greedy individual"
          allocate(insmall(bdim))
!          insmall(:)=0
!          insmall(col)=1
          allocate(smdimtmp(bdim))
          allocate( insmtmp(bdim))
          allocate(insmmin(bdim))

          sglob=0

          nsm=0
          do col=1,32

            insmtmp=0
            insmtmp(col)=1
            call submatrix_smdim_multi(bdim,meta,insmtmp,smdim)
            tsingle=(real(smdim,kind=8)**3)

            nsm=nsm+1
            insm=0
            insm(col)=1

            tmulti=huge(tmulti)
            do i=0,50
!            do i=0,int(bdim/40.0),int(bdim/40.0)-1
              if(i.ne.0)then
                smdimtmp(:)=huge(i)
                
                !$OMP parallel do private(insmtmp)
                do j=1,bdim
                  insmtmp=insm
                  if(insm(j).eq.1)then
                    cycle
                  endif
                  insmtmp(j)=1
                  call submatrix_smdim_multi(bdim,meta,insmtmp,smdimtmp(j))
                enddo

                smdim=huge(smdim)
                b2=0
                do j=1,bdim
                  if(smdimtmp(j).lt.smdim)then
                    b2=j
                    smdim=smdimtmp(j)
                  endif
                enddo
                if(b2.ne.0)insm(b2)=1
              else
                call submatrix_smdim_multi(bdim,meta,insm,smdim)
              endif
              svar=real(smdim,kind=8)**3/real(sum(insm),kind=8)/tsingle

!              print*,"greedy step",nsm,i,b2,svar,sum(insm)
              if(svar.lt.tmulti)then
                tmulti=svar
                insmmin=insm
              endif
            enddo

!            do i=1,bdim
!              insmall(i)=max(insmall(i),insmmin(i))
!            enddo
            print*,"greedy steps result",nsm,tmulti,sum(insmmin)

!            call submatrix_smdim_multi(bdim,meta,insmmin,smdim)
!            svar=real(smdim,kind=8)**3/real(sum(insmmin),kind=8)/tsingle
!            print*,"greedy cluster",nsm,svar,sum(insmmin),sum(insmall)
!            sglob=sglob+svar
          enddo
!          print*,"total",sglob/nsm,nsm
          deallocate(smdimtmp,insmall)
        endif


        stop
        

        
        








stop

        end subroutine


        subroutine dist(x1,x2,b,d)
        implicit none
        real(8),intent(in)    :: x1(3)
        real(8),intent(in)    :: x2(3)
        real(8),intent(in)    :: b(3)
        real(8),intent(out)   :: d
        real(8)               :: bd(3)
        integer(4)    :: i,j,k

        d=huge(d)
        do i=-1,1
          do j=-1,1
            do k=-1,1
              bd(1)=i*b(1)
              bd(2)=j*b(2)
              bd(3)=k*b(3)
              d=min(d,sqrt(sum((x1-x2+bd)**2)))
            enddo
          enddo
        enddo

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
!          print*,"progress",is,"of",smdim
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
        
        subroutine submatrix_smdim_multi(metadim,meta,insm,smdim)
        implicit none
        integer(4),intent(in)   :: metadim
        integer(4),intent(in)   :: meta(metadim,metadim)
        integer(4),intent(in)   :: insm(metadim)
        integer(4),intent(out)   :: smdim
        integer(4)    :: is,js,ib,jb,a,b,i,j,iatom,jatom

        smdim=0
        do i=1,metadim
          do j=1,metadim
            if(insm(j).eq.1)then
              if(meta(i,j).eq.1)then
                smdim=smdim+1
                exit
              endif
            endif
          enddo
        enddo
        end subroutine
        
        subroutine submatrix_create_multi(metadim,meta,insm,bs,smdim,sm)
        implicit none
        integer(4),intent(in)   :: metadim,bs,smdim
        integer(4),intent(in)   :: meta(metadim,metadim),insm(metadim)
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
          do j=1,metadim
            if(insm(j).eq.1)then
              if(meta(i,j).eq.1)then
                smdim2=smdim2+1
                map(smdim2)=i
                imap(i)=smdim2
                exit
              endif
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

        subroutine submatrix_create_block_multi(metadim,meta,insm,bs,smdim,sm)
        implicit none
        integer(4),intent(in)   :: metadim,bs,smdim
        integer(4),intent(in)   :: meta(metadim,metadim),insm(metadim)
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
          do j=1,metadim
            if(insm(j).eq.1)then
              if(meta(i,j).eq.1)then
                smdim2=smdim2+1
                map(smdim2)=i
                imap(i)=smdim2
                exit
              endif
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
!          !overlap to refmap
!          ov=0
!          do i=1,smdim
!            do j=1,refsmdim
!              if(map(i).eq.refmap(j))then
!                ov=ov+1
!                exit
!              endif
!            enddo
!          enddo
!          print*,"sm",col,smdim,ov,int(ov/real(smdim,kind=8)*100.0D0)
        enddo
        end subroutine
