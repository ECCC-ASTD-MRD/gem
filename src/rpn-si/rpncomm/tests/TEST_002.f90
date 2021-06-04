!
!    Test verifying that rpn_comm_xch_halo and rpn_comm_transpose
!    are giving correct results. Furthermore, it returns a timing
!    for rpn_comm_xch_halo
!
        subroutine rpn_comm_test_002
	logical :: asynchronous  ! use asynchronous code
!        print *,'I am version ',VeRsion
        call rpn_comm_test_002b
        stop
        end
        subroutine get_n_domains(ndomains, offset, err)
        integer n_domains
        common/ndomains/n_domains
        integer ndomains, offset, err
        character (len=128) SYSTEM_COMMAND
	SYSTEM_COMMAND="1"
        call RPN_COMM_env_var("TEST_DOMAINS",SYSTEM_COMMAND)
        if(SYSTEM_COMMAND == "" )SYSTEM_COMMAND="1"
        read(SYSTEM_COMMAND,*)ndomains
        n_domains=ndomains
        offset = 0
        err = 0
        return
        end
        subroutine rpn_comm_test_002b
        implicit none
!
!       dummy test program implemented as a subroutne
!       because some compilers do not like POINTER
!       statements with variable dimensions in a main program
!
        integer :: nptsx,nptsy,nptsz,ihalox,ihaloy
        common /the_problem/ nptsx,nptsy,nptsz,ihalox,ihaloy  ! common containing problem dimensions and halo size
        integer :: nodex, nodey
        common /pernode/ nodex, nodey

        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
        integer, allocatable, dimension(:)  ::iarr,jarr
        integer, allocatable, dimension(:,:,:) ::data,data2,glob
!
        external sss
        integer npex,npey
        integer min3,max3,ierr, istat
        integer mini,maxi,nil,nilmax,ni0,i0
        integer minj,maxj,njl,njlmax,nj0,j0
        integer nzl, nzlmzx, nz0,irep,iter,irep2
        real*8 time1,time2,MPI_wtime,time_min,time_max,time_tot
        external MPI_wtime
!
        integer RPN_COMM_topo
        logical RPN_COMM_ngrank, RPN_COMM_grank
        integer RPN_COMM_init_multi_level, mygrid, mygridcomm, mymultigridcomm,myworldcomm,myallgridcomm
        integer RPN_COMM_colors
        integer RPN_COMM_comm
        integer test_grids
	integer RPN_COMM_option_L
        external RPN_COMM_mydomain
        external get_n_domains
        external test_grids
	external RPN_COMM_option_L
        integer domains, mydomain,irank,isize
        integer n_domains,j
        integer peercomm, npeers
        common/ndomains/n_domains
        character *6 is_async(2)

        integer cache_flush(2000*2000)
        integer, parameter :: NEXCH=256
        real, parameter :: SCAL=1.0/NEXCH
        real time_ew(0:5*NEXCH), time_ns(0:5*NEXCH)
        integer nbytes, nodebytes
!
	call RPN_COMM_xch_halo(time_ew,-1024,0,0,0,1,0,0,1,0,.true.,.true.,0,-1) ! fake call, set EW timing array
	call RPN_COMM_xch_halo(time_ns,-1024,0,0,0,1,0,0,0,1,.true.,.true.,0,-1) ! fake call, set NS timing array

!
!       force a "callback" type initialization
!
        npex=0
        npey=0
        call RPN_COMM_mydomain (get_n_domains, mydomain,ierr)
!        print *,'This is member',mydomain+1,' of',n_domains,' domains'
!
        call RPN_COMM_set_petopo(1,1000)   ! force vertically striped distribution
        mygrid = RPN_COMM_init_multi_level(sss,Pelocal,Petotal,npex,npey,n_domains,1)
!
!       ============= determine resolution of MPI_wtine
!
        time1 = MPI_WTIME()
        time_min=1.0
        time_max=0.0
        time_tot = 0.0
        do irep = 1 , 100
          time2 = MPI_WTIME()
          if(time2-time1 > 0)time_min=min(time_min,time2-time1)
          time_max=max(time_max,time2-time1)
          time_tot = time_tot + time2-time1
          time1 = time2
        enddo
        if (Pelocal.eq.0) then
          print *,'MPI_WTIME cost(min,max,avg)=',real(time_min),real(time_max),real(time_tot)*.01
        endif
!
!	============= TEST for WORLD/MULTIGRID/GRID ====================
!
        if(test_grids(mygrid,Pelocal) /= 0) goto 9999
!
!	================= TEST for halo exchanger ===================
!
        call RPN_COMM_BCAST(nptsx,5,'MPI_INTEGER',0,'GRID',ierr)  ! broadcast problem dimensions to GRID
!
!       get data mapping (should really make 2 calls to
!       the 1D function RPN_COMM_topo)
!
        if(.not.RPN_COMM_grank('GRID')) goto 2222    ! if this PE is not a member of 'GRID', bail out
!
        istat= RPN_COMM_topo(nptsx,mini,maxi,nil,nilmax,ihalox,ni0,.TRUE.,.FALSE.)
        if(istat.ne.0) then
           write(*,*) 'Invalid distribution over x, abort',nptsx,npex
           goto 9999
        endif
!
        istat= RPN_COMM_topo(nptsy,minj,maxj,njl,njlmax,ihaloy,nj0,.FALSE.,.FALSE.)
        if(istat.ne.0) then
           write(*,*) 'Invalid distribution over y, abort',nptsy,npey
           goto 9999
        endif
!
        if (Pelocal.eq.0) then  ! print dimensions from PE 0
!          write(*,*) 'I AM ',Pelocal
           print *,'nptsx,nptsy,nptsz,mini,maxi,minj,maxj=',nptsx,nptsy,nptsz,mini,maxi,minj,maxj
           print *,'nil,njl,ihalox,ihaloy=',nil,njl,ihalox,ihaloy
           print *,'nilmax,njlmax=',nilmax,njlmax
        endif
!
!     allocate local arrays
!
        allocate(iarr(mini:maxi))
        allocate(jarr(minj:maxj))
        allocate(data(mini:maxi,minj:maxj,nptsz*2))
        allocate(glob(1-ihalox:nptsx+ihalox,minj:maxj,nptsz))
!
!     fill arrays with markers
!
        i0=ni0
        j0=nj0
        call set_ijarr(iarr,jarr,mini,maxi,minj,maxj,i0,j0,nptsx,nptsy)
        call setup_arr(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl)
        glob = -1
        glob(1:nil,minj:maxj,1:nptsz) = data(1:nil,minj:maxj,1:nptsz)
        call RPN_COMM_xch_halosl(glob,mini,nptsx+ihalox,minj,maxj, &
              & nil,njl,nptsz,ihalox,ihaloy,.true.,.true.,nptsx,0,nilmax)
        if (Pelocal.eq.0 .and. nptsx < 10 .and. njl < 5) then 
          do j=maxj,minj,-1
            print 101,j,glob(:,j,1)
          enddo
        endif
        call vfy_arr_glb(glob,jarr,1-ihalox,nptsx+ihalox,minj,maxj,nptsz,nptsx,njl,ihalox)
!
!     timing loop for halo exchange
!
	is_async(1) = 'ASYNC'
	is_async(2) = 'SYNC'
	do iter = 1 , 2
	ierr = RPN_COMM_option_L('async_exch',iter == 1)
101     format(I3,15I9)
!        print *,'----------------'
        time1 = MPI_WTIME()
        time_min=1.0
        time_max=0.0
        time_tot = 0.0
        nbytes=2*(ihaloy*nil + (maxj-minj+1)*ihalox) * nptsz * 2 * 4
        nodebytes=2*(ihaloy*nodex + nodey*ihalox) * nptsz * 2 * 4
        call rpn_comm_wtime_set(MPI_Wtime)
        do irep=1,NEXCH ! NEXCH repetitions of fully periodic halo exchange
          cache_flush = cache_flush + 1
          glob(1:nil,minj:maxj,1:nptsz) = data(1:nil,minj:maxj,1:nptsz)
          call RPN_COMM_barrier('GRID',ierr)
          time1 = MPI_WTIME()
          call RPN_COMM_xch_halo(data,mini,maxi,minj,maxj,nil,njl,nptsz,ihalox,ihaloy,.true.,.true.,nptsx,0)
          time2 = MPI_WTIME()
          time_min=min(time_min,time2-time1)
          time_max=max(time_max,time2-time1)
          time_tot = time_tot + time2-time1
          time1 = time2
        enddo
!
!          call affichage(data,mini,maxi,minj,maxj,nptsz) 
!
        if(Pelocal.eq.0 )then
           print *,is_async(iter)//'time_tot (min,max,avg)=',real(time_min),real(time_max),real(time_tot)/NEXCH
           print *,'pe=',Pelocal,                                          &
     &             ' Number of exchanges=', irep-1,                           &
     &             ' Time per exchange=',nint(1000000*time_tot/(irep-1)),            &
     &             ' microseconds', nint(.000001 * nbytes * NEXCH / time_tot),' MB/s',  &
     &             ' node BW=',nint(.000001 * nodebytes * NEXCH / time_tot),' MB/s'
        endif
        if(Pelocal.eq.0 )then
            irep=time_ew(0)-1
            irep2=time_ns(0)-1
            call timing_report(irep/3,time_ew(1),is_async(iter)//'time_ew')
            call timing_report(irep2/3,time_ns(1),is_async(iter)//'time_ns')
!            print *,is_async(iter)//'time_ew',irep
!            print 111,'T1 ',nint(time_ew(1:irep:3)*1000000)
!            print 111,'T2 ',nint(time_ew(2:irep:3)*1000000 - time_ew(1:irep:3)*1000000)
!            print 111,'T3 ',nint(time_ew(3:irep:3)*1000000 - time_ew(2:irep:3)*1000000)
111         format(A,15I6)
!           print *,'------COPIES------'
!           print *,times(4,1:20)-times(3,1:20)+times(2,1:20)-times(1,1:20)
!           print *,'-------MPI-------'
!           print *,times(3,1:20)-times(2,1:20)
!           print *,'------COPIES------'
!           print *,times2(4,1:20)-times2(3,1:20)+times2(2,1:20)-times2(1,1:20)
!           print *,'-------MPI-------'
!           print *,times2(3,1:20)-times2(2,1:20)
        endif
	call RPN_COMM_xch_halo(time_ew,-1024,0,0,0,1,0,0,1,0,.true.,.true.,0,-1) ! fake call, reset EW timing array
	call RPN_COMM_xch_halo(time_ns,-1024,0,0,0,1,0,0,0,1,.true.,.true.,0,-1) ! fake call, reset NS timing array
	enddo
!
!     verify that there were no errors
!
        call vfy_arr(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
!
!	call RPN_COMM_FINALIZE(ierr)
!	STOP
!
!	================= TEST for transpose ===================
!
           call setup_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl)
           istat =  RPN_COMM_topo(nptsz,min3,max3,nzl,nzlmzx,0,nz0,.true.,.true.)
           allocate(data2((max3-min3+1),(maxj-minj+1),(nptsx+10)*2))
           if(.true.) then
             if(Pelocal.eq.0 )then
              print *,' size of data =',                                      &
     &             (maxi-mini+1),(maxj-minj+1),nptsz,                         &
     &             (maxi-mini+1)*(maxj-minj+1)*nptsz*2
              print *,' size of data2 =',                                     &
     &             (maxj-minj+1),(max3-min3+1),nptsx,                         &
     &             (max3-min3+1)*(maxj-minj+1)*(nptsx+10)*2
              print *,'nptsz,min3,max3,nzl,nzlmzx,nz0',                       &
     &             nptsz,min3,max3,nzl,nzlmzx,nz0
             endif
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,1,2)
!
             call vfy_xpos(data2,jarr,minj,maxj,njl,nz0,nzl,min3,max3,nptsx,j0)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,-1,2)
!
             call vfy_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,1,2)
!
             call vfy_xpos(data2,jarr,minj,maxj,njl,nz0,nzl,min3,max3,nptsx,j0)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,-1,2)
!
             call vfy_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
           endif
!
!     THE END !!
!
!
 2222   Continue
        call RPN_COMM_BARRIER('WORLD',ierr)
 9999	continue
        call RPN_COMM_FINALIZE(ierr)
        stop
        end
!
        subroutine vfy_xpos(z,jtab,minj,maxj,njl,nz0,nzl,min3,max3,nptsx,j0)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers
!
!
        integer minj,maxj,njl,min3,max3,nptsx,nzl,nz0,j0
        integer z(2,minj:maxj,min3:max3,nptsx)
        integer jtab(minj:maxj)
        integer i,j,k
        integer error
        error=0
        do i=1,nptsx
        do k=1,nzl
        do j=1,njl
            if(z(1,j,k,i)/100000 .ne. i) error=error+1
            if(mod(z(1,j,k,i), 100000)/100 .ne. j0+j-1) error=error+1
            if(z(2,j,k,i).ne.( nz0-1+k )) error = error + 1
        enddo
        enddo
        enddo
        if(error.ne.0) then 
          print *,error,' VERIFICATION ERRORS for Pe ',Pelocal
        endif
      if(Pelocal.eq.0) print *,'vfy_xpos : Verification performed, number of errors was ',error
        return
        end
!
        subroutine vfy_arr2(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers ignoring k
!
        integer minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy
        integer z(2,minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k,k0
        integer error
        error=0
        do k=1,nk
         k0=mod(k,3)
          do j=1,njl
          do i=1,nil
            if(z(1,i,j,k).ne.jtab(j)*100+itab(i)*100000) error = error + 1
            if(z(2,i,j,k).ne.( k ))  error = error + 1
          enddo
          enddo
        enddo
        if(error.ne.0) then 
          print 1,'vfy_arr2:',error,' VERIFICATION ERRORS for Pe ',Pelocal,' (',(nil)*(njl)*nk,' points)'
1       format(A,I4,A,I4,A,I5,A)
        endif
        if(Pelocal.eq.0) print *,'vfy_arr2 : Verification performed,',error,' errors on PE 0'
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(1,i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        subroutine vfy_arr_glb(z,jtab,minx,maxx,miny,maxy,nk,gni,njl,ihalox)
!!!!        use rpn_comm
        implicit none
        integer :: minx,maxx,miny,maxy,nk,gni,njl,ihalox
        integer, dimension(minx:maxx,miny:maxy,nk) :: z
        integer, dimension(miny:maxy) :: jtab

        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify global (along X) array containing markers
!
        integer i,j,k,ref,error,ii

        error = 0
        do k=1,nk
         do j=1,njl
          do i=minx,maxx
            ii = i
            if(ii < 1) ii = ii + gni
            if(ii > gni) ii = ii - gni
            ref = k
            ref = ref + jtab(j)*100
            ref = ref + ii*100000
            if(z(i,j,k).ne.ref) then
              error = error + 1
            endif
          enddo
         enddo
        enddo
        if(error.ne.0) then 
          print 1,'vfy_arr_glb:',error,' VERIFICATION ERRORS for Pe ',Pelocal,' (',(maxx-minx+1)*njl*nk,' points)'
1       format(A,I4,A,I4,A,I5,A)
        endif
        if(Pelocal.eq.0) print *,'vfy_arr_glb : Verification performed,',error,' errors on PE 0'
        return
        end
        subroutine vfy_arr(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers
!
        integer minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy
        integer z(minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer error,k0
        integer i,j,k
        integer ref
        error=0
        do k=1,nk
          do j=1-ihaloy,njl+ihaloy
          do i=1-ihalox,nil+ihalox
            ref = k
            ref = ref + jtab(j)*100
            ref = ref + itab(i)*100000
            if(z(i,j,k).ne.ref) error = error + 1
          enddo
          enddo
        enddo
        if(error.ne.0) then 
          print 1,'vfy_arr:',error,' VERIFICATION ERRORS for Pe ',Pelocal,' (',(nil+2*ihalox)*(njl+2*ihaloy)*nk,' points)'
1       format(A,I4,A,I4,A,I5,A)
        endif
        if(Pelocal.eq.0) print *,'vfy_arr : Verification performed,',error,' errors on PE 0'
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        subroutine set_ijarr(itab,jtab,minx,maxx,miny,maxy,i0,j0,nptsx,nptsy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        precompute initial position tables
!
        integer minx,maxx,miny,maxy,nptsx,nptsy,i0,j0
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j
        do i=minx,maxx
          itab(i)=i0+i-1
          if(itab(i).gt.nptsx) itab(i)=itab(i)-nptsx
          if(itab(i).le.0    ) itab(i)=itab(i)+nptsx
        enddo
        do j=miny,maxy
          jtab(j)=j0+j-1
          if(jtab(j).gt.nptsy) jtab(j)=jtab(j)-nptsy
          if(jtab(j).le.0    ) jtab(j)=jtab(j)+nptsy
        enddo
        return
        end
!
        subroutine setup_arr2(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl)
!!!!        use rpn_comm
        implicit none 
!
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!        fill array with markers ignoring k subscript
!
        integer minx,maxx,miny,maxy,nk
        integer nil,njl,z(2,minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k
        z = -1
        do k=1,nk
          do j=1,njl
          do i=1,nil
              z(1,i,j,k) = jtab(j)*100
              z(1,i,j,k) = z(1,i,j,k) + itab(i)*100000
              z(2,i,j,k) = k
          enddo
          enddo
        enddo
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(1,i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        subroutine setup_arr(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        fill array with markers
!
        integer minx,maxx,miny,maxy,nk
        integer nil,njl,z(minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k,k0
        z = -1
        do k=1,nk
          do j=1,njl
            do i=1,nil
              z(i,j,k) = k
              z(i,j,k) = z(i,j,k) + jtab(j)*100
              z(i,j,k) = z(i,j,k) + itab(i)*100000
            enddo
          enddo
        enddo
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(i,j,1),i=minx,maxx)
          enddo
        endif
100     format(10z10)
        return
        end
!
        SUBROUTINE sss(nx,ny)
        implicit none 
        integer zouf,nx,ny
!
!        "callback routine" used to get initial topology
!        information
!
        integer :: nodex, nodey
        common /pernode/ nodex, nodey
        common /the_problem/ nptsx,nptsy,nptsz,ihalox,ihaloy
        integer nptsx,nptsy,nptsz,ihalox,ihaloy
        integer deltai,deltaj
        open(5,file='TEST_data_001',form='FORMATTED')
        print *,'PEs =',nx*ny
        read(5,*)nx,ny,nptsx,nptsy,nptsz,ihalox,ihaloy,deltai,deltaj,nodex,nodey
        print *, ' problem size is ',nptsx,' by ',nptsy,' by ',nptsz
        print *, ' halo size is ',ihalox,' by ',ihaloy
        print *, ' topology = ',nx,' by ',ny
        print *, ' PE block topology = ',deltai,' by ',deltaj
        print *, 'Node tiles = ',nodex,' by ',nodey
        call RPN_COMM_set_petopo(deltai,deltaj)
        return
        end
!
        SUBROUTINE affichage(g,minx,maxx,miny,maxy,nk)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
        integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
        real g(minx:maxx,miny:maxy,nk)
!
        integer i,j,k,m
           if(Pelocal.eq.0) then
              write(*,*) 'matrice',Pelocal
 100       format(30F12.0)
           do j=maxy,miny,-1
              write(*,100) (g(i,j,1),i=minx,maxx)
           enddo
           endif
!
!
        return
        end
!
	integer function test_grids(mygrid,Pelocal)
	implicit none
	integer mygrid, Pelocal

	integer RPN_COMM_comm, RPN_COMM_colors
	external RPN_COMM_comm, RPN_COMM_colors
	integer mygridcomm, mymultigridcomm,myworldcomm,myallgridcomm,peercomm
	integer irank, isize, ierr, npeers

        test_grids = -1
        mymultigridcomm=RPN_COMM_comm('MULTIGRID')
        myallgridcomm=RPN_COMM_comm('ALLGRIDS')
        mygridcomm=RPN_COMM_comm('GRID')
        peercomm=RPN_COMM_comm('GRIDPEERS')
        myworldcomm=RPN_COMM_comm('WORLD')
!
        call RPN_COMM_BARRIER('GRID',ierr)
        call RPN_COMM_rank( 'GRID', irank ,ierr )
        call RPN_COMM_size( 'GRID', isize ,ierr )
        if(irank.eq.0)     print *,'BARRIER on GRID',mygridcomm,isize

        call RPN_COMM_BARRIER('MULTIGRID',ierr)
        call RPN_COMM_rank( 'MULTIGRID', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on MULTIGRID',mymultigridcomm

        call RPN_COMM_BARRIER('ALLGRIDS',ierr)
        call RPN_COMM_rank( 'ALLGRIDS', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on ALLGRIDS',myallgridcomm

        call RPN_COMM_BARRIER('WORLD',ierr)
        call RPN_COMM_rank( 'WORLD', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on WORLD',myworldcomm

        call RPN_COMM_BARRIER('GRIDPEERS',ierr)
        call RPN_COMM_rank( 'GRIDPEERS', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on GRIDPEERS',peercomm

        call RPN_COMM_size( 'GRIDPEERS', npeers ,ierr )
        if(Pelocal.eq.0 )then
          print *,'mygrid=',mygrid,'myworld',myworldcomm,'mymultigridcomm=',mymultigridcomm,'mygridcomm=',mygridcomm
          print *,'ID=',RPN_COMM_colors('WORLD'),'/',RPN_COMM_colors('MULTIGRID'),'/',RPN_COMM_colors('GRID')
          print *,'peer to peer comm=',peercomm,' npeers=',npeers
        endif
        call RPN_COMM_BARRIER('WORLD',ierr)
        test_grids = 0
	return
	end
        subroutine timing_report(n,times,label)
        implicit none
        integer :: n
        character (len=*) :: label
        real, dimension(3,n) :: times
        real tmin(3), tmax(3)
        real *8 ttot(3), ttot2(3)
        integer i,j
        character (len=10) :: labels(3)

        labels(1)=' PRE COPY '
        labels(2)=' MPI CALL '
        labels(3)=' POST COPY'

!            print *,label,' N=',n
!            print 111,'T1 ',nint(times(1,:)*1000000)
!            print 111,'T2 ',nint(times(2,:)*1000000 - times(1,:)*1000000)
!            print 111,'T3 ',nint(times(3,:)*1000000 - times(2,:)*1000000)
!111         format(A,15I6)
!        return
        times = times * 1000000
        ttot = 0.0
        tmin = 1000000000.0
        tmax = 0.0
        do j = 1 , n
         times(3,j) = times(3,j) - times(2,j)
         times(2,j) = times(2,j) - times(1,j)
         do i = 1 , 3
          tmin(i)=min(tmin(i),times(i,j))
          tmax(i)=max(tmax(i),times(i,j))
          ttot(i)=ttot(i)+times(i,j)
         enddo
        enddo
        ttot = ttot / n
        do I = 1 , 3
          print 222,label//labels(i)//' (min,max,avg) microsec',nint(tmin(i)),nint(tmax(i)),nint(ttot(i))
        enddo
222     format(A,3I8)
        return
        end
