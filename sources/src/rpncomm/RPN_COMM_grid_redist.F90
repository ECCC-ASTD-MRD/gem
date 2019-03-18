!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
      integer function RPN_COMM_grid_redist_test(nparams,params)
      use rpn_comm
      implicit none
      integer :: nparams
      integer, dimension(100) :: params
        
      integer :: RPN_COMM_grid_redist
      external :: RPN_COMM_grid_redist
      integer :: RPN_COMM_grid_redist_n
      external :: RPN_COMM_grid_redist_n
      integer, pointer, dimension(:,:,:) :: localarray
      integer, pointer, dimension(:,:,:) :: globalarray
      integer*8, pointer, dimension(:,:,:) :: localarray2
      integer*8, pointer, dimension(:,:,:) :: globalarray2
      integer :: gni, gnj
      integer :: lminx,lmaxx,lminy,lmaxy
      integer, dimension(pe_nx) :: countx,offsetx
      integer, dimension(pe_ny) :: county,offsety
      integer :: i, ii, j, jj, k, ierr, ltok, status, status2
      integer :: ix0, ixn, jy0, jyn
      integer :: i0, in, j0, jn
      integer, parameter :: lni=3
      integer, parameter :: lnj=5
      integer :: nox
      integer :: noy
      integer :: maxz
      integer, pointer, dimension(:) :: outpe_x
      integer, pointer, dimension(:) :: outpe_y
      integer, pointer, dimension(:) :: zlist, zlist2
      integer :: nk, nz, k0
      real *8 :: t1,t2,t3,t4,t5
      integer, external :: RPN_COMM_limit
!
!      gni = params(1)
      gni = pe_nx*lni
!      gnj = params(2)
      gnj = pe_ny*lnj
!      nk = params(3)
      nk = 7 !nox*noy*3-1
!      nox = params(4)
      nox = 2
!      noy = params(5)
      noy = 2
!
      ix0 = 4 ! lni+(lni+1)/2-1
      ixn = 6 ! gni-2*lni+(lni+1)/2+1
      jy0 = 6 ! (lnj+1)/2 + lnj-1
      jyn = 10 ! gnj-2*lnj+(lnj+1)/2+1
!
      allocate(outpe_x(nox),outpe_y(noy))
      allocate(zlist(nk),zlist2(nk))
      maxz = (nk+nox+noy-1)/(nox*noy)
      if(pe_me==0) write(rpn_u,*)'grid redistribution test',&
     &    pe_tot_grid,pe_nx,pe_ny

      outpe_x(nox) = pe_nx-1
      outpe_y(noy) = pe_ny-1
      outpe_x(nox-1) = pe_nx-2
      outpe_y(noy-1) = pe_ny-2
      outpe_x(1) = 0
      outpe_y(1) = 0
      ierr = RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      ierr = RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
      allocate(localarray(lminx-1:lmaxx+2,lminy-2:lmaxy+1,nk))
      allocate(localarray2(lminx-1:lmaxx+2,lminy-2:lmaxy+1,nk))
      localarray = 0
      do k = 1,nk
      do j = lminy,lmaxy
      do i = lminx,lmaxx
        localarray(i,j,k) = k + 10*j + 1000*i
        localarray2(i,j,k) = k + 10*j + 1000*i
      enddo
      enddo
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(pe_me==0)then
        write(rpn_u,100)'outpe_x =',outpe_x
        write(rpn_u,100)'outpe_y =',outpe_y
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(pe_me==(pe_nx*pe_ny-1)) then
        write(rpn_u,100)'local array allocated',pe_mex,pe_mey,&
     &                lminx,lmaxx,lminy,lmaxy,nk
        do j = lmaxy,lminy,-1
!          write(rpn_u,100)' ',localarray(lminx:lmaxx,j,max(1,nk-1))
        enddo
100     format(A,10I6)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      nullify(globalarray)
      i0=1   ! max(lminx,ix0)
      if(ix0>lminx) i0 = i0 + (ix0-lminx)
      in=lni ! min(lmaxx,ixn)
      if(lmaxx>ixn) in = in - (lmaxx-ixn)
      j0=1   ! max(lminy,jy0)
      if(jy0>lminy) j0 = j0 + (jy0-lminy)
      jn=lnj ! min(lmaxy,jyn)
      if(lmaxy>jyn) jn = jn - (lmaxy-jyn)
      do j = 1,noy
      do i = 1,nox
        if(outpe_x(i)==pe_mex .and. outpe_y(j)==pe_mey)then
!          write(rpn_u,*)'global array allocated',ix0, ixn, jy0, jyn
          allocate(globalarray(ix0:ixn,jy0:jyn,maxz))
          allocate(globalarray2(ix0:ixn,jy0:jyn,maxz))
        endif
      enddo
      enddo
      if(.not. associated(globalarray)) allocate(globalarray(1,1,maxz))
      if(.not. associated(globalarray2))allocate(globalarray2(1,1,maxz))
      globalarray = 0
      do i=1,pe_nx*pe_ny
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(.false. .and. pe_me==i-1)&
     &  write(rpn_u,101)'global array dims',&
     &                pe_mex,pe_mey,shape(globalarray)
101     format(A,10I5)
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      zlist = -999999
      zlist2 = -999999
      nz = maxz
      ltok = 1
      t1=mpi_wtime()
      status = RPN_COMM_grid_redist_n(&
     & localarray,1-1,lni+2,i0,in,1-2,lnj+1,j0,jn,nk,&
     & globalarray,ix0,ixn,jy0,jyn,nz,zlist,&
     & gni,gnj,outpe_x,nox,outpe_y,noy,1)
      t2=mpi_wtime()
      RPN_COMM_grid_redist_test = status
      if(pe_me==0)write(rpn_u,*)'status=',status,' time=',t2-t1
      if(status==-1)goto 777
      
      goto 1111
      status2 = RPN_COMM_grid_redist_n(&
     & localarray2,1-1,lni+2,i0,in,1-2,lnj+1,j0,jn,nk,&
     & globalarray2,ix0,ixn,jy0,jyn,nz,zlist2,&
     & gni,gnj,outpe_x,nox,outpe_y,noy,2)
      RPN_COMM_grid_redist_test = status2
      if(pe_me==0)write(rpn_u,*)'status2=',status2
      if(status2==-1)goto 777
      
1111  continue
      do i=1,pe_nx*pe_ny
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(pe_me==i-1) then
          write(rpn_u,*)'pe=',pe_mex,pe_mey,' zlist=',zlist
          write(rpn_u,*)'pe=',pe_mex,pe_mey,' zlist2=',zlist2
          call flush(rpn_u)
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      do k0 = 1 , maxz
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(zlist(k0)>0) then
          k=zlist(k0)
          ierr=vfy_array(globalarray,ix0,ixn,jy0,jyn,2,k0,k)
          write(rpn_u,102)'pe=',pe_mex,pe_mey,&
     &         ',  level=',k,',  errors=',ierr
102       format(A,2I3,A,I3,a,I3)
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo
      goto 2222
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      do k0 = 1 , maxz
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(zlist2(k0)>0) then
          k=zlist2(k0)
          ierr=vfy_array2(globalarray2,ix0,ixn,jy0,jyn,2,k0,k)
          write(rpn_u,102)'pe2=',pe_mex,pe_mey,&
     &         ',  level=',k,',  errors=',ierr
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo
2222  continue
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      do k = 1,nk
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(any(k==zlist)) then
          write(rpn_u,103)'global array =',&
     &          pe_mex,pe_mey,ix0,ixn,jy0,jyn
103       format(A,20(1X,I5.5))
          do k0=1,2   ! size(zlist)
            if(zlist(k0)==k)then
              do jj = jyn,jy0,-1
                write(rpn_u,103)' ',globalarray(ix0:ixn,jj,k0)
!                write(rpn_u,103)'+',globalarray2(ix0:ixn,jj,k0)
              enddo
            endif
          enddo
          call flush(rpn_u)
        endif
      enddo
777   if(associated(globalarray)) deallocate(globalarray)
      deallocate(localarray)
      return
!
      contains
!
      integer function vfy_array(zin,mini,maxi,minj,maxj,nk,k,ref)
      implicit none
      integer, intent(IN) :: mini,maxi,minj,maxj,nk,k,ref
      integer, dimension(mini:maxi,minj:maxj,nk),&
     &         intent(IN) :: zin
!
      integer :: i,j,nerr
!
      nerr = 0
      do j=minj,maxj
      do i=mini,maxi
        if(zin(i,j,k)/=ref+10*j+1000*i)then
          nerr=nerr+1
!          print *,'error',i,j,k,ref,zin(i,j,k),ref+10*j+1000*i
        endif
      enddo
      enddo
      vfy_array = nerr
!
      return
      end function vfy_array
 !
      integer function vfy_array2(zin,mini,maxi,minj,maxj,nk,k,ref)
      implicit none
      integer, intent(IN) :: mini,maxi,minj,maxj,nk,k,ref
      integer*8, dimension(mini:maxi,minj:maxj,nk),&
     &         intent(IN) :: zin
!
      integer :: i,j,nerr
!
      nerr = 0
      do j=minj,maxj
      do i=mini,maxi
        if(zin(i,j,k)/=ref+10*j+1000*i)then
          nerr=nerr+1
!          print *,'error',i,j,k,ref,zin(i,j,k),ref+10*j+1000*i
        endif
      enddo
      enddo
      vfy_array2 = nerr
!
      return
      end function vfy_array2
     
      end function RPN_COMM_grid_redist_test
!=======================================================================
!=======================================================================
      integer function RPN_COMM_grid_redist(&
     &           zin,mini,maxi,i0,in,minj,maxj,j0,jn,nk,&
     &           zout,ix0,ixn,jy0,jyn,nz,zlist,&
     &           gni,gnj,outpe_x,noutpe_x,outpe_y,noutpe_y,ltok)
      use rpn_comm
      implicit none
      integer, intent(IN) :: mini,maxi,i0,in,minj,maxj,j0,jn,nk,ltok
      integer, intent(IN) :: ix0,ixn,jy0,jyn,nz
      integer, intent(IN) :: noutpe_x,noutpe_y,gni,gnj
      integer, dimension(mini:maxi,minj:maxj,nk),&   ! (ltok,mini:maxi,minj:maxj,nk)
     &         intent(IN) :: zin  ! zin(:,io:in,j0:jn,:) is useful
      integer, dimension(ix0:ixn,jy0:jyn,nz), intent(OUT) :: zout  ! (ltok,ix0:ixn,jy0:jyn,nz)
      integer, dimension(nz), intent(OUT) :: zlist  ! list of 2D fields returned to this processor
      integer, dimension(noutpe_x), intent(IN) :: outpe_x    ! list of columns where there are PEs doing IO
      integer, dimension(noutpe_y), intent(IN) :: outpe_y    ! list of rows where there are PEs doing IO
!
      integer, dimension(pe_nx) :: istart_g, iend_g, countx, offsetx
      integer, dimension(pe_ny) :: jstart_g, jend_g, county, offsety
      integer, dimension(pe_nx) :: gdispl_x, gcounts_x, gsize_x
      integer, dimension(pe_ny) :: gdispl_y, gcounts_y
      integer, dimension(noutpe_x) :: level_x
      integer, dimension(nz*noutpe_y) :: level_table
      integer :: istart, iend, jstart, jend
      integer, pointer, dimension(:)     :: dest, dest2, ptr1d
      integer, pointer, dimension(:,:)   :: source2
      integer, pointer, dimension(:,:,:) :: source
      integer :: lminx, lmaxx, lminy, lmaxy
      integer :: i, j, k, l, kbot, ktop, n2d, kout
      integer :: lev0, narrays, ierr, temp, my_out_col
      logical :: needed_for_pass2, size_error, no_holes
      integer, external :: RPN_COMM_limit
!
      RPN_COMM_grid_redist = -1  ! return -1 if an error occurred
      zlist = -1                 ! empty list
      if(nk>nz*noutpe_x*noutpe_y)return

      nullify(ptr1d)        ! nullify pointers to prevent mishaps
      allocate(dest(1))           ! possibly temporary allocation, to be increased later on
      allocate(source(1,1,1))     ! if necessary (data to collect on PE)
      allocate(dest2(1))
      allocate( source2(1,1))     ! array for gatherv pass along y
      dest = 0 ; source = 0 ; dest2 = 0 ; source2 = 0
!
      narrays = nk
      needed_for_pass2 = .false.
      size_error = .false.
      if(mod(gni,ltok)/=0) return  ! gni must be a multiple of tlok
!!!      if(mod(gni,pe_nx)/=0) return ! not safe yet if tiles not same size along X
!
!     may have to add some adjustments to account for ltok factor
!     ltok > 1 may be unsafe if all tile sizes along X are not equal
!     i.e. when mod(gni,pe_nx) is not zero
!     ltok is hidden inside mini,maxi,i0,in,ix0,ixn,gni
!     user must call RPN_COMM_grid_redist_n if ltok>1
!     as it will "fudge" mini,maxi,i0,in,ix0,ixn,gni appropriately
!     mod(gni,ltok) must be zero
!
      ierr = RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)  ! X decomposition
      ierr = RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)  ! Y decomposition
!
      level_x = 0
      level_table = -1
      temp = 0
      do k = 1 , narrays  ! distribute narrays 2D fields over noutpe_x columns
        temp = temp + 1
        if(temp > noutpe_x) temp = 1
        level_x(temp) = level_x(temp) + 1
      enddo
      n2d=-1
      temp=0
      my_out_col=-1
      do i = 1,noutpe_x
!        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(pe_mex==outpe_x(i)) then   ! my column, 
          n2d = level_x(i)      ! number of 2D fields for this PE column
          my_out_col=i
          needed_for_pass2 = .true.
          do j=1,level_x(i)
            level_table(j) = temp+j   ! list of levels processed by the PE column
          enddo
        endif
        temp=temp+level_x(i)
      enddo
!
      do j = 1 , pe_ny   ! start and end of valid data on PE(any,j) along y in global space
        jstart_g(j) = max(jy0,1+offsety(j))
        jend_g(j)   = min(jyn,offsety(j) + county(j))
      enddo
      jstart = jstart_g(pe_mey+1)  ! same as above but for my row
      jend   = jend_g(pe_mey+1)
      do i = 1 , pe_nx   ! start and end of valid data on PE(i,any) along x in global space
        istart_g(i) = max(ix0,1+offsetx(i))
        iend_g(i)   = min(ixn,offsetx(i) + countx(i))
      enddo
      istart = istart_g(pe_mex+1)  ! same as above but for my column
      iend   = iend_g(pe_mex+1)
      if((in-i0).ne.(iend-istart) .or. (jn-j0).ne.(jend-jstart)) then  ! size consistency problem ?
        if(istart<=iend) then
          size_error = .true.  ! add error message 
          write(rpn_u,100)'error on pe',pe_mex,pe_mey,' =',&
     &          i0,in,istart,iend,j0,jn,jstart,jend
100       format(A,2I2,A,10I8)
        endif
      endif
!     check if a size error occurred somewhere in pe_grid, if so return -1
      call MPI_allreduce(MPI_IN_PLACE,size_error,1,MPI_LOGICAL,MPI_LOR,&
     &                   pe_grid,ierr)
      if(size_error) goto 8888
!
      if(jstart <= jend) then  ! there is something to do during pass 1 for this PE row
        do i = 1 , pe_nx                                  ! horizontal size of 2D fragments
          gsize_x(i)  = max(0,iend_g(i)-istart_g(i)+1) *&  ! gsize_x may be zero for some PEs
     &                  max(0,jend-jstart+1)
        enddo
!
        kbot = 1
        do l = 1 , noutpe_x
          if(level_x(l)==0) cycle  ! no output to do on this column
          ktop = kbot + level_x(l) - 1
!         calculate counts and displacements for gatherv along x
!          gcounts_x = 1 + gsize_x * level_x(l)  ! 1 + "horizontal size" * "levels to gather"
          gcounts_x = gsize_x * level_x(l)  !  "horizontal size" * "levels to gather"
          gdispl_x(1) = 0
          do i = 2 , pe_nx
            gdispl_x(i) = gdispl_x(i-1) + gcounts_x(i-1)
          enddo
!
          if(istart<=iend) then ! there is valid data on this PE
            if(associated(source)) deallocate(source)
            allocate(source(istart:iend,jstart:jend,kbot:ktop))     ! source buffer
            source=77777
            source(istart:iend,jstart:jend,kbot:ktop) = &
     &         zin(i0:in,j0:jn,kbot:ktop)   ! extract subarray from input array
          endif
          if(pe_mex == outpe_x(l)) then  ! PEs doing gathering on this column
            if(associated(dest))deallocate(dest)
            allocate( dest( gdispl_x(pe_nx)+gcounts_x(pe_nx)+2 ) )   ! destination buffer 
            dest=99999
            if(associated(source2))deallocate(source2)
            allocate( source2(ix0:ixn,jstart:jend) )  ! reallocate with proper dimensions
            source2=88888
          endif
          call MPI_Gatherv(source,gcounts_x(pe_mex+1),MPI_INTEGER,&
     &                     dest,gcounts_x,gdispl_x,MPI_INTEGER,&
     &                     outpe_x(l),pe_myrow,ierr)
          kbot = kbot + level_x(l)
        enddo
        deallocate(source)
        nullify(source)
      endif   ! (jstart <= jend)
      if(associated(source)) deallocate(source) ! source no longer needed at this point
!
!      data gathering along X is done, we have the whole X axis in processor now
!
!      write(rpn_u,*)'end pass 1 :',pe_mex,pe_mey,needed_for_pass2,n2d
      if(needed_for_pass2 .and. n2d>0) then  ! nothing to do if n2d=0
!       get X counts and displacements for THIS column
!        gcounts_x = 1 + gsize_x * level_x(my_out_col)
        gcounts_x = gsize_x * level_x(my_out_col)
        gdispl_x(1) = 0
        do i = 2 , pe_nx
          gdispl_x(i) = gdispl_x(i-1) + gcounts_x(i-1)
        enddo
        lev0 = 0
        no_holes=.true.  ! will be true if one or more PE rows contribute no data
                         ! this rigamarole is needed because some MPI implementations
                         ! of gatherv fail with zero data counts. the alternative is 
                         ! to create a new communicator for part of the column.
        do j = 1 , pe_ny ! gcounts_y contains the size of a "full x" "partial y" slab
          gcounts_y(j) = (ixn-ix0+1) * max(0,jend_g(j)-jstart_g(j)+1)
!          no_holes = no_holes .and. gcounts_y(j)>0
          if(lev0==0 .and. gcounts_y(j)>0) lev0 = j ! find first non empty slab
!          if(gcounts_y(j)==0) gcounts_y(j)=1  ! minimum count set to 1 or gatherv might fail
        enddo
        gdispl_y(1) = 0
        do j = 2 , pe_ny
          gdispl_y(j) = gdispl_y(j-1) + gcounts_y(j-1)
        enddo
        do l = 1 , noutpe_y
          if(pe_mey==outpe_y(l)) then  ! this PE is a data gatherer
            zlist = -1
            if(.not. no_holes) then        ! will need an extra copy
              if(associated(dest2)) deallocate(dest2)
              allocate( dest2( gdispl_y(pe_ny)+gcounts_y(pe_ny)+1 ) )
            endif
          endif
        enddo
        do k = 1 , n2d  ! loop over number of 2D fields to distribute along y
!         reassemble data for this "level" from dest(localx,localy,localk,npex)
!                                          into source2(globalx,localy,localk)
          do i = 1 , pe_nx   ! pe_nx pieces to reassemble
            if(jstart<=jend) then
              ptr1d => dest(gdispl_x(i)+1:gdispl_x(i)+gcounts_x(i))
              call place(&
     &           ptr1d,istart_g(i),iend_g(i),jstart,jend,n2d,&
     &           source2,ix0,ixn,k)
              endif
          enddo
          kout = 1 + (k-1)/noutpe_y
          l = 1 + mod(k-1,noutpe_y)
          if(pe_mey == outpe_y(l)) zlist(kout) = level_table(k)
          if(no_holes) then  ! all PEs contribute, can gather directly into zout
            call MPI_Gatherv(&
     &         source2,gcounts_y(pe_mey+1),MPI_INTEGER,&
     &         zout(ix0,jy0,kout),gcounts_y,gdispl_y,MPI_INTEGER,&
     &         outpe_y(l),pe_mycol,ierr)
          else  ! not all PEs contribute, nedd an intermediate array
            call MPI_Gatherv(&         !  gatherv, then copy into zout
     &         source2,gcounts_y(pe_mey+1),MPI_INTEGER,&
     &         dest2,gcounts_y,gdispl_y,MPI_INTEGER,&
     &         outpe_y(l),pe_mycol,ierr)
            if(pe_mey == outpe_y(l)) then  ! we are on root PE of gather
              temp = 1 + gdispl_y(lev0) ! point to start of useful data
              do j = jy0 , jyn          ! copy reassembled data into zout
              do i = ix0 , ixn
                zout(i,j,kout) = dest2(temp)
                temp = temp + 1
              enddo  ! j
              enddo  ! i
            endif  ! pe_mey == outpe_y(l)
          endif  ! no_holes
        enddo  ! k = 1 , n2d
      endif  ! needed_for_pass2 .and. n2d>0
7777  continue
      RPN_COMM_grid_redist = 1  ! number of 2D reassembled fields returned
8888  if(associated(dest2))   deallocate(dest2)
      if(associated(dest))    deallocate(dest)
      if(associated(source2)) deallocate(source2)
!
!      write(rpn_u,*)'end of pass 1+2, pe=',pe_me
      return
!
      contains
      subroutine place(src,i0,in,j0,jn,nk,dest,ix0,ixn,k)
      implicit none
      integer :: i0,in,j0,jn,nk,ix0,ixn,k
      integer, dimension(i0:in,j0:jn,nk), intent(IN) :: src
      integer, dimension(ix0:ixn,j0:jn), intent(OUT) :: dest
      dest(i0:in,j0:jn) = src(i0:in,j0:jn,k)
      end subroutine place
!
      end function RPN_COMM_grid_redist
!=======================================================================
!=======================================================================
      integer function RPN_COMM_grid_redist_n(&
     &           zin,mini,maxi,i0,in,minj,maxj,j0,jn,nk,&
     &           zout,ix0,ixn,jy0,jyn,nz,zlist,&
     &           gni,gnj,outpe_x,noutpe_x,outpe_y,noutpe_y,&
     &           ltok)
      use rpn_comm
      implicit none
      integer, intent(IN) :: mini,maxi,i0,in,minj,maxj,j0,jn,nk,ltok
      integer, intent(IN) :: ix0,ixn,jy0,jyn,nz
      integer, intent(IN) :: noutpe_x,noutpe_y,gni,gnj
      integer, dimension(ltok,mini:maxi,minj:maxj,nk), &
     &         intent(IN) :: zin  ! zin(:,io:in,j0:jn,:) is useful
      integer, dimension(ltok,ix0:ixn,jy0:jyn,nz), intent(OUT) :: zout
      integer, dimension(nz), intent(OUT) :: zlist  ! list of 2D fields returned to this processor
      integer, dimension(noutpe_x), intent(IN) :: outpe_x    ! list of columns where there are PEs doing IO
      integer, dimension(noutpe_y), intent(IN) :: outpe_y    ! list of rows where there are PEs doing IO
!
      external :: RPN_COMM_grid_redist
      integer :: RPN_COMM_grid_redist
!
      RPN_COMM_grid_redist_n = -1  ! return -1 if an error occurred
      zlist = -1                 ! empty list
!      if(ltok .ne. 1 ) then
!        if(pe_me==0) write(rpn_u,*)'ERROR: ltok >1 not supported yet'
!        return    ! for now, ltok MUST be = 1 (feature not implemented yet)
!      endif
      if(nk>nz*noutpe_x*noutpe_y)return
!
!     "fudge" mini,maxi,i0,in,ix0,ixn,gni to account for ltok before calling 
!     RPN_COMM_grid_redist that will perform the actual work
!     ltok is passed to RPN_COMM_grid_redist tht may need it to compute data
!     distribution correctly (where all tiles do not have the same size along X)
!
      RPN_COMM_grid_redist_n=RPN_COMM_grid_redist(zin,&
     &           1+(mini-1)*ltok,1+(maxi-1)*ltok+ltok-1,&
     &           1+(i0-1)*ltok,1+(in-1)*ltok+ltok-1,&
     &           minj,maxj,j0,jn,nk,&
     &           zout,1+(ix0-1)*ltok,1+(ixn-1)*ltok+ltok-1,&
     &           jy0,jyn,nz,zlist,&
     &           gni*ltok,gnj,&
     &           outpe_x,noutpe_x,outpe_y,noutpe_y,ltok)
!
      return
      end function RPN_COMM_grid_redist_n
