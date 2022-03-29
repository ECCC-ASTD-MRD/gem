      subroutine gem_xch_halo ( f, lminx,lmaxx,lminy,lmaxy,nk)
      use iso_c_binding
      use gem_options
      use inp_mod
      use ptopo
      use glb_ld
      implicit none

      integer, intent (IN) :: lminx,lmaxx,lminy,lmaxy,nk
      real, intent (INOUT) :: f(lminx:lmaxx,lminy:lmaxy,nk)

      include 'mpif.h'
      include "rpn_comm.inc"
      integer i, j, k, len, dsp, err
      integer ORIGIN_DATATYPE, TARGET_DATATYPE, TARGET_RANK
      integer ORIGIN_COUNT, TARGET_COUNT
      integer(C_INTPTR_T) :: TARGET_DISP
      real, dimension((l_maxx-l_minx+1)*(l_maxy-l_miny+1)*nk/4,8) :: buf
!
!---------------------------------------------------------------------
!
      dsp= (max(G_lnimax,G_lnjmax)+2*max(G_halox,G_haloy))*max(G_halox,G_haloy)*nk

      ORIGIN_DATATYPE= RPN_COMM_datyp('MPI_REAL')
      TARGET_DATATYPE= ORIGIN_DATATYPE
      call MPI_Win_fence(0, Inp_window, err)
      if (.not. l_south) then
         TARGET_DISP= 4*dsp
         TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol,Ptopo_myrow-1)
         len = 0
         do k = 1, nk
            do j = 1, G_haloy
               do i = 1-G_halox*west, l_ni+G_halox*east
                  len      = len + 1
                  buf(len,1) = f(i,j,k)
               end do
            end do
         end do
         ORIGIN_COUNT = len
         TARGET_COUNT = ORIGIN_COUNT
         call MPI_Put( buf(1,1), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Inp_window, err)
         if (.not. l_west) then
            TARGET_DISP= 5*dsp
            TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol-1,Ptopo_myrow-1)
            len = 0
            do k = 1, nk
               do j = 1, G_haloy
                  do i = 1, G_halox
                     len      = len + 1
                     buf(len,2) = f(i,j,k)
                  end do
               end do
            end do
            ORIGIN_COUNT = len
            TARGET_COUNT = ORIGIN_COUNT
            call MPI_Put( buf(1,2), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                          TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                          TARGET_DATATYPE, Inp_window, err)
         endif
         if (.not. l_east) then
            TARGET_DISP= 3*dsp
            TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol+1,Ptopo_myrow-1)
            len = 0
            do k = 1, nk
               do j = 1, G_haloy
                  do i = l_ni-G_halox+1, l_ni
                     len      = len + 1
                     buf(len,3) = f(i,j,k)
                  end do
               end do
            end do
            ORIGIN_COUNT = len
            TARGET_COUNT = ORIGIN_COUNT
            call MPI_Put( buf(1,3), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                          TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                          TARGET_DATATYPE, Inp_window, err)
         endif
      endif
      if (.not. l_north) then
         TARGET_DISP= dsp
         TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol,Ptopo_myrow+1)
         len = 0
         do k = 1, nk
            do j = l_nj-G_haloy+1, l_nj
               do i = 1-G_halox*west, l_ni+G_halox*east
                  len      = len + 1
                  buf(len,4) = f(i,j,k)
               end do
            end do
         end do
         ORIGIN_COUNT = len
         TARGET_COUNT = ORIGIN_COUNT
         call MPI_Put( buf(1,4), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Inp_window, err)
         if (.not. l_west) then
            TARGET_DISP= 2*dsp
            TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol-1,Ptopo_myrow+1)
            len = 0
            do k = 1, nk
               do j = l_nj-G_haloy+1, l_nj
                  do i = 1, G_halox
                     len      = len + 1
                     buf(len,5) = f(i,j,k)
                  end do
               end do
            end do
            ORIGIN_COUNT = len
            TARGET_COUNT = ORIGIN_COUNT
            call MPI_Put( buf(1,5), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                          TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                          TARGET_DATATYPE, Inp_window, err)
         endif
         if (.not. l_east) then
            TARGET_DISP= 0
            TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol+1,Ptopo_myrow+1)
            len = 0
            do k = 1, nk
               do j = l_nj-G_haloy+1, l_nj
                  do i = l_ni-G_halox+1, l_ni
                     len      = len + 1
                     buf(len,6) = f(i,j,k)
                  end do
               end do
            end do
            ORIGIN_COUNT = len
            TARGET_COUNT = ORIGIN_COUNT
            call MPI_Put( buf(1,6), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                          TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                          TARGET_DATATYPE, Inp_window, err)
         endif
      endif
      
      if (.not. l_west) then
         TARGET_DISP= 7*dsp
         TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol-1,Ptopo_myrow)
         len = 0
         do k = 1, nk
            do j = 1-G_haloy*south, l_nj+G_haloy*north
               do i = 1, G_halox
                  len      = len + 1
                  buf(len,7) = f(i,j,k)
               end do
            end do
         end do
         ORIGIN_COUNT = len
         TARGET_COUNT = ORIGIN_COUNT
         call MPI_Put( buf(1,7), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Inp_window, err)
      endif
      if (.not. l_east) then
         TARGET_DISP= 6*dsp
         TARGET_RANK= Ptopo_colrow(0,Ptopo_mycol+1,Ptopo_myrow)
         len = 0
         do k = 1, nk
            do j = 1-G_haloy*south, l_nj+G_haloy*north
               do i = l_ni-G_halox+1, l_ni
                  len      = len + 1
                  buf(len,8) = f(i,j,k)
               end do
            end do
         end do
         ORIGIN_COUNT = len
         TARGET_COUNT = ORIGIN_COUNT
         call MPI_Put( buf(1,8), ORIGIN_COUNT, ORIGIN_DATATYPE   ,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Inp_window, err)
      endif
      
      call MPI_Win_fence(0, Inp_window, err)
      
      if (.not. l_south) then
         len = dsp
         do k = 1, nk
            do j = 1-G_haloy, 0
               do i = 1-G_halox*west, l_ni+G_halox*east
                  len      = len + 1
                  f(i,j,k) = Inp_recv(len)
               end do
            end do
         end do
         if (.not. l_west) then
            len = 0
            do k = 1, nk
               do j = 1-G_haloy, 0
                  do i = 1-G_halox, 0
                     len      = len + 1
                     f(i,j,k) = Inp_recv(len)
                  end do
               end do
            end do
         endif
         if (.not. l_east) then
            len = 2*dsp
            do k = 1, nk
               do j = 1-G_haloy, 0
                  do i = l_ni+1, l_ni+G_halox
                     len      = len + 1
                     f(i,j,k) = Inp_recv(len)
                  end do
               end do
            end do
         endif
      endif
      if (.not. l_north) then
         len = 4*dsp
         do k = 1, nk
            do j = l_nj+1, l_nj+G_haloy
               do i = 1-G_halox*west, l_ni+G_halox*east
                  len      = len + 1
                  f(i,j,k) = Inp_recv(len)
               end do
            end do
         end do
         if (.not. l_west) then
            len = 3*dsp
            do k = 1, nk
               do j = l_nj+1, l_nj+G_haloy
                  do i = 1-G_halox, 0
                     len      = len + 1
                     f(i,j,k) = Inp_recv(len)
                  end do
               end do
            end do
         endif
         if (.not. l_east) then
            len = 5*dsp
            do k = 1, nk
               do j = l_nj+1, l_nj+G_haloy
                  do i = l_ni+1, l_ni+G_halox
                     len      = len + 1
                     f(i,j,k) = Inp_recv(len)
                  end do
               end do
            end do
         endif
      endif
      if (.not. l_west) then
         len = 6*dsp
         do k = 1, nk
            do j = 1-G_haloy*south, l_nj+G_haloy*north
               do i = 1-G_halox, 0
                  len      = len + 1
                  f(i,j,k) = Inp_recv(len)
               end do
            end do
         end do
      endif
      if (.not. l_east) then
         len = 7*dsp
         do k = 1, nk
            do j = 1-G_haloy*south, l_nj+G_haloy*north
               do i = l_ni+1, l_ni+G_halox
                  len      = len + 1
                  f(i,j,k) = Inp_recv(len)
               end do
            end do
         end do
      endif
!     
!---------------------------------------------------------------------
!
      return
      end subroutine gem_xch_halo
