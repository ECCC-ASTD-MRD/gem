      subroutine glbdist_os ( f2bc,f2rc,lminx,lmaxx,lminy,lmaxy,&
                              nka, glbmaxx,glbmaxy, zlist,nz, mult,add )
      use iso_c_binding
      use gem_options
      use HORgrid_options
      use inp_mod
      use ptopo
      use glb_ld
      implicit none
      
      integer, intent (IN) :: lminx,lmaxx,lminy,lmaxy,nka,&
                              glbmaxx,glbmaxy,nz
      integer, intent(IN) :: zlist(nz)
      real, intent (IN)  :: f2bc(lminx:glbmaxx,lminy:glbmaxy,*)
      real, intent (OUT) :: f2rc(lminx:lmaxx  ,lminy:lmaxy  ,nka)
      real(kind=REAL64), intent (IN) :: add, mult

      include 'mpif.h'
      include "rpn_comm.inc"
      logical mem_L
      integer dst_level, i, j, k, iproc, err, i0, in, j0, jn, len, cnt
      integer ORIGIN_DATATYPE, TARGET_DATATYPE, TARGET_RANK
      integer ORIGIN_COUNT, TARGET_COUNT
      integer(C_INTPTR_T) :: TARGET_DISP
      real, dimension(:,:), allocatable :: buf
!
!---------------------------------------------------------------------
!     
      ORIGIN_DATATYPE= RPN_COMM_datyp('MPI_REAL')
      TARGET_DATATYPE= ORIGIN_DATATYPE

      mem_L= any(zlist(1:nz) > -1)
      if (mem_L) allocate (buf((l_maxx-l_minx+1)*(l_maxy-l_miny+1),Ptopo_numproc))
      cnt= 0

      do dst_level= 1, nz
         call MPI_Win_fence(0, Inp_window, err)
         if (zlist(dst_level) > -1) then
            k= zlist(dst_level) ; cnt=cnt+1
            do iproc= 1, Ptopo_numproc
               i0= Ptopo_gindx(1,iproc)-G_halox
               in= Ptopo_gindx(2,iproc)+G_halox
               j0= Ptopo_gindx(3,iproc)-G_haloy
               jn= Ptopo_gindx(4,iproc)+G_haloy
               len = 0
               do j = j0, jn
                  do i = i0, in
                     len      = len + 1
                     buf(len,iproc) = f2bc(i,j,cnt)
                  end do
               end do
               ORIGIN_COUNT = len
               TARGET_COUNT = ORIGIN_COUNT
               TARGET_RANK  = iproc-1
               TARGET_DISP  = (k-1)*(in-i0+1)*(jn-j0+1)
               call MPI_Put( buf(1,iproc),ORIGIN_COUNT,ORIGIN_DATATYPE,&
                             TARGET_RANK,TARGET_DISP, TARGET_COUNT    ,&
                             TARGET_DATATYPE, Inp_window, err)
            end do
         endif
         call MPI_Win_fence(0, Inp_window, err)
      end do
      if (mem_L) deallocate(buf)

      cnt=0
      do k=1,nka
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               cnt= cnt+1
               f2rc(i,j,k)= Inp_recv(cnt) * mult + add
            end do
         end do
      end do
!     
!---------------------------------------------------------------------
!
      return
      end subroutine glbdist_os
