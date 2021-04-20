      subroutine adz_od_corflux (F_flk,F_cptr,F_nptr,F_geom,F_export,F_tr_num,F_flkn)

      use adz_mem
      use ptopo
      implicit none

      integer, intent(in) :: F_nptr,F_flkn
      type(C_PTR), intent(in) :: F_geom
      type(C_PTR), dimension (2*F_nptr), intent(in) :: F_cptr
      type(ADZ_SLOD), intent(INOUT) :: F_export
      type(flux), dimension(F_flkn) :: F_flk
      integer, intent(in) :: F_tr_num(F_flkn)

      include "tricublin_f90.inc"
      integer :: i,j,k,n, cnt,dim,nr,jj,km,err,ijk, req, err_dest
      integer, dimension(Ptopo_world_numproc) :: rank, nreq, dsp
      real, dimension(Adz_MAX_MPI_OS_SIZE*F_nptr*2) :: send
      real, dimension(Adz_MAX_MPI_OS_SIZE*3       ) :: request_from
      real, dimension(:,:,:,:), allocatable :: cor
!
!     ---------------------------------------------------------------
!
! get requested positions from others and post return adresses
      dim= Adz_MAX_MPI_OS_SIZE*3
      call adz_traj_from (request_from,nr,F_export,F_nptr*2,dim)

! compute the requested interpolated values and put the results in
! Adz_cor which is exposed to all PEs through Adz_wincor window

      call tricublin_zyx1_m_n (send, F_cptr, request_from, &
                               F_geom, nr, F_nptr*2)
      err_dest=0
      if (2*F_nptr*nr>size(Adz_cor)) then
         err_dest= -1
      else
         call flux_2_win (Adz_cor,send,nr,F_nptr)
      endif
      call gem_error(err_dest,'adz_od_corflux',&
           'NOT ENOUGH MEMORY - Adz_MAX_MPI_OS_SIZE for Adz_cor')

!     get requested positions from others

      ORIGIN_DATATYPE = REAL_DATATYPE
      TARGET_DATATYPE = REAL_DATATYPE
      dim= l_ni*l_nj

      call MPI_Win_fence(0, Adz_offs_win, err)
      call rpn_comm_barrier ("MULTIGRID", err)

      cnt= 0
      do req=1,2*Ptopo_world_numproc,2
         if (Adz_offs(req) > 0) then
            cnt=cnt+1
            rank(cnt)= req/2
            nreq(cnt)= Adz_offs(req)
            dsp (cnt)= Adz_offs(req+1)
         endif
         nr= cnt
      end do

! fetch and plug interpolated requests from neighbors
      allocate (cor(2,F_nptr,Adz_MAX_MPI_OS_SIZE,nr))
      do req=1,nr
         ORIGIN_COUNT= nreq(req) * 2*F_nptr
         TARGET_RANK = rank(req)
         TARGET_DISP = dsp (req)
         TARGET_COUNT= ORIGIN_COUNT
         call MPI_Get( cor(1,1,1,req), ORIGIN_COUNT, ORIGIN_DATATYPE,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Adz_wincor, err )
      end do
      call MPI_Win_fence(0, Adz_wincor, err)
      do req=1,nr
         do km=1, nreq(req)
            ijk= F_export%dest(km,rank(req))
            k = ijk/dim - 1 + min(mod(ijk,dim),1)
            jj= ijk-k*dim
            j = jj/l_ni - 1 + min(mod(jj,l_ni),1)
            i = jj - j*l_ni ; j=j+1 ; k=k+1
            do n=1, F_nptr
               F_flk(F_tr_num(n))%fo(i,j,k) = cor(1,n,km,req)
               F_flk(F_tr_num(n))%fi(i,j,k) = cor(2,n,km,req)
            end do
         end do
      enddo
      deallocate (cor)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_od_corflux

      subroutine flux_2_win (F_cor,send,nr,F_nptr)
      implicit none
      integer nr, F_nptr
      real, dimension(*) :: send
      real, dimension(2,F_nptr,nr) :: F_cor

      integer i,j
      do i=1,nr
         do j=1,F_nptr
            F_cor(1,j,i)= send(i+2*(j-1)*nr   )
            F_cor(2,j,i)= send(i+2*(j-1)*nr+nr)
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine flux_2_win
