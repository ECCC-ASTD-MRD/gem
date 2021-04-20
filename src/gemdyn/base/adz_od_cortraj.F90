!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

      subroutine adz_od_cortraj ()
      use adz_mem
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      include "tricublin_f90.inc"
      integer :: i,j,k, cnt,dim,nr,jj,ijk,km,err,req
      integer, dimension(Ptopo_world_numproc) :: rank, nreq, dsp
      real, dimension  (Adz_MAX_MPI_OS_SIZE*3) :: request_from
      real, dimension(:,:,:), allocatable :: cor
!     
!     ---------------------------------------------------------------
!
      if (.not.ADZ_OD_L) return
      
      dim= Adz_MAX_MPI_OS_SIZE*3
      call adz_traj_from (request_from,nr,Adz_expq,3,dim)

! compute the requested interpolated values and put the results in
! Adz_cor which is exposed to all PEs through Adz_wincor window
      
      call tricublin_zyx3_n ( Adz_cor,Adz_uvw_d(1,1,1,1), &
                              request_from,Adz_cpntr_q,nr )

      call MPI_Win_fence(0, Adz_offs_win, err)
      call rpn_comm_barrier ("MULTIGRID", err)
      
      ORIGIN_DATATYPE = REAL_DATATYPE
      TARGET_DATATYPE = REAL_DATATYPE
      dim = l_ni * l_nj
      
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
      allocate (cor(3,Adz_MAX_MPI_OS_SIZE,nr))
      do req=1,nr
         ORIGIN_COUNT= nreq(req)*3
         TARGET_RANK = rank(req)
         TARGET_DISP = dsp (req)
         TARGET_COUNT= ORIGIN_COUNT
         call MPI_Get( cor(1,1,req), ORIGIN_COUNT, ORIGIN_DATATYPE,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Adz_wincor, err )
      end do
      call MPI_Win_fence(0, Adz_wincor, err)
      
      do req=1,nr
         do km=1, nreq(req)
            ijk= Adz_expq%dest(km,rank(req))
            k = ijk/dim - 1 + min(mod(ijk,dim),1)
            jj= ijk-k*dim
            j = jj/l_ni - 1 + min(mod(jj,l_ni),1)
            i = jj - j*l_ni ; j=j+1 ; k=k+1
            Adz_uvw_dep(1,i,j,k) = cor(1,km,req)
            Adz_uvw_dep(2,i,j,k) = cor(2,km,req)
            Adz_uvw_dep(3,i,j,k) = cor(3,km,req)
         end do
      end do
      
      deallocate (cor)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_od_cortraj
