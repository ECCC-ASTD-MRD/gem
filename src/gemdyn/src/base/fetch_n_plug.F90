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

      subroutine fetch_n_plug (F_stk,F_export,F_nptr)
      use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
      use adz_mem
      use ptopo
      implicit none

      integer, intent(IN ) :: F_nptr
      type(ADZ_SLOD), intent(IN) :: F_export
      type(Adz_pntr_stack), dimension(F_nptr), intent(INOUT) :: F_stk

      integer :: i,j,k,n,jj,ijk,km,req,dim,cnt,nr,err
      integer, dimension(Ptopo_world_numproc) :: rank, nreq, dsp
      real, dimension(:,:,:), allocatable :: cor
!
!     ---------------------------------------------------------------
!
      call MPI_Win_fence(0, Adz_offs_win, err)
      call rpn_comm_barrier ("MULTIGRID", err)
      call MPI_Win_fence(0, Adz_wincor, err)
      
      ORIGIN_DATATYPE = REAL_DATATYPE
      TARGET_DATATYPE = REAL_DATATYPE
      dim= l_ni*l_nj
      
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
      allocate (cor(F_nptr,Adz_MAX_MPI_OS_SIZE,nr))
      do req=1,nr
         ORIGIN_COUNT= nreq(req) * F_nptr
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
            ijk= F_export%dest(km,rank(req))
            k = ijk/dim - 1 + min(mod(ijk,dim),1)
            jj= ijk-k*dim
            j = jj/l_ni - 1 + min(mod(jj,l_ni),1)
            i = jj - j*l_ni ; j=j+1 ; k=k+1
            do n=1, F_nptr
               F_stk(n)%dst(i,j,k) = cor(n,km,req)
            end do
         end do
      end do
      
      deallocate (cor)
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine fetch_n_plug 
