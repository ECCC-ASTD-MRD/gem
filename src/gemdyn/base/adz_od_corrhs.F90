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

      subroutine adz_od_corrhs (F_stk,F_cptr,F_nptr,F_geom,F_export)
      use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
      use adz_mem
      use ptopo
      implicit none

      integer, intent(in) :: F_nptr
      type(C_PTR), intent(in) :: F_geom
      type(C_PTR), dimension (F_nptr), intent(in) :: F_cptr
      type(ADZ_SLOD), intent(INOUT) :: F_export
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk

      include "tricublin_f90.inc"
      integer :: nr, dim, err_dest
      real, dimension(Adz_MAX_MPI_OS_SIZE*3) :: request_from
      real, dimension(Adz_MAX_MPI_OS_SIZE*F_nptr) :: send
!
!     ---------------------------------------------------------------
!
! get requested positions from others and post return adresses

      dim= Adz_MAX_MPI_OS_SIZE*3
      call adz_traj_from (request_from,nr,F_export,F_nptr,dim)

! compute the requested interpolated values and put the results in
! Adz_cor which is exposed to all PEs through Adz_wincor window

      call tricublin_zyx1_m_n ( send, F_cptr, request_from,&
                                F_geom ,nr, F_nptr )

      err_dest=0
      if (F_nptr*nr>size(Adz_cor)) then
         err_dest= -1
      else
         call reshape_2_win (Adz_cor,send,nr, F_nptr)
      endif
      call gem_error(err_dest,'adz_od_corrhs',&
        'NOT ENNOUGH MEMORY - Adz_MAX_MPI_OS_SIZE for Adz_cor')
      
! fetch and plug interpolated requests from neighbors

      call fetch_n_plug (F_stk,F_export,F_nptr)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_od_corrhs
      
      subroutine reshape_2_win (F_cor,send,nr, F_nptr)
      implicit none
      integer nr, F_nptr
      real, dimension(*) :: send
      real, dimension(F_nptr,nr) :: F_cor

      integer i,j
      do i=1,nr
         do j=1,F_nptr
            F_cor(j,i)= send(i+(j-1)*nr)
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine reshape_2_win
