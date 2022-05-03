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

!**s/r  pre_jacobi2D -  Jacobi additive-Schwarz preconditioning
!
      subroutine pre_jacobi2D ( lhs, rhs, evec_local, ni, nj, ai, bi, ci, level )

      use glb_ld, only: l_nk ! Number of vertical levels
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
      ! csubich: as computed at initialization, ai, bi, and ci are three-dimensional
      ! arrays of size (ni,nj,nk) that are for unfortunate historical reasons allocated
      ! as giant one-dimensional arrays.  We need access to the appropriate vertical
      ! index to properly include the effects attribuable to the vertical eigenvalues.
      integer, intent(in) :: ni,nj, level
      real(kind=REAL64), dimension(ni,nj), intent(out) :: lhs
      real(kind=REAL64), dimension(ni,nj), intent(in) :: rhs
      real(kind=REAL64), dimension(ni,nj,l_nk), intent(in) :: ai, bi, ci
      real(kind=REAL64), dimension(ni,ni), intent(in) :: evec_local
!
!author
!       Abdessamad Qaddouri - December  2006
!
      integer :: i, j, jr
      real(kind=REAL64), dimension(ni,nj) :: fdg
!
!     ---------------------------------------------------------------
!
      call dgemm('T', 'N', ni, nj, ni, 1.d0, evec_local, ni, &
                    rhs(1,1), ni, 0.d0, fdg(1,1), ni)

      do j =2, Nj
         jr =  j - 1
         do i=1,Ni
            fdg(i,j) = fdg(i,j) - ai(i,j,level) * fdg(i,jr)
         end do
      end do

      j = Nj

      do i=1,Ni
         fdg(i,j) = fdg(i,j) / bi(i,j,level)
      end do

      do j = Nj-1, 1, -1
         jr =  j + 1
         do i=1, Ni
            fdg(i,j) = ( fdg(i,j)-ci(i,j,level)*fdg(i,jr) ) / bi(i,j,level)
         end do
      end do

      call dgemm('N','N', ni, nj, ni, 1.d0, evec_local, ni, &
                 fdg(1,1), ni, 0.d0, lhs(1,1), ni)
!
!     ---------------------------------------------------------------
!
      return
      end

