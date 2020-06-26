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
!**s/r set_poic - preparation of projection matrix in the east-west,
!                   C grid model

      subroutine set_poic ( F_eval_8, F_xp0_8, F_xp2_8, &
                            F_npts, NSTOR )
      use HORgrid_options
      use tdpack
      use glb_ld
      use glb_pil
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_npts, NSTOR
      real(kind=REAL64), dimension(NSTOR), intent(out) ::  F_eval_8
      real(kind=REAL64), dimension(NSTOR,3) :: F_xp0_8, F_xp2_8
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_eval_8     O    - eigenvalue vector
! F_xp0_8      I    - alongx projector diagonals
! F_xp2_8      I    - alongx second derivatives diagonals
! F_npts       I    - number of points to operate on
!


      integer :: i, j
      real(kind=REAL64) :: a_8(F_npts,F_npts), b_8(F_npts,F_npts), r_8(F_npts)
      real(kind=REAL64) :: pdfaz
      real(kind=REAL64), parameter :: one = 1.d0, two = 2.d0
!
! --------------------------------------------------------------------
!
!     put input arguments in real(kind=REAL64) arrays

      F_eval_8 = 0.d0
      a_8 = 0.d0
      b_8 = 0.d0

      do i = 1, F_npts-1
         a_8(i,i+1) = F_xp2_8(i+Lam_pil_w,3)
         a_8(i,i  ) = F_xp2_8(i+Lam_pil_w,2)
         a_8(i+1,i) = a_8(i,i+1)
         b_8(i,i+1) = F_xp0_8(i+Lam_pil_w,3)
         b_8(i,i  ) = F_xp0_8(i+Lam_pil_w,2)
         b_8(i+1,i) = b_8(i,i+1)
      end do

      a_8(F_npts,F_npts)= F_xp2_8(NSTOR-Lam_pil_e,2)
      b_8(F_npts,F_npts)= F_xp0_8(NSTOR-Lam_pil_e,2)

      call eigen ( r_8, a_8, b_8, F_npts, F_npts, 3*F_npts-1 )

      if (.not.Grd_yinyang_L) then
         pdfaz= ONE/sqrt((G_xg_8(G_ni-glb_pil_e)-G_xg_8(glb_pil_w)))
         r_8(1)          = 0.0d0
         a_8(1:F_npts,1) = pdfaz
      end if

!     put real(kind=REAL64) results in output arguments arrays
!     NOTE: for non-LAM case, the loops are F_npts=NSTOR

      do j = 1,F_npts
         F_eval_8(j+Lam_pil_w) = r_8(j)
      end do
!
! --------------------------------------------------------------------
!
      return
      end
