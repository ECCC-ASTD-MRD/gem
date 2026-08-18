!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_incld_steady_mod.ftn90
! Creation       : S. Menard, S. Gravel, GEM-MACH, June 2008
! Description    : Diagnostically calculate H+, OH- AND FEMN
!
! Extra info     : ADOM  VERSION: ADOMIIB(CLEAN)  LEVEL: 06/09/89  STEADY  ENSR(AES)
!
!                : Winter 2020 (P. Makar, W. Gong, and A. Akingunola)
!                  Remove the FE/MN update. For aqueous phase chemistry Fe and Mn are
!                  there to evaluate the metal catalysed SIV-> SVI oxidation reaction
!                  rates only
!
! Arguments 
!
!            INOUT
!              X         --> H+ and OH-
!
!=============================================================================
!
 module mach_incld_steady_mod
   use mach_cam_utils_mod, only: maxnsaq
   private
   public :: mach_incld_steady
!
   interface mach_incld_steady
      module procedure mach_incld_steady1
      module procedure mach_incld_steady2
   end interface mach_incld_steady

 contains

   subroutine mach_incld_steady1(x, nptsnz, nq)
      implicit none
      integer(kind=4), intent   (in) :: nq
      integer(kind=4), intent   (in) :: nptsnz
      real(kind=4),    intent(inout) :: x(nptsnz, maxnsaq, 2)
!
! Local variables
!
      integer(kind=4) :: i
      real(kind=4)    :: xvec(maxnsaq, 2)

      do i = 1, nptsnz
         xvec = x(i, :, :)
         call mach_incld_steady2(xvec, nq)
         x(i, :, :) = xvec
      end do

      return
   end subroutine mach_incld_steady1

   subroutine mach_incld_steady2(x, nq)
      implicit none
      integer(kind=4), intent   (in) :: nq
      real(kind=4),    intent(inout) :: x(maxnsaq, 2)
!
! Local variables
!
      integer(kind=4) :: i
      real(kind=4)    :: xhi, xoh, xh

      xhi = x(1, nq) + 2.0 * x(4, nq) + x(5, nq) - x(6, nq) + &
            x(8, nq) - x(7, nq)
      xoh = 0.0
      do i = 1, 2
         xh = max(3.0e-7, xhi + xoh)
         xoh = 1.e-14 / xh
      end do
!  h+
      x(9, nq) = xh
!  oh-
      x(10, nq) = xoh

      return
   end subroutine mach_incld_steady2

 end module mach_incld_steady_mod

