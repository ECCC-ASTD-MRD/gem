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
! Fichier/File   : mach_incld_concmp.ftn90
! Creation       : S. Menard, W. Gong, S. Gravel, GEM-MACH, June 2008.
!
! Description    : Compute gas and aqueous species concentration calculated from [H+]
!
! Extra info     :
!
! Arguments:
!           IN
!
!             nptsnz        -->   Total number of grids to integrate
!             T             -->   Component Total in
!                                     T(ivec,1) : Total S(IV)
!                                     T(ivec,2) : Total NO3
!                                     T(ivec,3) : Total AMMONIA
!                                     T(ivec,4) : Total CARBON DIOXIDE
!                                     T(ivec,5) : Net charge due to SO4
!                                     and CAT1
!             R            -->    Aqueous rate constant array
!             B            -->    Aqueous variable coefficients
!
!           OUT
!             C            -->   Gas and Aq. species concentration calculated from [H+]
!
!=============================================================================
!
!!if_on
subroutine mach_incld_concmp (t, C, r, b, nptsnz)
!!if_off

   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nptsnz
   real(kind=4),    intent   (in) :: t(nptsnz,5)
   real(kind=4),    intent(inout) :: c(nptsnz,12)
   real(kind=4),    intent   (in) :: r(nptsnz,25)
   real(kind=4),    intent   (in) :: b(nptsnz,5)
!!if_off
!
!local variables
!
   integer(kind=4) :: ivec

   do ivec = 1, nptsnz
      c(ivec, 1)  = t(ivec, 1) * r(ivec, 3) * c(ivec, 5) / (r(ivec, 2) + r(ivec, 3) * c(ivec, 5))
      c(ivec, 2)  = t(ivec, 2) * r(ivec, 9) * c(ivec, 5) / (r(ivec, 8) + r(ivec, 9) * c(ivec, 5))
      c(ivec, 4)  = t(ivec, 4) * r(ivec, 16) * c(ivec, 5) / (r(ivec, 15) + r(ivec, 16) * c(ivec, 5))
      c(ivec, 10) = r(ivec, 18) / (r(ivec, 17) * c(ivec, 5))
      c(ivec, 3)  = t(ivec, 3) * r(ivec, 13) * c(ivec, 10) / (r(ivec, 12) + r(ivec, 13) * c(ivec, 10))
      c(ivec, 6)  = (t(ivec, 1) - c(ivec, 1)) / b(ivec, 2)
      c(ivec, 7)  = (t(ivec, 2) - c(ivec, 2)) / b(ivec, 2)
      c(ivec, 8)  = (t(ivec, 3) - c(ivec, 3)) / b(ivec, 2)
      c(ivec, 9)  = r(ivec, 15) * c(ivec, 4) / (b(ivec, 2) * r(ivec, 16) * c(ivec, 5))
  end do

  return
end
