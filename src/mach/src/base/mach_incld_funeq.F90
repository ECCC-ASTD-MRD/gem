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
! Fichier/File   : mach_incld_funeq.ftn90
! Creation       : S. Menard,  W. Gong, S. Gravel, V. Bouchet, B. Pabla, GEM-MACH,  June 2008.
!
! Description    : Compute Total S(IV), Total NO3, Total AMMONIA, Total CARBON DIOXIDE
!                  and Net charge due to SO4=
!
! Extra info     :
!
! Arguments:
!           IN
!              c      -->    ALocal Array of Aq. and Gas species
!              b      -->    Aqueous variable coefficients
!              nptsnz -->    Total number of grids to integrate
!
!           OUT
!              Y      -->    Component Total in
!                            Y(ivec,1) : Total S(IV)
!                            Y(ivec,2) : Total NO3
!                            Y(ivec,3) : Total AMMONIA
!                            Y(ivec,4) : Total CARBON DIOXIDE
!                            Y(ivec,5) : Net charge due to SO4=
!=============================================================================
!
!!if_on
subroutine mach_incld_funeq (y, c, b, nptsnz)
!!if_off

   implicit none
!!if_on
   integer(kind=4), intent (in)  :: nptsnz
   real(kind=4),    intent (in)  :: b(nptsnz, 5)
   real(kind=4),    intent (in)  :: c(nptsnz, 12)
   real(kind=4),    intent(out)  :: y(nptsnz, 5)
!!if_off
!
! Local variables
!
   integer(kind=4) :: ivec

   do ivec = 1, nptsnz
!  total s(iv)
      y(ivec, 1) = c(ivec, 1) + b(ivec, 2) * c(ivec, 6)

!  total no3
      y(ivec, 2) = c(ivec, 2) + b(ivec, 2) * c(ivec, 7)

!  total ammonia
      y(ivec, 3) = c(ivec, 3) + b(ivec, 2) * c(ivec, 8)

!  total carbon dioxide
      y(ivec, 4) = c(ivec, 4) + b(ivec, 2) * c(ivec, 9)

!  net charge due to so4= and cat1
      y(ivec, 5) = 2.0 * c(ivec, 11) - c(ivec, 12)

   end do

   return
end
