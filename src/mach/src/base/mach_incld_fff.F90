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
! Fichier/File   : mach_incld_fff.ftn90
! Creation       : S. Menard,  W. Gong, S. Gravel, V. Bouchet, B. Pabla, GEM-MACH,  June 2008.
!
! Description    : Compute the charge balance function
!
! Extra info     :
!
! Arguments:
!           IN
!              T      -->    Initial component totals
!                            1-Total S(IV)
!                            2-Total NO3
!                            3-Total Ammonia
!                            4-Total Carbon Dioxide
!                            5-Net charge due to SO4= and CAT1
!              C      -->    Aq. proton concentration
!                            C= 0.1*[H+] or 10.0 * [H+]
!                            (units :  M/L)
!              R      -->    Aqueous rate constant array
!              B      -->    Aqueous variable coefficients
!              nleft  -->    Grids Integration start from "nleft" into memory location
!              nptsnz -->    Total number of grids to integrate
!           OUT
!              Fnew   -->    Charge balance function
!
!=============================================================================
!
!!if_on
subroutine mach_incld_fff (FNEW, t, c, r, b, nptsnz, nleft)
!!if_off

   implicit none
!!if_on
   
   integer(kind=4), intent (in) :: nptsnz
   real(kind=4),    intent(out) :: fnew(nptsnz)
   real(kind=4),    intent (in) :: t(nptsnz,5)
   real(kind=4),    intent (in) :: c(nptsnz)
   real(kind=4),    intent (in) :: r(nptsnz,25)
   real(kind=4),    intent (in) :: b(nptsnz,5)
   integer(kind=4), intent (in) :: nleft
!!if_off
!
!  Local variables
!
   integer(kind=4)  :: ivec
   real(kind=4)     :: a1, a2, a3, a4

   do ivec = nleft, nptsnz
      a1 = t(ivec, 1) * r(ivec, 2) / b(ivec, 2)
      a2 = t(ivec, 2) * r(ivec, 8) / b(ivec, 2)
      a4 = t(ivec, 4) * r(ivec, 15) / b(ivec, 2)
      a3 = t(ivec, 3) * r(ivec, 12) * r(ivec, 17) / b(ivec, 2)

      fnew(ivec) = c(ivec) - t(ivec, 5) - a1 / (r(ivec, 2) + r(ivec, 3) * c(ivec))         &
                   - a2 / (r(ivec, 8) + r(ivec, 9) * c(ivec))                              &
                   + a3 * c(ivec) / (r(ivec, 13) * r(ivec, 18) + r(ivec, 12) * r(ivec, 17) &
                   * c(ivec)) - a4 / (r(ivec, 15) + r(ivec, 16) * c(ivec))                 &
                   - r(ivec, 18) / (r(ivec, 17) * c(ivec))
   end do

return
end
