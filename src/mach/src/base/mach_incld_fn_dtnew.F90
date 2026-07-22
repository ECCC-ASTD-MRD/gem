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
! Fichier/File   : mach_incld_fn_dtnew.ftn90
! Creation       : W. Gong, S. Gravel, GEM-MACH,  June 2008.
!
! Description    : Calculate new time step
!
! Extra info     :  ADOM    VERSION: ADOMIIA     LEVEL: 04/18/88     DTNEW        ERT
!
! Arguments:
!           IN
!
!              DELT  -->   Time step
!              K     -->   Iteration number
!
!=============================================================================
!
!!if_on
real function mach_incld_fn_dtnew (delt, k)
!!if_off
   implicit none
!!if_on
   integer(kind=4), intent(in) :: k
   real(kind=4),    intent(in) :: delt
!!if_off

   if (k > 2) then
      mach_incld_fn_dtnew = delt
   else if (delt >= 3.) then
      if (k == 1) mach_incld_fn_dtnew = 1.50 * delt
      if (k == 2) mach_incld_fn_dtnew = 1.25 * delt
   else if (delt >= .10) then
      if (k == 1) mach_incld_fn_dtnew = 3. * delt
      if (k == 2) mach_incld_fn_dtnew = 1.5 * delt
   else
      if (k == 1) mach_incld_fn_dtnew = 10. * delt
      if (k == 2) mach_incld_fn_dtnew = 3. * delt
   end if
!****
!    if (k > 2) then
!       mach_incld_fn_dtnew = delt
!    else if (k == 2) then
!       if (delt >= 3.) then
!          mach_incld_fn_dtnew = 1.25 * delt
!       else if (delt >= .10) then
!          mach_incld_fn_dtnew = 1.5 * delt
!       else
!          mach_incld_fn_dtnew = 3. * delt
!       end if
!    else if (k == 1) then
!       if (delt >= 3.) then
!          mach_incld_fn_dtnew = 1.50 * delt
!       else if (delt >= .10) then
!          mach_incld_fn_dtnew = 3.0 * delt
!       else
!          mach_incld_fn_dtnew = 10.0 * delt
!       end if
!    end if

   return
end function mach_incld_fn_dtnew
