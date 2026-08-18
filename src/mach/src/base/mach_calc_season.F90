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
! Projet / Project : GEM-MACH
! Fichier / File   : mach_calc_season.ftn90
! Creation         : Leiming Zhang - Oct 2000
! Description      : Define seasonal category according to latitude and month
!                    based on monthly averaged temperature.
!
! Extra Info       : Adapted for GEM-MACH by A.Kallaur and P.A. Beaulieu, Jan 2008
!                    Adding compatibility with Southern hemisphere by A. Robichaud Nov 2008.

!                    Seasonal categories
!                    1  midsummer with lush vegetation
!                    2  autumn with cropland before harvest
!                    3  lter autumn after frost, no snow
!                    4  winter, snow on ground and subfreezing
!                    5  transitional spring with partially green short annuals
!
! Arguments:  IN
!                metvar2d -> local storage (chemistry lib only, see chm_exe) for met
!
!             OUT
!                iseasn   -> Assigned season descriptors
!================================================================================================
!
!!if_on
subroutine mach_calc_season(metvar2d, iseasn)
   use chm_ptopo_grid_mod,   only: chm_ni
   use chm_metvar_mod,       only: SIZE_MV2D
!!if_off
   use chm_metvar_mod,       only: MV2D_DLAT
   use chm_consphychm_mod,   only: pi
   use chm_datime_mod,       only: ijul_day !         Julian day (1-365)
   implicit none
!!if_on
   integer(kind=4), intent(out) :: iseasn  (chm_ni)
   real(kind=4),    intent (in) :: metvar2d(chm_ni, SIZE_MV2D)
!!if_off
!
! Local variables
!
   integer(kind=4) :: i, jjul
   real(kind=4)    :: rad_to_deg
   real(kind=4)    :: glt(chm_ni)
!
!  Making a local conversion of latitudes obtained from the permanent
!  bus from radians to decimal degrees.
   rad_to_deg   = 180.0 / pi
   glt(1 : chm_ni) = metvar2d(1:chm_ni, MV2D_DLAT) * rad_to_deg

!  Set seasons over the whole domain
   do  i = 1, chm_ni

!  Modification for validity in SH (shift of 180 days vs NH)
      if (glt(i) < 0.0) then
         glt(i) = abs(glt(i))
         jjul = ijul_day + 180
         if (jjul >= 366) then
            jjul = jjul - 365
         end if
      else 
         jjul = ijul_day
      end if

! December January February
      if (jjul <= 59 .or. jjul >= 335) then

         if (glt(i) <= 30.0) then
            iseasn(i) = 1
         else if (glt(i) > 30.0 .and. glt(i) <= 35.0) then
            iseasn(i) = 5
         else if (glt(i) > 35.0) then
            iseasn(i) = 4
         end if

! March
      else if (jjul >= 60 .and. jjul <= 90)  then

         if (glt(i) <= 35.0) then
            iseasn(i) = 1
         else if (glt(i) > 35.0 .and. glt(i) <= 40.0) then
            iseasn(i) = 5
         else if (glt(i) > 40.0) then
            iseasn(i) = 4
         end if

! April
      else if (jjul >= 91 .and. jjul <= 120) then

         if (glt(i) <= 35.0) then
            iseasn(i) = 1
         else if (glt(i) > 35.0 .and. glt(i) <= 50.0) then
            iseasn(i) = 5
         else if (glt(i) > 50.0) then
            iseasn(i) = 4
         end if

! May
      else if (jjul >= 121 .and. jjul <= 151) then

         if (glt(i) <= 45.0) then
            iseasn(i) = 1
         else if (glt(i) > 45.0 .and. glt(i) <= 60.0) then
            iseasn(i) = 5
         else if (glt(i) > 60.0) then
            iseasn(i) = 4
         end if

! June
      else if (jjul >= 152 .and. jjul <= 181) then

         if (glt(i) <= 55.0) then
            iseasn(i) = 1
         else if (glt(i) > 55.0 .and. glt(i) <= 70.0) then
            iseasn(i) = 5
         else if (glt(i) > 70.0) then
            iseasn(i) = 4
         end if

! July
      else if (jjul >= 182 .and. jjul <= 212) then

         if (glt(i) <= 70.0) then
            iseasn(i) = 1
         else if (glt(i) > 70.0 .and. glt(i) <= 80.0) then
            iseasn(i) = 5
         else if (glt(i) > 80.0) then
            iseasn(i) = 4
         end if

! August
      else if (jjul >= 213 .and. jjul <= 243) then

         if (glt(i) <= 55.0) then
            iseasn(i) = 1
         else if (glt(i) > 55.0 .and. glt(i) <= 80.0) then
            iseasn(i) = 2
         else if (glt(i) > 80.0) then
            iseasn(i) = 4
         end if

! September
      else if (jjul >= 244 .and. jjul <= 273) then

         if (glt(i) <= 35.0) then
            iseasn(i) = 1
         else if (glt(i) > 35.0 .and. glt(i) <= 65.0) then
            iseasn(i) = 2
         else if (glt(i) > 65.0 .and. glt(i) <= 80.0) then
            iseasn(i) = 3
         else if (glt(i) > 80.0) then
            iseasn(i) = 4
         end if

! October
      else if (jjul >= 274 .and. jjul <= 304) then

         if (glt(i) <= 35.0) then
            iseasn(i) = 1
         else if (glt(i) > 35.0 .and. glt(i) <= 50.0) then
            iseasn(i) = 2
         else if (glt(i) > 50.0 .and. glt(i) <= 65.0) then
            iseasn(i) = 3
         else if (glt(i) > 65.0) then
            iseasn(i) = 4
         end if

! November
      else if (jjul >= 305 .and. jjul <= 334) then

         if (glt(i) <= 35.0) then
            iseasn(i) = 1
         else if (glt(i) > 35.0 .and. glt(i) <= 40.0) then
            iseasn(i) = 2
         else if (glt(i) > 40.0 .and. glt(i) <= 55.0) then
            iseasn(i) = 3
         else if (glt(i) > 55.0) then
            iseasn(i) = 4
         end if

      end if

   end do

   return

end subroutine mach_calc_season
