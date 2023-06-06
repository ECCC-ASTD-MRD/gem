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
! Fichier/File   : mach_landuse.ftn90
! Creation       : Paul Makar, Deji Akingunola, and Junhua Zhang - February 2020
! Description    : Set the month's LAI from the input 12-month satellite LAI, and
!                  calculate a seasonal mask from the annual LAI variation for possible
!                  use in later biogenic emissions
!
! Arguments:
!            IN
!              landuse  --> 15 landuse for dry deposition
!
!==============================================================================
!
!!if_on
 subroutine mach_lai_adjust(busper, landuse, trnch)
   use chm_ptopo_grid_mod,   only: chm_ni
   use mach_drydep_mod,      only: lucprm
!!if_off
   use mach_drydep_mod,      only: laindex_sat, nmth
   use chm_utils_mod,        only: chm_error_l, chm_msg_debug
   use chm_datime_mod,       only: imonth
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_LAI, sp_STSE
   use phymem,               only: phymeta, phyvar, phymem_find
!
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: trnch
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: landuse(chm_ni, lucprm)
!!if_off
!
!  Local variables:
   real(kind=4), dimension(nmth, chm_ni) :: lai_mnt
   real(kind=4)                          :: laimin, laimax
   integer(kind=4)                       :: ii, im, jj, istat
!
   type(phyvar) :: myvar(1)
   type(phymeta), pointer :: vmeta
!

   call msg_toall(chm_msg_debug, 'mach_lai_adjust [BEGIN]')

   istat = phymem_find(myvar, 'LAI_ENT', F_npath='V', F_bpath='E', &
                          F_quiet=.false., F_shortmatch=.false.)
   if (istat < 0) then
      write(*, *) 'Error in retrieving LAI from the entry bus'
      chm_error_l = .true.
      return
   end if
   vmeta => myvar(1)%meta
   jj = 0
   do im = 1, nmth
      do ii = 1, chm_ni
         jj = jj + 1
         lai_mnt(im, ii) = max(vmeta%bptr(jj, trnch), 0.0)
      end do
   end do
!
!  Determine seasonality mask for later use in biogenic emissions following
!  Junhua Zhang's suggestion
!
!  Resulting field sat_season_mask values:
!  Default, and for non-lai containing grid cells: 1.0; "summer"
!  lai-containing grid cells: 0.0 to 1.0, depending on the amount of foliage.
!    Evergreen forests will always have the "summer" value of 1.0.  This
!    follows the winter_std_p and summer_std_p designations in the BEIS algorithms,
!    which are 1.0 for summer, 0.5 for winter for deciduous species, with
!    evergreens and season-independent foliage set to 1.0.
!
!  For all other vegetation, the summer or winter value will depend on the foliage amount
!  with the maximum foliage being assumed to occur in summer, the minimum in winter.  The
!  linear interpolation between the two replaces the previous "season" mask in the biogenic
!  emissions processing setup
!
   do ii = 1, chm_ni
      laimin = 1000.
      laimax = 0.
      do im = 1, nmth
         do jj = 1, lucprm
            if (landuse(ii, jj) > 0.0001 .and. laindex_sat(jj, im) > 0.0) then
               laimax = max(laimax, lai_mnt(im, ii))
               laimin = min(laimin, lai_mnt(im, ii))
            end if
         end do
      end do
!
!  Determine current LAI relative to grid-cell maximum and minimum
         ! default value = summer (or constant LAI)
      busper(sm(sp_STSE) % per_offset + ii - 1) = 1.0
      do jj = 1, lucprm
         if (laimax > 0.0 .and. (laimax - laimin) > 0.0001) then
            busper(sm(sp_STSE) % per_offset + ii - 1) = &
                             (lai_mnt(imonth, ii) - laimin) / (laimax - laimin)
         end if
      end do
!
!  save the current month's LAI field in the permanent bus for subsequent calculations
      busper(sm(sp_LAI) % per_offset + ii - 1) = lai_mnt(imonth, ii)
!
   end do
   call msg_toall(chm_msg_debug, 'mach_lai_adjust [END]')
   return
 end subroutine mach_lai_adjust
