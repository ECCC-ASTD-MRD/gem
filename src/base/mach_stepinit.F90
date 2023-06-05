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
! Fichier/File   : mach_stepinit.ftn90
! Creation       : Paul Makar, Deji Akingunola, Feb 2016
! Description    : Setup for chemistry (GEM-MACH) step
!
! Extra info     :
!
! Arguments:
!           IN
!              step              -> Timestep number
!              ni_can            -> Number of surface grids with canopy shading
!
!           IN/OUT
!              busper            -> Permanent bus
!
!=============================================================================
!
!!if_on
subroutine mach_stepinit(busper, step, trnch, ni_can)
!!if_off
   use chm_utils_mod,        only: chm_error_l, LONG_VARNAME, chm_timestep, &
                                   chm_msg_debug
   use chm_nml_mod,          only: chm_canopy_shading_l
   use chm_species_info_mod, only: sm, chem_ent_vars
   use chm_species_idx_mod,  only: sp_HC, sp_FRT, sp_CLU, sp_LAI, var_step
   use chm_datime_mod,       only: imonth
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use mach_pkg_misc_mod,    only: ent_vars_num
   use phymem,               only: phymeta, phyvar, phymem_find
   use chm_phyvar_mod,       only: dxdy

   implicit none
!!if_on
   integer(kind=4), intent   (in) :: step, trnch
   integer(kind=4), intent  (out) :: ni_can
   real(kind=4),    dimension(:), pointer, contiguous :: busper
!!if_off
!
!  Declaration of local variables
!
   integer(kind=4)  :: ii, istat, jj, busid
   real(kind=4), dimension(:), pointer :: pop, frt, lai
   real(kind=4)     :: flandnonforest !  all non-forested land use types
   real(kind=4)     :: beam_penetration, pdens
! minimum value of the leaf area index for canopy shading to happen
   real(kind=4), parameter :: laimin   = 0.1
! minimum height of plant canopy for shading to happen
   real(kind=4), parameter :: hcmin    = 0.5
! minimum height of plant canopy for beam penetration to matter
   real(kind=4), parameter :: hcmin2   = 18.0
! To make the canopy code more consistent across different grid cell
! sizes, maximum population density in mks (people/m2) is utilized.
! This maximum value is derived from the 10km-resolution value for the
! maximum population of 10000 people within the grid cell for which the
! canopy code is applied.
! (i.e. 10000 people / (1E4 x 1E4) m^2)
! maximum population density (canopy will be applied if less than this value)
   real(kind=4), parameter :: pdensmax   = 1.0E-4
! minimum value of nonforested land fraction within grid square
! (if >= this value, no canopy is applied)
   real(kind=4), parameter :: noforestcutoff = 0.50
!  If beam penetration is > this value, no canopy
   real(kind=4), parameter :: beam_min = 0.45
!
   character(len=LONG_VARNAME) :: vname
   type(phyvar) :: myvar(1)
   type(phymeta), pointer :: cmeta

   call msg_toall(chm_msg_debug, 'mach_stepinit [BEGIN]')
!
!==============================================================================
!  Initialize the KPP solver internal timestep
!==============================================================================
   if (step < 1 .and. var_step > 0) then
      busid = sm(var_step) % per_offset
      do jj = 1, chm_nk
         do ii = 1, chm_ni
            busper(busid) = chm_timestep
            busid = busid + 1
         end do
      end do
   end if
!
   ni_can = 0
   if (.not. chm_canopy_shading_l) return
!
   if (step < 1) then
!
!==============================================================================
!  Initialize the KPP solver internal timestep for the canopy columns
!==============================================================================
      if (var_step > 0) then
         busid = sm(var_step) % dd_offset
         do jj = 1, 3
            do ii = 1, chm_ni
               busper(busid) = chm_timestep
               busid = busid + 1
            end do
         end do
      end if
!
!==============================================================================
!  Determine canopy columns for the forest canopy model
!==============================================================================
      nullify(lai, pop, frt)
      do ii = 1, ent_vars_num
         vname = chem_ent_vars(ii) % ent_name
         istat = phymem_find(myvar, trim(vname), F_npath='V', F_bpath='E', &
                             F_quiet=.false., F_shortmatch=.false.)
         if (istat < 0) then
            write(*, *) 'Error in retrieving ', trim(vname), ' from the entry bus'
            chm_error_l = .true.
            return
         end if
         cmeta => myvar(1)%meta
         if (trim(vname) == 'LAI_ENT') then
            jj = (imonth - 1) * chm_ni
            lai(1:chm_ni) => cmeta%bptr(1+jj:chm_ni+jj, trnch)
         end if
         if (trim(vname) == 'POPU_ENT') &
            pop(1:chm_ni) => cmeta%bptr(1:chm_ni, trnch)
         if (trim(vname) == 'FRT_BELD3') &
            frt(1:chm_ni) => cmeta%bptr(1:chm_ni, trnch)
      end do
!
!    If not using the monthly LAI input;
      if (.not. associated(lai)) then
         lai(1:chm_ni) => busper(sm(sp_LAI) % per_offset:)
      end if
!
!  Determine the number of columns in the current slice containing a forest canopy
      do ii = 1, chm_ni
!
!  use BELD3 data for flandnonforest:
         flandnonforest = 1.0 - frt(ii)
         flandnonforest = max(0.0, flandnonforest)
!
!  Calculate probability of direct downward beam of light reaching surface:
!
         beam_penetration = exp(-0.5 * lai(ii) * &
                                busper(sm(sp_CLU) % per_offset + ii - 1))
!
! Determine which model columns will be treated with a forest canopy, and
! store the information in sm(sp_FRT) as -1.0 or 1.0
! Population aspect of the criteria is based on the population density.
! Population field, used in calculation of population density, is limited
! to minimum of 1.E-03 to prevent possibility of underflow in the
! population density decision in the subsequent IF statement.
         pdens = max(pop(ii),1.E-03) / busper(dxdy + ii - 1)
         if (lai(ii) < laimin .or.      &
             busper(sm(sp_HC) % per_offset + ii - 1) < hcmin .or. &
             pdens > pdensmax .or.      &
             flandnonforest >= noforestcutoff            .or.    &
             ((beam_penetration > beam_min) .and.                &
             (busper(sm(sp_HC) % per_offset + ii - 1) < hcmin2))) then

            busper(sm(sp_FRT) % per_offset + ii - 1) = -1.0
         else
            busper(sm(sp_FRT) % per_offset + ii - 1) = 1.0
            ni_can = ni_can + 1
         end if
      end do
!
   else
!
!===============================================================================
!  Determine canopy columns for the forest canopy model
!===============================================================================
!
      do ii = 1, chm_ni
         if (busper(sm(sp_FRT) % per_offset + ii - 1) > 0.0) then
            ni_can = ni_can + 1
         end if
      end do
!
   end if
!
   call msg_toall(chm_msg_debug, 'mach_stepinit [END]')
   return
end subroutine mach_stepinit
