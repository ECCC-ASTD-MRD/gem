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
! Fichier/File   : chm_load_store_tracers.ftn90
! Creation       : A. Akingunola  ,  GEM-MACH, Feb 2020
! Description    : Load/store the advected chemical tracers from/to the
!                  dynamic bus
!
! Arguments:
!
!  busdyn      --> Adress of DYN bus
!  chem_tr     --> Chemical species' array.
!  flag        --> Flag to load or store tracers from/to the dynamic bus
!
!==============================================================================
!
!!if_on
subroutine chm_load_store_tracers(busdyn, chem_tr, flag)
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
!!if_off

   use chm_species_info_mod, only: species_master
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: flag
   real(kind=4),    dimension(:), pointer, contiguous :: busdyn
   real(kind=4),    intent(inout) :: chem_tr(chm_ni, chm_nk + 1, nb_dyn_tracers)
!!if_off
!
! Local Variables
!
   integer(kind=4) :: isp, ii, kk, indx, busid
!
   if (flag == 0) then
    ! Copy tracers concentrations from the dynamic bus
      do isp = 1, nb_dyn_tracers
         if (species_master(isp) % dyn_offset > 0) then
            busid = species_master(isp) % dyn_offset
            indx = 0
            do kk = 1, chm_nk + 1
               do ii = 1, chm_ni
                  chem_tr(ii, kk, isp) = busdyn(busid + indx)
                  indx = indx + 1
               end do
            end do
         else
            chem_tr(:, :, isp) = 0.0
         end if
      end do
   else
      do isp = 1, nb_dyn_tracers
         if (species_master(isp) % dyn_offset > 0) then
            busid = species_master(isp) % dyn_offset
            indx = 0
            do kk = 1, chm_nk + 1
               do ii = 1, chm_ni
                  busdyn(busid + indx) = chem_tr(ii, kk, isp)
                  indx = indx + 1
               end do
            end do
         end if
      end do
   end if

   return
end subroutine chm_load_store_tracers
