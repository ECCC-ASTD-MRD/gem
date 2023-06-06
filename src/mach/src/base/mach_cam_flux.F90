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
! Fichier/File   : mach_cam_flux.ftn90
! Creation       : S. Gravel
!                  Jun. 2012
!
! Description    : Entry point for aerosol flux calculations
!
! Extra info     : Includes calculations of Sea-salt surface flux and/or modulation
!                  of anthropogenic dust emissions by meteorological conditions.
!
!                  The modulation factor "fmet" to be applied to fugitive dust
!                  emissions has value 0 or 1 depending on meteorological conditions.
!                  fmet scaling factor is applied on fugitive dust emissions in
!                  mach_diffusion.ftn90 (D. Akingunola, Fall 2018) - Mantis #2055.
!
! Arguments:  IN
!               metvar2d(ix, ...)     -> 2D met variables
!
!             OUT
!               busvol                -> Volatile bus
!
!             IN/OUT
!               busper                -> Permanent bus
!
!============================================================================
!
!!if_on
subroutine mach_cam_flux(busper, busvol, metvar2d, fland)
   use chm_metvar_mod,       only: SIZE_MV2D
   use chm_ptopo_grid_mod,   only: chm_ni
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_metvar_mod,       only: MV2D_TDIAG, MV2D_WSDIAG, MV2D_DXDY, &
                                   MV2D_WSOIL, MV2D_SNOF
   use chm_utils_mod,        only: global_debug, chm_lun_out
   use mach_cam_headers_mod, only: mach_cam_sfss
   use chm_nml_mod,          only: chm_seaflux_s, chm_met_modulation_s
   use mach_cam_utils_mod,   only: isize, iae_SS
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_AERO, sp_FMET
   implicit none
!!if_on
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: fland   (chm_ni, lucprm)
!!if_off
!
!  Local variables
!
   integer(kind=4)         :: i, k, nsp_SS
   real(kind=4), parameter :: kg2g = 1.0e3
   real(kind=4)            :: surfwd(chm_ni)
   real(kind=4)            :: tdiag (chm_ni)
   real(kind=4)            :: rsfrow(chm_ni, isize)
   real(kind=4)            :: snof
   logical(kind=4)         :: local_dbg

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
!
   if (chm_seaflux_s == 'GONG_MONAHAN_F') then
!
!    2D fields
      do i = 1, chm_ni
!    for seal-salt
         surfwd(i) = metvar2d(i, MV2D_WSDIAG)
         tdiag(i)  = metvar2d(i, MV2D_TDIAG)
      end do

!     Sea-salt surface flux
      if (local_dbg) then
         write (chm_lun_out, *) 'Compute sea-salt surface flux by cam scheme: ', chm_seaflux_s
      end if

!     Sea-salt surface flux
      call mach_cam_sfss(tdiag, surfwd, rsfrow, fland, chm_ni)

!    Emissions of sea-salt calculated in kg and account for grid area.
!    Change units and divide by dxdy to harmonize with other area emissions before
!    use by diffusion operator
      do k = 1, isize
         nsp_SS = (iae_SS - 1) * isize + k + sp_AERO - 1
         do i = 1, chm_ni
            busper(sm(nsp_SS) % ae_offset + i - 1) = rsfrow(i, k) * &
                                                     metvar2d(i, MV2D_DXDY) * kg2g
         end do
      end do

   end if

! Evaluation of the modulation factor (fmet) to be applied later to the
! fugitive emissions of aerosol species. See subroutine mach_diffusion.ftn90
! for the application of (fmet) prior to vertical diffusion operator
   if (local_dbg) then
      write (chm_lun_out, *) 'Fugitive aerosol emissions modulated by meteorological conditions ', &
                              chm_met_modulation_s
   end if

   if ((chm_met_modulation_s == 'ON') .or. &
       (chm_met_modulation_s == 'CM_ONLY')) then
      do i = 1, chm_ni
         snof = metvar2d(i, MV2D_SNOF)
      ! If the area-averaged grid is predicted to be more than 10% wet
         if (metvar2d(i, MV2D_WSOIL) >= 0.1) then
            busvol(sm(sp_FMET) % out_offset + i - 1) = 0.0
      ! or else allow for fractional modulation, if part of the grid is covered
      ! by snow
         else if (snof > 0.0) then
            busvol(sm(sp_FMET) % out_offset + i - 1) = 1.0 - snof
         else
            busvol(sm(sp_FMET) % out_offset + i - 1) = 1.0
         end if
      end do
   else
      do i = 1, chm_ni
         busvol(sm(sp_FMET) % out_offset + i - 1) = 1.0
      end do
   end if

   return
end subroutine mach_cam_flux
