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
! Fichier/File   : mach_plumerise.ftn90
! Creation       : H. Landry, A. Kallaur, S. Menard, March 2008 adaptation for GEM-MACH
!                  Wanmin Gong, for AURAMS 2004
!                  Original version from Janusz Pudykiewicz for CHRONOS 1995
!
! Description    : Assembly of major-point-source emission data for calculation
!                  with Eulerian transport model simulating atmospheric
!                  oxidants. Plume rise formula based on Briggs (1984).
! Arguments:
!            IN
!                slab_index    -> Slab number
!                busvol        -> Physics volatile bus.
!                metvar{2d,3d} -> met. fields (chemistry lib only, see chm_exe)
!                                 copied from Physics buses (see chm_load_metvar).
!
!             IN/OUT
!                chem_tr       -> Chemical species concentrations (ug/kg).
!
! Extra info     :
!
!============================================================================
!
!!if_on
subroutine mach_plumerise(busvol, chem_tr, metvar2d, metvar3d, landuse, &
                          slab_index)
   use chm_metvar_mod,          only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,      only: chm_ni, chm_nk
   use chm_species_info_mod,    only: nb_dyn_tracers
   use mach_drydep_mod,         only: lucprm
!!if_off
   use chm_nml_mod,             only: wf_case, chm_timings_L
   use chm_metvar_mod,          only: MV2D_ILMO, MV2D_H, MV2D_DXDY, MV2D_UE, &
                                      MV2D_TDIAG, MV2D_WSDIAG,               &
                                      MV3D_TPLUS, MV3D_WS, MV3D_ZMOM, MV3D_RHO
   use chm_utils_mod,           only: chm_lun_out, global_debug, ik,         &
                                      chm_error_l, CHM_MSG_DEBUG, chm_timestep
   use chm_species_info_mod,    only: species_master
   use chm_mjrpts_sortinfo_mod, only: lnb_sources, lstack_info, Lstack_emis,   &
                                      me_species_index, nb_me_species, i_gilc, &
                                      i_gjlc, i_typ
   use mach_headers_mod,        only: mach_plumerise_weight, &
                                      mach_plumerise_weight4fire
   use chm_phyvar_mod,          only: shear2, rig
   implicit none
!
!!if_on
   integer(kind=4), intent   (in) :: slab_index
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: landuse (chm_ni, lucprm)
!!if_off
   external msg_toall, timing_start_omp, timing_stop_omp
!
!  Local variables
!========================
!
!  Automatic arrays
!
   integer(kind=4)                   :: index_above_pbl
   ! Weight of each cell in the column over the point source
   real(kind=4), dimension(chm_nk)   :: weight
   ! Safe inverse of monin-obukhov_length used only to avoid the possibility
   ! of a division by zero
   real(kind=4)                      :: safe_inv_mo_length
   ! height of the momentum model levels in meter above ground
   real(kind=4), dimension(chm_nk+1) :: z_magl
   ! Air density
   real(kind=4), dimension(chm_nk)   :: rho
   ! Emission rate to concentration factor
   real(kind=4) , dimension(chm_nk)  :: rate_factor
   ! Area of a grid cell (m^2)
   real(kind=4)                      :: cell_area
   ! Boundary layer depth (m)
   real(kind=4)                      :: pbl_hgt
   ! Surface frition velocity
   real(kind=4)                      :: ustar
   ! Temperature vertical profile
   real(kind=4), dimension(chm_nk+1) :: tplus
   ! Wind speed magnitude profile
   real(kind=4), dimension(chm_nk+1) :: wnd_spd
   ! Atmospheric stability parameter
   real(kind=4), dimension(chm_nk)   :: stab
   ! Vegetation fraction
   real(kind=4), dimension(4)        :: F_bio
   ! Atmospheric stability parameter
   real(kind=4), dimension(:, :), pointer :: zrig, zshear2
!
!  Scalar variables
!
   logical(kind=4) :: local_dbg
   integer(kind=4) :: isp
   integer(kind=4) :: k, cur_source, stack_i, i_me, i_prev  ! loop indices
   real(kind=4)    :: emiss_rate
!
! Constants
!
   real(kind=4), parameter :: GRAMMES_TO_MICROGRAMMES = 1.0e+06
   real(kind=4), parameter :: THRESHOLD = 1.e06

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_plumerise [BEGIN]')
   if (chm_timings_L) call timing_start_omp(315, 'mach_plumerise', 480)
!
   if (local_dbg) then
      write(chm_lun_out, *) ' '
      write(chm_lun_out, *) '*************************************************'
      write(chm_lun_out, *) 'ENTER MACH PLUMERISE SUBROUTINE'
      write(chm_lun_out, *) 'Processing slab number: ', slab_index
      write(chm_lun_out, *) 'Number of sources is lnb_sources: ', lnb_sources
   end if

!
   nullify(zrig, zshear2)
   zrig(1:chm_ni, 1:chm_nk)    => busvol(rig:)
   zshear2(1:chm_ni, 1:chm_nk) => busvol(shear2:)
!
! Loop over the total number of stacks (point sources) found:
! calculate the weight factor for the entire column for each stack,
! and then apply it to the emission rate of all the species.
!
   i_prev = 0
   do cur_source = 1, lnb_sources
      if (nint(lstack_info(i_gjlc, cur_source)) == slab_index) then
         stack_i = nint(lstack_info(i_gilc, cur_source))
      else
         cycle
      end if

      if (i_prev /= stack_i) then
         cell_area = metvar2d(stack_i, MV2D_DXDY)
         pbl_hgt   = metvar2d(stack_i, MV2D_H)
         ustar     = metvar2d(stack_i, MV2D_UE)
!
      ! Retrieve height of the model level interface in meter above ground-level
         do k = 1, chm_nk
            z_magl(k)  = metvar3d(stack_i, k, MV3D_ZMOM)
            rho(k)     = metvar3d(stack_i, k, MV3D_RHO)
            tplus(k)   = metvar3d(stack_i, k, MV3D_TPLUS)
            wnd_spd(k) = metvar3d(stack_i, k, MV3D_WS)
            stab(k)    = zrig(stack_i, k) * zshear2(stack_i, k)
         end do
         z_magl(chm_nk + 1)  = 0.0
         tplus(chm_nk + 1)   = metvar2d(stack_i, MV2D_TDIAG)
         wnd_spd(chm_nk + 1) = metvar2d(stack_i, MV2D_WSDIAG)

         if (abs(metvar2d(stack_i, MV2D_ILMO)) > THRESHOLD) then
            safe_inv_mo_length = sign(THRESHOLD, metvar2d(stack_i, MV2D_ILMO))
         else
            safe_inv_mo_length = metvar2d(stack_i, MV2D_ILMO)
         end if

    ! Evaluate position of boundary layer from the ground up, for this source
         do k = chm_nk, 1, -1
            index_above_pbl = k
            if (z_magl(k) >= pbl_hgt) then
               exit
            end if
         end do

         i_prev = stack_i
      end if
!
! Calculate values for weight function (in the vertical, values 0.0<=w(k)<=1.0,
! with sigma(k=1,chm_nk)w(k)=1.0)
! i_typ = flag for forest fire (fire=500, others=100)
! wf_case=0 (Brigg's plumerise) cannot be applied to CFFEPS fire emissions,
! but is utilized for antropogenic and FEPS fire emissions
!
      weight = 0.0
      if ((lstack_info(i_typ,cur_source) < 500.) .or. (wf_case == 0)) then
         call mach_plumerise_weight(cur_source, z_magl, index_above_pbl, &
                                    safe_inv_mo_length, weight, pbl_hgt, &
                                    ustar, tplus, wnd_spd, stab)
         if (chm_error_l) return
      else
!
! landuse 1 boreal forest
! landuse 2 temperate forest
! landuse 3 shrubs
! landuse 4 grass and cropland
         F_bio(1) = landuse(stack_i, 1)
         F_bio(2) = landuse(stack_i, 4) + landuse(stack_i, 5)
         F_bio(3) = landuse(stack_i, 9) + landuse(stack_i, 10)
         F_bio(4) = landuse(stack_i, 6) + landuse(stack_i, 7)
!
         call mach_plumerise_weight4fire(cur_source, z_magl, index_above_pbl, &
                                         pbl_hgt, safe_inv_mo_length, weight, &
                                         F_bio, rho)
      end if

      if (local_dbg) then
!     if(mod(cur_source,100) == 0) then
         write(chm_lun_out, *) 'current source number : ', cur_source
         write(chm_lun_out, *) 'current source coordinates ', stack_i, slab_index
         write(chm_lun_out, *) 'height of momentum levels  ', z_magl
         write(chm_lun_out, *) 'current source weight      ', weight
!
         if (any(weight < 0.0)) then
            write(0, *) '### Error in mach_plumerise ###'
            write(0, *) '# Negative weight major point emission', weight
            write(0, *) '###         ABORT         ###'
            chm_error_l = .true.
            return
         end if
      end if

      do k = 1, chm_nk
         rate_factor(k) = GRAMMES_TO_MICROGRAMMES * weight(k) * chm_timestep / &
                           (cell_area * rho(k))
      end do

      do i_me = 1, nb_me_species

!  Apply plumerise injection for species to the chemistry dynamical bus

         emiss_rate = Lstack_emis(i_me, cur_source)

         if (emiss_rate > 0.0) then
            isp = me_species_index(i_me)

            if (local_dbg) then
!              if (mod(cur_source,100) == 0) then
               write(chm_lun_out, *) 'FOUND SPECIES {DYN,mjpt} NAME : ', &
                                      species_master(isp) % dyn_name,    &
                                      species_master(isp) % me_name
               write(chm_lun_out, *) 'current source value  : ', emiss_rate
            end if

            do k = 1, chm_nk
               chem_tr(stack_i, k, isp) = chem_tr(stack_i, k, isp) + &
                                          emiss_rate * rate_factor(k)
            end do
         end if
      end do
!
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_plumerise [END]')
   if (chm_timings_L) call timing_stop_omp(315)
   !-----------------------------------------------------------------

   return
end subroutine mach_plumerise
