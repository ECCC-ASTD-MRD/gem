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
! Fichier/File   : mach_main.ftn90
! Creation       : S. Menard  ,  GEM-MACH, Feb 2007.
! Description    : GEM-MACH chemistry processes are called from mach_main.ftn90. Each chemistry
!                  process involved from mach_main.ftn90 need to have it's own key in
!                  gem_settings.nml namelist.
!
! Extra info     :  Type of bus
! -----------
!                     ** permanant bus
!                           variables on this bus keep constant values throughout model integration
!                           steps unless you change those at
!                           any model time step
!                     ** volatile bus
!                           variables stay throughout model time steps, but their values change
!                           every time step e.g.
!
! Arguments:
!            IN
!               slab_index   -> Slice number
!               step         -> Timestep number
!               ni_can       -> Number of horizontal grid cells with canopy
!               ni_nocan     -> Number of horizontal grid cells without canopy
!
!           IN/OUT
!               chem_tr      -> Chemical tracers concentrations (ug/kg)
!               busper       -> Permanent bus for physics
!               busvol       -> Volatile bus for physics
!
!===============================================================================
!
!!if_on
subroutine mach_main(busper, busvol, chem_tr, metvar2d, metvar3d, &
                     slab_index, step, ni_can, ni_nocan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
!!if_off
   use chm_utils_mod,        only: global_debug, chm_lun_out, chm_error_l, ik
   use chm_nml_mod,          only: chm_biog_s, chm_gas_drydep_s,           &
                                   chm_do_mjpts_l, chm_mj_treatment_s,     &
                                   chm_diag_accum_L, chm_diag_colum_L,     &
                                   chm_sat_seasons_l, chm_do_cffeps_l,     &
                                   chm_pkg_gas_s, chm_pkg_pm_s, chm_vert_diff_s
   use chm_species_idx_mod,  only: sp_SO4
   use mach_drydep_mod,      only: lucprm
   use mach_pkg_gas_mod,     only: num_be_sp
   use mach_headers_mod,     only: mach_input_check, mach_calc_diag  , &
                                   mach_calc_season, mach_landuse    , &
                                   mach_plumerise  , mach_biog_main  , &
                                   mach_diff_prep  , mach_output     , &
                                   mach_lai_adjust , mach_cffeps_main, &
                                   mach_canopy_levels, mach_canopy_transfer
   use mach_gas_headers_mod, only: mach_gas_drydep_main, mach_gas_main, &
                                   mach_gas_canopy
   use mach_cam_headers_mod, only: mach_cam_flux, mach_pm_chem, gocart_so2so4
   use chm_ptopo_grid_mod,   only: pm_nk, pm_nkc, nkt, nkc

   use mach_pkg_tendencies_mod

   implicit none
!
!  Declaration of subroutine arguments
!
!!if_on
   integer(kind=4), intent   (in) :: slab_index
   integer(kind=4), intent   (in) :: step
   integer(kind=4), intent   (in) :: ni_can, ni_nocan
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr(chm_ni, chm_nk + 1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
!!if_off
!
!  Declaration of local variables.
!
!  * Seasons index
   integer(kind=4) :: iseasn(chm_ni)
!  * Land-use fractions modified by meteorological conditions (snow cover, etc)
   real(kind=4)    :: landuse(chm_ni, lucprm)
   integer(kind=4) :: i, k, ic, kc, sp
   logical(kind=4) :: local_dbg

   integer(kind=4)   :: kmod(ni_can, chm_nk), kcan(ni_can, nkc)
   integer(kind=4)   :: kmod2(ni_nocan, chm_nk)
   integer(kind=4)   :: imod(ni_can), imod2(ni_nocan)
!
   real(kind=4)      :: oldso4_nocan(ni_nocan, chm_nk)
   real(kind=4)      :: oldso4_can(ni_can, nkt)
   real(kind=4)      :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   real(kind=4)      :: metvar3dnocan(ni_nocan, chm_nk, SIZE_MV3D)
   real(kind=4)      :: emisbio_can(ni_can, nkt, num_be_sp)
   real(kind=4)      :: tracers_can(ni_can, nkc, nb_dyn_tracers)

!===============================================================================
! Code statements begin here
!===============================================================================
!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
!
   if (local_dbg) then
      write (chm_lun_out, *) 'in mach_main'
      write (chm_lun_out, *) 'slab_index, step ', slab_index, step
      write (chm_lun_out, *) 'chm_ni, chm_nk ', chm_ni, chm_nk
   end if
!
!================================================================================
! Start of landuse and season calculations
!================================================================================
!
   call mach_landuse(busper, metvar2d, landuse)
   if (chm_error_l) return
!
   if (step <  1) then
      if (chm_sat_seasons_l) then
         call mach_lai_adjust(busper, landuse, slab_index)
         if (chm_error_l) return
      end if
!
      call mach_output(busper, busvol, chem_tr, metvar2d, metvar3d, landuse)
!
      if (local_dbg) then
         write (chm_lun_out, *) 'Initialization for timestep 0 completed. Exiting mach_main '
      end if
      return
   end if
!
   call mach_calc_season(metvar2d, iseasn)
!
!===============================================================================
! Start of dry deposition for gases
!===============================================================================
!
   select case (chm_gas_drydep_s)
!
      case ('ROBICHAUD', 'ROBICHAUD2', 'ROBICHAUD3')
         if (local_dbg) then
            write (chm_lun_out, *) 'Compute the dry deposition for gas: ', chm_gas_drydep_s
         end if
         call mach_gas_drydep_main(busper, busvol, metvar2d, landuse, iseasn)
!
      case default
         if (local_dbg) then
            write (chm_lun_out, *) '> Warning '
            write (chm_lun_out, *) '> No dry deposition for gas: ', chm_gas_drydep_s
         end if
!
   end select
!
!===============================================================================
! Start of major points sources emissions injection
!===============================================================================
!
   if (chm_do_mjpts_l) then
      select case (chm_mj_treatment_s)
!
         case ('PLUMERISE', 'PLUMERISE2')
            if (local_dbg) then
               write (chm_lun_out, *) 'Compute the major points plumerise ', chm_mj_treatment_s, &
                                      ' chm_do_mjpts_l ', chm_do_mjpts_l
            end if
            call mach_plumerise(busvol, chem_tr, metvar2d, metvar3d, landuse, &
                                slab_index)
            if (chm_error_l) return
!
         case default
            if (local_dbg) then
               write (chm_lun_out, *) '> Warning '
               write (chm_lun_out, *) '> No major point source treatment: ', chm_mj_treatment_s
            end if
!
      end select
   end if
!
   if (chm_do_cffeps_l) then
      call mach_cffeps_main(chem_tr, metvar2d, metvar3d, slab_index, step)
      if (chm_error_l) return
   end if
!
!===============================================================================
! Set up canopy sub-model
!===============================================================================
!
   if (ni_can > 0) then
      call mach_canopy_levels(busper, metvar3dcan, metvar3dnocan, metvar3d, &
                         metvar2d, imod, imod2, kmod, kcan, ni_can, ni_nocan)

      if (chm_error_l) return
   else
      metvar3dnocan = metvar3d
      imod  = 0
      imod2 = (/ (i, i = 1, chm_ni) /)
      metvar3dcan = 0.0
      kcan        = 0
      kmod        = 0
   end if
!
!===============================================================================
! Start of the biogenic emissions scheme
!===============================================================================
!
   if (trim(chm_biog_s) /= 'NIL') then
      if (local_dbg) then
         write (chm_lun_out, *) 'Compute the biogenic emissions: ', chm_biog_s
      end if
      call mach_biog_main(busper, metvar2d, metvar3d, iseasn, &
                          metvar3dcan, emisbio_can, kcan, ni_can)
      if (chm_error_l) return
!
   else
      if (local_dbg) then
         write (chm_lun_out, *) '> Warning '
         write (chm_lun_out, *) '> No biogenic emissions: ', chm_biog_s
      end if
      emisbio_can = 0.0
!
   end if
!
!===============================================================================
! Sea-salt emissions and meteorological modulation of fugitive dust
!===============================================================================
!
   if (chm_pkg_pm_s(1:3) == 'CAM') then
      call mach_cam_flux(busper, busvol, metvar2d, landuse)
   end if
!
!===============================================================================
! Distribute tracer concentration from model resolved layers into canopy layers
!===============================================================================
   if (ni_can > 0) then
      call mach_canopy_transfer(chem_tr, tracers_can, metvar2d, metvar3d, &
                                metvar3dcan, kmod, kcan, imod, ni_can, 0)
   end if
!
!===============================================================================
! Start of the vertical diffusion
!===============================================================================
!
   if (chm_vert_diff_s /= 'NIL') then
      if (local_dbg) then
         write (chm_lun_out, *) 'calling diffusion ', chm_vert_diff_s
         write (chm_lun_out, *) 'setting min and max on vertical diffusion coefficient'
      end if
      call mach_input_check (busvol, metvar2d, metvar3d, landuse)
!
#if defined(MACH_TENDENCIES)
      call tendency_store(tend_diffusion, chem_tr, busvol)
#endif
!
! Calling diffusion for all species
      call mach_diff_prep(busper, busvol, chem_tr, metvar2d, metvar3d, &
                          metvar3dnocan, metvar3dcan, emisbio_can, tracers_can, &
                          kmod, kcan, imod, imod2, ni_can, ni_nocan)
!
#if defined(MACH_TENDENCIES)
      call tendency_delta(tend_diffusion, chem_tr, busvol)
#endif
   else
      if (local_dbg) write (chm_lun_out, *) '> Warning: No mach_diffusion'
   end if
!
!===============================================================================
! Start of Gas Phase Chemistry (only between level (nk_start and level chm_nk)
!===============================================================================
!
!  Keep [SO4]g before gas chemistry for use in CAM
!
   if (chm_pkg_pm_s(1:3) == 'CAM') then
      do k = 1, chm_nk
         do i = 1, ni_nocan
            oldso4_nocan(i, k) = chem_tr(imod2(i), k, sp_SO4)
            kmod2(i, k) = k
         end do
      end do
      if (ni_can > 0) then
         do k = 1, chm_nk
            do ic = 1, ni_can
               oldso4_can(ic, kmod(ic, k)) = chem_tr(imod(ic), k, sp_SO4)
            end do
         end do
         do kc = 1, nkc
            do ic = 1, ni_can
               oldso4_can(ic, kcan(ic, kc)) = tracers_can(ic, kc, sp_SO4)
            end do
         end do
      end if
   end if
!
   if (chm_pkg_gas_s /= 'NIL') then
!
      if (local_dbg) then
         write (chm_lun_out, *) 'Compute the gas phase chemistry: ', chm_pkg_gas_s
      end if
#if defined(MACH_TENDENCIES)
      call tendency_store(tend_gaschem, chem_tr, busvol)
#endif
!
      if (ni_nocan > 0) then
         call mach_gas_main(busper, busvol, chem_tr, metvar2d, metvar3dnocan, &
              step, ni_nocan, imod2, landuse)
         if (chm_error_l) return
      end if
      if (ni_can > 0) then
         call mach_gas_canopy(busper, busvol, chem_tr, tracers_can, metvar2d, &
              metvar3dcan, step, ni_can, imod, kcan, kmod, landuse)
         if (chm_error_l) return
      end if
#if defined(MACH_TENDENCIES)
      call tendency_delta(tend_gaschem, chem_tr, busvol)
#endif
!
   else
      if (local_dbg) then
         write (chm_lun_out, *) '> Warning '
         write (chm_lun_out, *) '> No gas phase chemistry: ', chm_pkg_gas_s
      end if
!
   end if
!
!===============================================================================
! Start of Aerosol Processes
!===============================================================================
!
   select case (chm_pkg_pm_s)
!
       case ('CAM12BINS', 'CAM2BINS')
          if (local_dbg) then
             write (chm_lun_out, *) 'Compute processes for CAM aerosol scheme: ', chm_pkg_pm_s
          end if
#if defined(MACH_TENDENCIES)
          call tendency_store(tend_pmchem, chem_tr, busvol)
#endif
!  Columns without vegetation canopies
          if (ni_nocan > 0) then
             call mach_pm_chem(busvol, busper, chem_tr, metvar2d, metvar3dnocan,&
                  landuse, oldso4_nocan, iseasn, step, slab_index, ni_nocan, &
                  pm_nk, chm_nk, imod2, kmod2)
             if (chm_error_l) return
          end if
!  Columns with vegetation canopies
          if (ni_can > 0) then
             call mach_pm_chem(busvol, busper, chem_tr, metvar2d, metvar3dcan, &
                  landuse, oldso4_can, iseasn, step, slab_index, ni_can, &
                  pm_nkc, nkt, imod, kmod, tracers_can, kcan)
             if (chm_error_l) return
          end if
#if defined(MACH_TENDENCIES)
          call tendency_delta(tend_pmchem, chem_tr, busvol)
#endif
!
       case ('GOCART_SO2')
          if (local_dbg) then
             write (chm_lun_out, *) 'Simple COGART process of SO2 oxydation ', chm_pkg_pm_s
          end if
          call gocart_so2so4(chem_tr, metvar3d)
!
       case default
          if (local_dbg) then
             write (chm_lun_out, *) '> Warning '
             write (chm_lun_out, *) '> No aerosol process: ', chm_pkg_pm_s
          end if
!
   end select
!
!===============================================================================
!   For all enabled chemistry tracers, copy prognostic level (chm_nk)
!   into diagnostic level (chm_nk + 1)
!===============================================================================
!
   if (ni_nocan > 0) then
      do sp = 1, nb_dyn_tracers
         do i = 1, ni_nocan
            chem_tr(imod2(i), chm_nk + 1, sp) = chem_tr(imod2(i), chm_nk, sp)
         end do
      end do
   end if

   if (ni_can > 0) then
      call mach_canopy_transfer(chem_tr, tracers_can, metvar2d, metvar3d, &
                                metvar3dcan, kmod, kcan, imod, ni_can, 1)
      if (chm_error_l) return
   end if
!
   if (chm_diag_accum_L .or. chm_diag_colum_L) then
!
      call mach_calc_diag(busvol, busper, chem_tr, metvar2d, metvar3d, step)
   end if
!
!===============================================================================
!  Units conversion for a subset of PM and gases species before outputting them.
!===============================================================================
!
   call mach_output (busper, busvol, chem_tr, metvar2d, metvar3d, landuse)
!
!===============================================================================

return
end subroutine mach_main
