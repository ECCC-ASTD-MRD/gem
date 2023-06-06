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
! Fichier/File   : mach_gas_main.ftn90
! Creation       : P. Makar, B. Pabla, S. Menard Feb 2007 for GEM-MACH
!                  Janusz Pudykiewicz for CHRONOS 1995
! Description    : Entry point for Gas-phase chemistry process describing NOX/VOC system
!
! Extra info     : - Modification to calculate NON-mass conserving condensable
!                    organic mass for later SOA calculations. Method of Odum
!                    et al., and the two product yield method from Seinfeld's
!                    lab, along with corrections to account for lumping, have
!                    been used. (P. Makar, July 1999).
!                  - Changed test variable for first time step from KT to
!                    IFIRST for initialization of initial total SOA mass (M. Moran, Jan 2000)
!                  - Implement the general form of the instantaneous aerosol yield,
!                    split OC into primary and secondary (Jiang, July 2004)
!
! Arguments:  IN
!               metvar2d -> 2D meteorology variables
!               metvar3d -> 3D meteorology variables
!               step     -> Model step
!               lfu      -> Land-use fractions
!
!             INOUT
!               busvol   -> Volatile bus
!               chem_tr  -> Tracers concentrations (ug/kg)
!
!             LOCAL
!               p2d      -> pressure (pa) on thermodynamic levels
!               rho      -> 3-d air density array (kg/m^3)
!
!============================================================================
!
!!if_on
subroutine mach_gas_main(busper, busvol, chem_tr, metvar2d, metvar3dnocan, &
                         step, ni_nocan, imod2, lfu)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_metvar_mod,       only: MV2D_CANG, MV2D_MT, MV3D_SIGT, MV3D_TPLUS, &
                                   MV3D_HUPLUS, MV3D_ZPLUS, MV3D_ZMOM,        &
                                   MV3D_QCPLUS, MV3D_FTOT, MV2D_PPLUS,        &
                                   MV3D_RHO, MV3D_O3L
   use chm_utils_mod,        only: ik, chm_error_l, chm_timestep, CHM_MSG_DEBUG
   use chm_consphychm_mod,   only: consth, rgasi, mwt_air
   use chm_nml_mod,          only: chm_soa_s, nk_start, chm_strato_s, &
                                   chm_pkg_gas_s, chm_messy_jval_l,   &
                                   chm_debug_3d_i, chm_timings_l,     &
                                   chm_diag_colum_L, chm_mar_halo_l,  &
                                   chm_active_ch4_l
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_O3, sp_CH4, sp_AERO, sp_NWOC, sp_JNO2, &
                                   sp_diff_o3, sp_tend_o3, var_step !, DBG_3D
   use mach_gas_headers_mod, only: mach_soa_jiang,                        &
                                   mach_gas_messy, mach_gas_strato,       &
                                   mach_adom2yb_main, mach_adom2kpp_main, &
                                   mach_adom2kppb_main, mach_saprc_main,  &
                                   mach_marine_halo
   use mach_pkg_gas_mod,     only: gas_species, nspec, nvar, nsp_soa_gases, &
                                   njrxs, ind_ch4, ch4
   use mach_adom2_rates_mod, only: mach_gas_lookup_jvals, ntype2, ntype4, &
                                   mach_adom2_jcorr
   use chm_phyvar_mod,       only: o3lplus
   use linoz_param,          only: p_linoz_tropo, hu_linoz_tropo
   use mach_cam_utils_mod,   only: isize, iae_OC, iae_PC
   use mach_headers_mod,     only: mach_cffeps_calc
   implicit none
!
! Declare subroutine arguments
!
!!if_on
   integer(kind=4), intent   (in) :: step
   integer(kind=4), intent   (in) :: ni_nocan
   integer(kind=4), intent   (in) :: imod2   (ni_nocan)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3dnocan(ni_nocan, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: lfu     (chm_ni, lucprm)
!!if_off
!
!  Declare local variables
!
   real(kind=4)    :: p2d      (ni_nocan, chm_nk)
   real(kind=4)    :: tplus    (ni_nocan, chm_nk)
   real(kind=4)    :: hu_ppm   (ni_nocan, chm_nk)
   real(kind=4)    :: rho      (ni_nocan, chm_nk)
   real(kind=4)    :: sigt     (ni_nocan, chm_nk)
   real(kind=4)    :: zplus    (ni_nocan, chm_nk)
   real(kind=4)    :: zmom     (ni_nocan, chm_nk)
   real(kind=4)    :: cldw     (ni_nocan, chm_nk)
   real(kind=4)    :: cldf     (ni_nocan, chm_nk)
   real(kind=4)    :: skyo     (ni_nocan, chm_nk)
   real(kind=4)    :: csza     (ni_nocan)
   real(kind=4)    :: psurf    (ni_nocan)
   real(kind=4)    :: trppm    (nspec, ni_nocan, chm_nk)
!
! SOA yield calculation related
   real(kind=8)    :: soa      (ni_nocan, chm_nk)
   real(kind=8)    :: poa      (ni_nocan, chm_nk)
   real(kind=8)    :: moi      (ni_nocan, chm_nk)
   real(kind=8)    :: voc_diff (nsp_soa_gases, ni_nocan, chm_nk)
   real(kind=4)    :: dsoa     (ni_nocan, chm_nk)
!
!  conversion factor array to convert 1/kg air to (1/moles air).
!  = 0.08314 * T(kelvin) *rho(kg/m3) / p(mb) note that the constant
! includes conversion
   real(kind=4)    :: conv1(ni_nocan, chm_nk)
!
! MESSY photolysis
   real(kind=8)    :: rjval  (ni_nocan * chm_nk, njrxs)
! o3 column in mol/mol
   real(kind=4)    :: o3vmr2d(ni_nocan, chm_nk)
   real(kind=4)    :: o3_old (ni_nocan, chm_nk)
   real(kind=4)    :: zo3    (ni_nocan, chm_nk)
!
! Look-up table photolysis
   real(kind=8)    :: rjval_lookup(ni_nocan * chm_nk, ntype2+ntype4)

   integer(kind=4) :: i, k, ii, jj, isp, isp_PC, isp_OC
!
! Linoz Related
   real(kind=4)    :: lino3_new(ni_nocan, chm_nk)
   real(kind=4)    :: o3_tend
!
   real(kind=8)    :: hstart   (ni_nocan, chm_nk)
!
!  index remapping 2-d arrays to 1-d
   integer(kind=4) :: this_ik, ik0
!
!  Declare external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp
!
!  CODE BEGINS HERE
!
   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_gas [BEGIN]')
   if (chm_timings_L) call timing_start_omp(330, 'mach_gas', 480)

!* Saved original ozone profile (on resolved model levels), for used by
!  the MESSY photolysis and Linoz
   if (chm_messy_jval_l .or. chm_strato_s == 'LINOZ') then
      do k = 1, chm_nk
         do i = 1, ni_nocan
            ! ozone mmr in ug/kg at T-dT
            o3_old(i, k) = chem_tr(imod2(i), k, sp_O3)
         end do
      end do
   end if
!
   do i = 1, ni_nocan
      ii = imod2(i)
      csza(i)  = metvar2d(ii, MV2D_CANG)
      psurf(i) = metvar2d(ii, MV2D_PPLUS)
   end do
   do k = 1, chm_nk
      do i = 1, ni_nocan
         ii = imod2(i)
         rho(i, k)   = metvar3dnocan(i, k, MV3D_RHO)

         tplus(i, k)  = metvar3dnocan(i, k, MV3D_TPLUS)
         hu_ppm(i, k) = consth * metvar3dnocan(i, k, MV3D_HUPLUS)  ! units ppmv
         sigt(i, k)   = metvar3dnocan(i, k, MV3D_SIGT)
         zplus(i, k)  = metvar3dnocan(i, k, MV3D_ZPLUS) &! Therm. hgt. above ground
                          + max(0.0, metvar2d(ii, MV2D_MT))
!
         p2d(i, k)    = sigt(i, k) * psurf(i)
! rgasi is gas constant, used in conversion of ug/kg to ppmv
         conv1(i, k)  = rgasi * tplus(i, k) * rho(i, k) / p2d(i, k)
!
!  Convert units from mass mixing ratio to ppmv and place into separate array
!  for gas-phase chemical integration.
         do isp = 1, nvar
            trppm(isp, i, k) = chem_tr(ii, k, gas_species(isp)) * &
                               conv1(i, k) / sm(gas_species(isp)) % mol_wt
         end do
         do isp = nvar + 1, nspec
            trppm(isp, i, k) = 0.0
         end do

      end do
   end do
   if (chm_active_ch4_l) then
      do k = 1, chm_nk
         do i = 1, ni_nocan
            ii = imod2(i)
            trppm(ind_ch4, i, k) = chem_tr(ii, k, sp_CH4) * conv1(i, k) / &
                                   sm(sp_CH4) % mol_wt
         end do
      end do
   else
      trppm(ind_ch4, :, :) = ch4
   end if
!
! parameterized marine halogen chemistry
   MARINE_HALO: if (chm_mar_halo_l) then
      jj = findloc(gas_species, sp_O3, dim = 1)
      do k = 1, chm_nk
         do i = 1, ni_nocan
            zo3(i, k) = trppm(jj, i, k)
         end do
      end do
      call mach_marine_halo(zo3, p2d, lfu, imod2, ni_nocan, chm_nk)
      do k = 1, chm_nk
         do i = 1, ni_nocan
            trppm(jj, i, k) = zo3(i, k)
         end do
      end do
   end if MARINE_HALO
!
!
   MESSY_photolysis: if (chm_messy_jval_l) then
            ! dynamic ozone tracer ug/kg -> mole/mole
      o3vmr2d = o3_old * mwt_air / (sm(sp_o3) % mol_wt) * 1.E-9

      if (maxval(o3vmr2d) < 0.1E-10 ) then
         write(*, *) 'Possible error in o3 data for MESSY max(o3vmr) < 0.1E-10'
         chm_error_l = .true.
         return
      end if
!
      call mach_gas_messy(metvar3dnocan, metvar2d, o3vmr2d, rjval, imod2, &
                          ni_nocan, chm_nk)
!
! Skip photolysis cloud correction if MESSY is selected
      skyo = 1.0
!
   else
!
!    determine the cloudy/clear sky correction factor for photolysis rates:
!    cosine of solar zenith angle is positive; sun is above the horizon.
!    calculate photolysis correction factor for clouds.
      do k = 1, chm_nk
         do i = 1, ni_nocan
            zmom(i, k) = metvar3dnocan(i, k, MV3D_ZMOM)   ! Mom. hgt. above ground
            cldf(i, k) = metvar3dnocan(i, k, MV3D_FTOT)   ! Cloud fraction
             ! Cloud liquid water content in g/m3
            cldw(i, k) = metvar3dnocan(i, k, MV3D_QCPLUS) * 1000.0 * rho(i, k)
         end do
      end do
      call mach_adom2_jcorr(zmom, cldw, cldf, csza, skyo, ni_nocan, chm_nk)
!
!     Evaluate photolysis rates from the (ADOM-II) look-up table
      call mach_gas_lookup_jvals(zplus, skyo, csza, rjval_lookup, &
                                 ni_nocan, chm_nk, ni_nocan*chm_nk)
!
      ! Placeholder for non-MESSY photolysis, useful for non-canopy case as well.
      rjval = 1.0d0

      ! Use the lookup table JVAL rates for ADOM2KPP is MESSy_Jval is disabled
      if (trim(chm_pkg_gas_s) == 'ADOM2KPP') rjval = rjval_lookup
   end if MESSY_photolysis
!
   moi = 0.0d0
   if (chm_soa_s /= 'NIL') then
!     Initialize the total organic aerosols
      soa = 0.0d0
      poa = 0.0d0
      do k = 1, chm_nk
         do i = 1, ni_nocan
            ii = imod2(i)
            do jj = 1, isize
               isp_OC = (iae_OC - 1) * isize + jj + sp_AERO - 1
               isp_PC = (iae_PC - 1) * isize + jj + sp_AERO - 1
               soa(i, k) = soa(i, k) + dble(chem_tr(ii, k, isp_OC))
               poa(i, k) = poa(i, k) + dble(chem_tr(ii, k, isp_PC))
               ! add total condensible organics, convert from ug/kg to ug/m3
               moi(i, k) = (soa(i, k) + poa(i, k)) * dble(rho(i, k))
            end do
         end do
      end do
   end if
!
! KPP solver internal time-step
   if (var_step > 0) then
      do k = 1, chm_nk
         do i = 1, ni_nocan
            ii = imod2(i)
            this_ik = ik(ii, k, chm_ni)
            hstart(i, k) = dble(busper(sm(var_step) % per_offset + this_ik))
         end do
      end do
   end if
!
   select case (trim(chm_pkg_gas_s))
      case ('ADOM2')
      !  ADOM-II mechanism with the Young & Boris Solver
         call mach_adom2yb_main(p2d, tplus, hu_ppm, sigt, rjval, trppm, &
                                rjval_lookup, csza, voc_diff, step,     &
                                ni_nocan, chm_nk)

      !  ADOM-II mechanism with the KPP Solver
      case ('ADOM2KPP')
         call mach_adom2kpp_main(p2d, tplus, hu_ppm, sigt, rjval, trppm, &
                                 hstart, voc_diff, ni_nocan, chm_nk)
!
      !  ADOM-II mechanism with the KPP Solver, using rates from adom2_uprate
      case ('ADOM2KPPB')
         call mach_adom2kppb_main(p2d, tplus, hu_ppm, sigt, rjval, trppm,   &
                                  hstart, rjval_lookup, voc_diff, ni_nocan, &
                                  chm_nk)
!
      !  SAPRC07C mechanisms with the KPP Solver
      case ('SAPRC07C', 'SAPRC07CS')
         call mach_saprc_main(p2d, tplus, hu_ppm, rjval, trppm, hstart, &
                              moi, dsoa, ni_nocan, chm_nk)
   end select

   if (chm_error_l) return
!
!============================================================================
! Return updated gas concentrations into the dynamic bus
! (Convert units back from ppmv to mass mixing ratio)
!============================================================================
   do k = nk_start, chm_nk
      do i = 1, ni_nocan
         ii = imod2(i)
         do isp = 1, nvar
            chem_tr(ii, k, gas_species(isp)) = trppm(isp, i, k) / &
                                 conv1(i, k) * sm(gas_species(isp)) % mol_wt
         end do
         if (chm_active_ch4_l) then
            chem_tr(ii, k, sp_CH4) = trppm(ind_ch4, i, k) / &
                                     conv1(i, k) * sm(sp_CH4) % mol_wt
         end if
      end do
   end do
!
! Return the KPP solver internal time-step to the bus
   if (var_step > 0) then
      do k = 1, chm_nk
         do i = 1, ni_nocan
            ii = imod2(i)
            this_ik = ik(ii, k, chm_ni)
            busper(sm(var_step) % per_offset + this_ik) = real(hstart(i, k))
         end do
      end do
   end if
!
 !---- write Jx value to file ---------
   if (sp_JNO2 > 0) then
      ! Output JNO2 at the lowest model layer
      do i = 1, ni_nocan
         do k = 1, chm_nk
            ik0 = ik(i+1, k, ni_nocan)
            busvol(sm(sp_JNO2) % out_offset + imod2(i) - 1) = real(rjval(ik0, 1))
         end do
      end do
   end if
!
!!!  need to set "chm_debug_3d_i = 6" in gem_settings.nml and
!!!  add "sortie_p([3DB1,3DB2,3DB3,3DB4,3DB5,3DB6])" in outcfg.out
!   if (chm_debug_3d_i >= 6) then
!      do i = 1, ni_nocan
!         ii = imod2(i)
!         do k = 1, chm_nk
!            this_ik = ik(ii, k, chm_ni)
!            ik0 = ik(i+1, k, ni_nocan)
!            busvol(sm(dbg_3d(1)) % out_offset + this_ik) = real(rjval(ik0, 1))  !jcor JNO2
!            busvol(sm(dbg_3d(2)) % out_offset + this_ik) = real(rjval(ik0, 5))  !jcor JO3=JO1d
!            busvol(sm(dbg_3d(3)) % out_offset + this_ik) = real(rjval(ik0, 2))  !jcor JNO3
!            busvol(sm(dbg_3d(4)) % out_offset + this_ik) = real(rjval(ik0, 7))  !jcor JHNO3
!            busvol(sm(dbg_3d(5)) % out_offset + this_ik) = real(rjval(ik0, 9))  !jcor JH2O2
!            busvol(sm(dbg_3d(6)) % out_offset + this_ik) = real(rjval(ik0, 4))  !jcor JO3p
!         end do
!      end do
!   end if
!
!==============================================================================
! SOA yield calculation
!==============================================================================
!
   if (chm_soa_s /= 'NIL') then
! Calculate the change in concentration of soa using JIANG et al. yield approach
      if (chm_soa_s == 'JIANG' .and. chm_pkg_gas_s(1:4) == 'ADOM') then
         call mach_soa_jiang(moi, dsoa, voc_diff, tplus, p2d, ni_nocan, chm_nk)
      end if
!
      do k = nk_start, chm_nk
         do i = 1, ni_nocan
        ! convert dsoa from ug/m3 to ug/kg and output to vol bus
            busvol(sm(sp_NWOC) % out_offset + ik(imod2(i), k, chm_ni)) = &
                                                    dsoa(i, k) / rho(i, k)
         end do
      end do
!
   end if
!
!==============================================================================
! Start of Gas Phase Stratospheric Chemistry
!==============================================================================
!
   if (chm_strato_s == 'LINOZ') then
!
      if (o3lplus > 0) then
         lino3_new(:, :) = metvar3dnocan(:, :, MV3D_O3L)
      else
!      Evaluate LINOZ O3 tendency locally
         call mach_gas_strato(busper, hu_ppm, tplus, sigt, o3_old, psurf, &
                              lino3_new, ni_nocan, chm_nk, imod2)
      end if
!
!     Apply LINOZ tendencies to ozone (O3) above hu_linoz_tropo specific humidity or
!     for pressure lower than p_linoz_tropo,
      do k = 1, chm_nk
         do i = 1, ni_nocan
            if ((hu_ppm(i, k) < hu_linoz_tropo) .or. &
                (p2d(i, k) < p_linoz_tropo)) then
               ii = imod2(i)
               chem_tr(ii, k, sp_O3) = lino3_new(i, k)
            end if
         end do
      end do
!
!    Output diagnostics
!
      IF_DIAG: if (chm_diag_colum_L) then
         do k = 1, chm_nk
            do i = 1, ni_nocan
               ii = imod2(i)
               this_ik = ik(ii, k, chm_ni)
         ! Ozone tendency
               o3_tend = (chem_tr(ii, k, sp_O3) - o3_old(i, k)) / chm_timestep    !ug /kg /sec
               busvol(sm(sp_tend_o3) % out_offset + this_ik) = o3_tend
               busvol(sm(sp_diff_o3) % out_offset + this_ik) = o3_tend
            end do
         end do

      end if IF_DIAG
!
   end if ! End chm_strato_s

   call msg_toall(CHM_MSG_DEBUG, 'mach_gas [END]')
   if (chm_timings_L) call timing_stop_omp(330)
   !-----------------------------------------------------------------
   return
end subroutine mach_gas_main
