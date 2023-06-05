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
!               busvol      -> Volatile bus
!               chem_tr     -> Tracers concentrations (ug/kg)
!               tracers_can -> Tracers concentrations within the canopy
!
!             LOCAL
!               p2d      -> pressure (pa) on thermodynamic levels
!               rho      -> 3-d air density array (kg/m^3)
!
!============================================================================
!
!!if_on
subroutine mach_gas_canopy(busper, busvol, chem_tr, tracers_can, metvar2d, &
                           metvar3dcan, step, ni_can, imod, kcan, kmod, lfu)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
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
   use chm_species_idx_mod,  only: sp_O3, sp_CH4, sp_AERO, sp_NWOC, sp_diff_o3, &
                                   sp_tend_o3, sp_CRL, sp_LAI, sp_CLU, sp_HC,   &
                                   sp_JNO2, var_step !, dbg_3d
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
   implicit none
!
! Declare subroutine arguments
!
!!if_on
   integer(kind=4), intent   (in) :: step
   integer(kind=4), intent   (in) :: ni_can
   integer(kind=4), intent   (in) :: imod(ni_can)
   integer(kind=4), intent   (in) :: kcan(ni_can, nkc), kmod(ni_can, chm_nk)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent(inout) :: tracers_can(ni_can, nkc, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   real(kind=4),    intent   (in) :: lfu     (chm_ni, lucprm)
!!if_off
!
!  Declare local variables
!
   real(kind=4)    :: p2d      (ni_can, nkt)
   real(kind=4)    :: tplus    (ni_can, nkt)
   real(kind=4)    :: hu_ppm   (ni_can, nkt)
   real(kind=4)    :: rho      (ni_can, nkt)
   real(kind=4)    :: sigt     (ni_can, nkt)
   real(kind=4)    :: zplus    (ni_can, nkt)
   real(kind=4)    :: zmom     (ni_can, nkt)
   real(kind=4)    :: cldw     (ni_can, nkt)
   real(kind=4)    :: cldf     (ni_can, nkt)
   real(kind=4)    :: skyo     (ni_can, nkt)
   real(kind=4)    :: csza     (ni_can)
   real(kind=4)    :: psurf    (ni_can)
   real(kind=4)    :: trppm    (nspec, ni_can, nkt)
!
! SOA yield calculation related
   real(kind=8)    :: soa      (ni_can, nkt)
   real(kind=8)    :: poa      (ni_can, nkt)
   real(kind=8)    :: moi      (ni_can, nkt)
   real(kind=8)    :: voc_diff (nsp_soa_gases, ni_can, nkt)
   real(kind=4)    :: dsoa     (ni_can, nkt)
!
! conversion factor array to convert 1/kg air to (1/moles air).
! = 0.08314 * T(kelvin) *rho(kg/m3) / p(mb) note that the constant
! includes conversion
   real(kind=4)    :: conv1(ni_can, nkt)
!
! MESSY photolysis
   real(kind=8)    :: rjval(ni_can * nkt, njrxs)
!
! o3 column in mol/mol
   real(kind=4)    :: o3vmr2d(ni_can, nkt)
   real(kind=4)    :: o3_old (ni_can, nkt)
   real(kind=4)    :: zo3    (ni_can, nkt)
!
! Look-up table photolysis
   real(kind=8)    :: rjval_lookup(ni_can * nkt, ntype2+ntype4)
!
   real(kind=8)    :: hstart(ni_can, nkt)

   integer(kind=4) :: k, jj, isp, isp_PC, isp_OC
   integer(kind=4) :: ii, ic, kc, kk
!
! Linoz Related
   real(kind=4)    :: lino3_new(ni_can, nkt)
   real(kind=4)    :: o3_tend
!
! index remapping 2-d arrays to 1-d
   integer(kind=4) :: this_ik, ik0

   real(kind=4)    :: coeff, klight
   real(kind=4)    :: laicumfrac(4)
   real(kind=4)    :: corr_bio(6), zbio(6)
   real(kind=4)    :: corr(ni_can, nkt), corrmom(nkt + 1)
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
!
   do ic = 1, ni_can
      ii = imod(ic)
      csza(ic)  = metvar2d(ii, MV2D_CANG)
      psurf(ic) = metvar2d(ii, MV2D_PPLUS)
   end do
   do k = 1, nkt
      do ic = 1, ni_can
         ii = imod(ic)

         rho(ic, k)    = metvar3dcan(ic, k, MV3D_RHO)
         tplus(ic, k)  = metvar3dcan(ic, k, MV3D_TPLUS)
         hu_ppm(ic, k) = consth * metvar3dcan(ic, k, MV3D_HUPLUS)      ! units ppmv
         sigt(ic, k)   = metvar3dcan(ic, k, MV3D_SIGT)
         zplus(ic, k)  = metvar3dcan(ic, k, MV3D_ZPLUS) &! Therm. hgt. above ground
                          + max(0.0, metvar2d(ii, MV2D_MT))
         zmom(ic, k) = metvar3dcan(ic, k, MV3D_ZMOM)     ! Mom. hgt. above ground
!
         p2d(ic, k)   = sigt(ic, k) * psurf(ic)
! rgasi is gas constant, used in conversion of ug/kg to ppmv
         conv1(ic, k) = rgasi * tplus(ic, k) * rho(ic, k) / p2d(ic, k)
      end do
   end do
!
!  Convert units from mass mixing ratio to ppmv and place into separate array
!  for gas-phase chemical integration.
   trppm = 0.0
   do ic = 1, ni_can
      ii = imod(ic)
!
! (i): Model resolved layers:
      do kk = 1, chm_nk
         k = kmod(ic, kk)
         do isp = 1, nvar
            trppm(isp, ic, k) = chem_tr(ii, kk, gas_species(isp)) * &
                               conv1(ic, k) / sm(gas_species(isp)) % mol_wt
         end do
         if (chm_active_ch4_l) then
            trppm(ind_ch4, ic, k) = chem_tr(ii, kk, sp_CH4) * &
                                    conv1(ic, k) / sm(sp_CH4) % mol_wt
         else
            trppm(ind_ch4, ic, k) = ch4
         end if
!
         ! Old O3 On non-canopy grids.
         o3_old(ic, k) = chem_tr(ii, kk, sp_O3)
      end do
!
! (ii): Canopy shaded layers:
      do kc = 1, nkc
         k = kcan(ic, kc)
         do isp = 1, nvar
            jj = gas_species(isp)
            trppm(isp, ic, k) = tracers_can(ic, kc, jj) * conv1(ic, k) / &
                                sm(jj) % mol_wt
         end do
         if (chm_active_ch4_l) then
            trppm(ind_ch4, ic, k) = tracers_can(ic, kc, sp_CH4) * &
                                    conv1(ic, k) / sm(sp_CH4) % mol_wt
         else
            trppm(ind_ch4, ic, k) = ch4
         end if
!
         o3_old(ic, k) = tracers_can(ic, kc, sp_O3)
      end do
!
   end do
!
! parameterized marine halogen chemistry
   MARINE_HALO: if (chm_mar_halo_l) then
      jj = findloc(gas_species, sp_O3, dim = 1)
      do k = 1, nkt
         do ic = 1, ni_can
            zo3(ic, k) = trppm(jj, ic, k)
         end do
      end do
      call mach_marine_halo(zo3, p2d, lfu, imod, ni_can, nkt)
      do k = 1, nkt
         do ic = 1, ni_can
            trppm(jj, ic, k) = zo3(ic, k)
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
      call mach_gas_messy(metvar3dcan, metvar2d, o3vmr2d, rjval, imod, &
                          ni_can, nkt)
!
! Skip photolysis cloud correction if MESSY is selected
      skyo = 1.0
!
   else
!    determine the cloudy/clear sky correction factor for photolysis rates:
!    cosine of solar zenith angle is positive; sun is above the horizon.
!    calculate photolysis correction factor for clouds.
      do k = 1, nkt
         do ic = 1, ni_can
            cldf(ic, k) = metvar3dcan(ic, k, MV3D_FTOT)  ! Cloud fraction
             ! Cloud liquid water content in g/m3
            cldw(ic, k) = metvar3dcan(ic, k, MV3D_QCPLUS) * 1000.0 * rho(ic, k)
         end do
      end do
      call mach_adom2_jcorr(zmom, cldw, cldf, csza, skyo, ni_can, nkt)
!
!     Evaluate photolysis rates from the (ADOM-II) look-up table
      call mach_gas_lookup_jvals(zplus, skyo, csza, rjval_lookup, ni_can, nkt, &
                                 ni_can*nkt)
!
      ! Placeholder for non-MESSY photolysis, useful for non-canopy case as well.
      rjval = 1.0d0

      ! Use the lookup table JVAL rates for ADOM2KPP is MESSy_Jval is disabled
      if (trim(chm_pkg_gas_s) == 'ADOM2KPP') rjval = rjval_lookup
   end if MESSY_photolysis
!*
! Forest canopy shading:
   corr = 1.0
   do ic = 1, ni_can
      ii = imod(ic) - 1
! Fraction of cumulative LAI above each canopy layer interface, counting down from top:
      laicumfrac(1) = busper(sm(sp_CRL) % per_offset + ii)
      laicumfrac(2) = busper(sm(sp_CRL + 1) % per_offset + ii)
      laicumfrac(3) = busper(sm(sp_CRL + 2) % per_offset + ii)
      laicumfrac(4) = busper(sm(sp_CRL + 3) % per_offset + ii)
!
!  Use BELD3 data for clumping index instead of the 15 land use categories:
      klight = 0.5 * busper(sm(sp_CLU) % per_offset + ii)
!
! The cumulative LAI is zero for the region above canopy layer 1, and is the
! full value for the region below canopy layer 2 (that is, the canopy is
! assumed to exist between those levels).
!
! Calculate the attenuation under the canopy (corr_f) and at various fractions
! of hc under the canopy.
!
      coeff = klight * busper(sm(sp_LAI) % per_offset + ii) / &
              max(0.05, csza(ic))
      corr_bio(1) = 1.0 !  top of canopy; z = hc
      corr_bio(2) = max(1.0E-10, exp(- laicumfrac(1) * coeff)) ! value at z = 0.75 hc
      corr_bio(3) = max(1.0E-10, exp(- laicumfrac(2) * coeff)) ! value at z = 0.50 hc
      corr_bio(4) = max(1.0E-10, exp(- laicumfrac(3) * coeff)) ! value at z = 0.35 hc
      corr_bio(5) = max(1.0E-10, exp(- laicumfrac(4) * coeff)) ! value at z = 0.20 hc
      corr_bio(6) = max(1.0E-10, exp(-coeff))
      zbio (1) = busper(sm(sp_HC) % per_offset + ii)
      zbio (2) = zbio(1) * 0.75
      zbio (3) = zbio(1) * 0.50
      zbio (4) = zbio(1) * 0.35
      zbio (5) = zbio(1) * 0.20
      zbio (6) = 0.0
! work out correction factors for all momentum levels:
      do k = 1, nkt
         if (zmom(ic, k) >= zbio(1)) then
            corrmom(k) = 1.0   ! no attenuation; canopy model layer is above hc
         else
            do kk = 1, 5
               if (zmom(ic, k) < zbio(kk) .and. zmom(ic, k) >= zbio(kk+1)) then
                  corrmom(k) = corr_bio(kk+1) + (corr_bio(kk) - corr_bio(kk+1)) / &
                               (zbio(kk) - zbio(kk+1)) * (zmom(ic, k) - zbio(kk+1))
               end if
            end do
         end if
      end do
      corrmom(nkt + 1) = corr_bio(6)
! thermodynamic level correction:  average of momentum level corrections:
      do k = 1, nkt
         corr(ic, k) = (corrmom(k) + corrmom(k+1)) * 0.5
!
!  Scale the model rate constants downward due to shading:
         ik0 = ik(ic+1, k, ni_can)
         do jj = 1, njrxs
            rjval(ik0, jj) = rjval(ik0, jj) * corr(ic, k)
         end do
      end do
   end do
!**
!
!  Initialize the total organic aerosols
   moi = 0.0d0
   if (chm_soa_s /= 'NIL') then
      soa = 0.0d0
      poa = 0.0d0
      do ic = 1, ni_can
         ii = imod(ic)
!
! (i): Model resolved layers:
         do kk = 1, chm_nk
            k = kmod(ic, kk)
            do jj = 1, isize
               isp_OC = (iae_OC - 1) * isize + jj + sp_AERO - 1
               soa(ic, k) = soa(ic, k) + dble(chem_tr(ii, kk, isp_OC))
               isp_PC = (iae_PC - 1) * isize + jj + sp_AERO - 1
               poa(ic, k) = poa(ic, k) + dble(chem_tr(ii, kk, isp_PC))
               ! add total condensible organics, convert from ug/kg to ug/m3
               moi(ic, k) = (soa(ic, k) + poa(ic, k)) * dble(rho(ic, k))
            end do
         end do
!
! (ii): Canopy shaded layers:
         do kc = 1, nkc
            k = kcan(ic, kc)
            do jj = 1, isize
               isp_OC = (iae_OC - 1) * isize + jj + sp_AERO - 1
               soa(ic, k) = soa(ic, k) + dble(tracers_can(ic, kc, isp_OC))
               isp_PC = (iae_PC - 1) * isize + jj + sp_AERO - 1
               poa(ic, k) = poa(ic, k) + dble(tracers_can(ic, kc, isp_PC))
               moi(ic, k) = (soa(ic, k) + poa(ic, k)) * dble(rho(ic, k))
            end do
         end do
      end do
   end if
!
! KPP solver internal time-step to the bus
   if (var_step > 0) then
      do ic = 1, ni_can
         ii = imod(ic)
! (i): Model resolved layers:
         do kk = 1, chm_nk
            k = kmod(ic, kk)
            this_ik = ik(ii, kk, chm_ni)
            hstart(ic, k) = dble(busper(sm(var_step) % per_offset + this_ik))
         end do
! (ii): Canopy shaded layers:
         do kc = 1, nkc
            k = kcan(ic, kc)
            this_ik = ik(ii, kc, chm_ni)
            hstart(ic, k) = dble(busper(sm(var_step) % dd_offset + this_ik))
         end do
      end do
   end if
!
   select case (trim(chm_pkg_gas_s))
      case ('ADOM2')
      !  ADOM-II mechanism with the Young & Boris Solver
         call mach_adom2yb_main(p2d, tplus, hu_ppm, sigt, rjval, trppm, &
                                rjval_lookup, csza, voc_diff, step, ni_can, nkt)

      !  ADOM-II mechanism with the KPP Solver
      case ('ADOM2KPP')
         call mach_adom2kpp_main(p2d, tplus, hu_ppm, sigt, rjval, trppm, &
                                 hstart, voc_diff, ni_can, nkt)
!
      !  ADOM-II mechanism with the KPP Solver, using rates from adom2_uprate
      case ('ADOM2KPPB')
         call mach_adom2kppb_main(p2d, tplus, hu_ppm, sigt, rjval, trppm, &
                                  hstart, rjval_lookup, voc_diff, ni_can, nkt)
!
      !  SAPRC07C mechanisms with the KPP Solver
      case ('SAPRC07C', 'SAPRC07CS')
         call mach_saprc_main(p2d, tplus, hu_ppm, rjval, trppm, hstart, &
                              moi, dsoa, ni_can, nkt)
   end select

   if (chm_error_l) return
!
!============================================================================
! Return updated gas concentrations into the dynamic bus
! (Convert units back from ppmv to mass mixing ratio)
!============================================================================
   do ic = 1, ni_can
      ii = imod(ic)
!
! (i): Model resolved layers:
      do kk = nk_start, chm_nk
         k = kmod(ic, kk)
         do isp = 1, nvar
            jj = gas_species(isp)
            chem_tr(ii, kk, jj) = trppm(isp, ic, k) / conv1(ic, k) * sm(jj) % mol_wt
         end do
         if (chm_active_ch4_l) then
            chem_tr(ii, kk, sp_CH4) = trppm(ind_ch4, ic, k) / &
                                      conv1(ic, k) * sm(sp_CH4) % mol_wt
         end if
      end do
! (ii): Canopy shaded layers:
      do kc = 1, nkc
         k = kcan(ic, kc)
         do isp = 1, nvar
            jj = gas_species(isp)
            tracers_can(ic, kc, jj) = trppm(isp, ic, k) / conv1(ic, k) * &
                                      sm(jj) % mol_wt
         end do
         if (chm_active_ch4_l) then
            tracers_can(ic, kc, sp_CH4) = trppm(ind_ch4, ic, k) / &
                                          conv1(ic, k) * sm(sp_CH4) % mol_wt
         end if
      end do
!
   end do
!
! Return the KPP solver internal time-step to the bus
   if (var_step > 0) then
      do ic = 1, ni_can
         ii = imod(ic)
! (i): Model resolved layers:
         do kk = 1, chm_nk
            k = kmod(ic, kk)
            this_ik = ik(ii, kk, chm_ni)
            busper(sm(var_step) % per_offset + this_ik) = real(hstart(ic, k))
         end do
! (ii): Canopy shaded layers:
         do kc = 1, nkc
            k = kcan(ic, kc)
            this_ik = ik(ii, kc, chm_ni)
            busper(sm(var_step) % dd_offset + this_ik) = real(hstart(ic, kc))
         end do
      end do
   end if
!
 !---- write Jx value to output ---------
   if (sp_JNO2 > 0) then
      ! Output JNO2 at the lowest model layer
      do ic = 1, ni_can
         do kk = 1, chm_nk
            ik0 = ik(ic+1, kk, ni_can)
            busvol(sm(sp_JNO2) % out_offset + imod(ic) - 1) = real(rjval(ik0, 1))
         end do
      end do
   end if
!
!!!  need to set "chm_debug_3d_i = 6" in gem_settings.nml and
!!!  add "sortie_p([3DB1,3DB2,3DB3,3DB4,3DB5,3DB6])" in outcfg.out
!   if (chm_debug_3d_i >= 6) then
!      do ic = 1, ni_can
!         ii = imod(ic)
!         do kk = 1, chm_nk
!            k = kmod(ic, kk)
!            this_ik = ik(ii, kk, chm_ni)
!            ik0 = ik(ic+1, k, ni_can)
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
         call mach_soa_jiang(moi, dsoa, voc_diff, tplus, p2d, ni_can, nkt)
      end if
!
      do ic = 1, ni_can
         ii = imod(ic)
!
! (i): Model resolved layers:
         do kk = nk_start, chm_nk
            k = kmod(ic, kk)
            this_ik = ik(ii, kk, chm_ni)
        ! convert dsoa from ug/m3 to ug/kg and output to vol bus
            busvol(sm(sp_NWOC) % out_offset + this_ik) = dsoa(ic, k) / rho(ic, k)
         end do
! (ii): Canopy shaded layers:
! Note the use of sm(sp_NWOC) % vd_offset here, to hold the SOA yields within canopy
         do kc = 1, nkc
            k = kcan(ic, kc)
            this_ik = ik(ii, kc, chm_ni)
            busvol(sm(sp_NWOC) % vd_offset + this_ik) = dsoa(ic, k) / rho(ic, k)
         end do
!
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
         lino3_new(:, :) = metvar3dcan(:, :, MV3D_O3L)
      else
!      Evaluate LINOZ O3 tendency locally
         call mach_gas_strato(busper, hu_ppm, tplus, sigt, o3_old, psurf, &
                              lino3_new, ni_can, nkt, imod)
      end if
!
!     Apply LINOZ tendencies to ozone (O3) above hu_linoz_tropo specific humidity or
!     for pressure lower than p_linoz_tropo,
      do kk = 1, chm_nk
         do ic = 1, ni_can
            k = kmod(ic, kk)
            if ((hu_ppm(ic, k) < hu_linoz_tropo) .or. &
                (p2d(ic, k) < p_linoz_tropo)) then
               ii = imod(ic)
               chem_tr(ii, kk, sp_O3) = lino3_new(ic, k)
            end if
         end do
      end do
!
!    Output diagnostics
!
      IF_DIAG: if (chm_diag_colum_L) then
         do kk = 1, chm_nk
            do ic = 1, ni_can
               k = kmod(ic, kk)
               ii = imod(ic)
               this_ik = ik(ii, kk, chm_ni)
         ! Ozone tendency
               o3_tend = (chem_tr(ii, kk, sp_O3) - o3_old(ic, k)) / chm_timestep    !ug /kg /sec
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
end subroutine mach_gas_canopy
