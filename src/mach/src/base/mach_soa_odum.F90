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
! Fichier/File   : mach_soa_odum.ftn90
! Creation       : J. Chen, C. Stroud, and D. Akingunola, Oct. 2020
! Description    : Calculates SOA formed for input dVOC from precursor
!                  gases. SOA based on 5-bin volativility of Donahue et al.
!                  accounting for 6 oxidation pathways and NOx conditions
! rewrite of the Odum equation by Neil Donahue (5 term volatility set)
! Odum equation has 2 terms in summation each with an alpha and K value
!  so there are 4 fit variables.
! The Donahue equation has 5 terms in summation with each term having
!  an alpha and Ca value. Ca values are fixed in the fitting to chamber data
!  only the alpha variables that change and are optimized in fit to chamber data.
! Note that K = 1/C* in comparing Odum and Donahue equations,
!(K is gas/particle equilibrium coefficient and C* is saturation concentration).
!
! Yield(i)  = Alpha(1,i)*(1.0d0/(1.0d0+((Ca(1,i))/Mo))) +
!             Alpha(2,i)*(1.0d0/(1.0d0+((Ca(2,i))/Mo))) +
!             Alpha(3,i)*(1.0d0/(1.0d0+((Ca(3,i))/Mo))) +
!             Alpha(4,i)*(1.0d0/(1.0d0+((Ca(4,i))/Mo))) +
!             Alpha(5,i)*(1.0d0/(1.0d0+((Ca(5,i))/Mo)))
! The alpha values are mass stoichiometric coefficients for products of a
!  precursor reaction and the C values are saturation vapour concentrations
!  (binned into C values, one order of magnitude apart).
! It is not the VBS method. Here we calculate the mass of SOA formed in each
!  chemistry time step and it is assumed to be irreversibly formed.
! The true VBS method assumes real equilibrium for all the organic aerosol
!  mass at each time step and the mass can evaporate; also assumes that
!  gas-phase organic mass in each volatility bin can chemically age:
! e.g. mass moves from bin 1->2 in a SVOC+OH reaction, bin 2->3...
!
! Input dvoc is in ppmv and is convert to ug/m3 when converted to
!  soa precurosr species in dvoc_soa
!
! Modified       :
!
! Arguments:
!       IN
!                conx1      -> ppmv to molec/m3 conversion factor
!                tplus      -> temperature (k)
!                ksoa       -> reaction rates
!                kNO_r2o2   -> reaction rates
!                kHO2_r2o2  -> reaction rates
!                dvoc       -> soa precursor gas increment (ppmv)
!                moi        -> total condensible organics (ug/m3)
!                trppm      -> species conc. array (ppmv)
!
!       OUT
!                dsoa       -> dSOA formed from input dVOC (ug/m3)
!
!============================================================================
!
!!if_on
subroutine mach_soa_odum(moi, dsoa, dvoc, ksoa, kNO_r2o2, kHO2_r2o2, tplus, &
                         conx1, trppm, gni, gnk)
   use chm_utils_mod,      only: dp, sp, i4
   use mach_pkg_gas_mod,   only: nsp_soa_gases, nspec, nvsoa, minconc
!!if_off
   use chm_utils_mod,      only: CHM_MSG_DEBUG
   use chm_nml_mod,        only: chm_pkg_gas_s, chm_timings_l, nk_start
   use chm_consphychm_mod, only: avno, rgasi ! gas constant (J K-1 mol-1)
   use mach_pkg_gas_mod,   only: ind_o3, ind_no3, ind_oh, ind_no, ind_ho2
   use mach_soa_odum_mod,  only: alpha, cmatx, h_vap, ncond, nvbin, mstep,   &
                                 lo_OH, hi_OH, lo_O3, hi_O3, lo_NO3, hi_NO3, &
                                 i_kNO3, i_kO3, i_kOH, map_soa_dvoc_svoc_saprc

   implicit none
!!if_on
   integer(i4), intent (in) :: gni, gnk
   real(dp),    intent (in) :: dvoc (nsp_soa_gases, gni, gnk)
   real(dp),    intent (in) :: moi  (gni, gnk)
   real(dp),    intent (in) :: ksoa (3, nvsoa, gni, gnk)
   real(dp),    intent (in) :: kHO2_r2o2(gni, gnk)
   real(dp),    intent (in) :: kNO_r2o2 (gni, gnk)
   real(sp),    intent (in) :: tplus(gni, gnk)
   real(dp),    intent (in) :: conx1(gni * gnk)
   real(sp),    intent (in) :: trppm(nspec, gni, gnk)
   real(sp),    intent(out) :: dsoa (gni, gnk)
!!if_off
!
!  Declaration of local variables
   real(dp), dimension(nvbin, ncond, nvsoa) :: cstar
   real(dp), dimension(ncond, nvsoa)        :: yield
   real(dp), dimension(nvsoa, gni, gnk)     :: dvoc_soa
   real(dp), dimension(nvsoa)               :: OHfrac, O3frac, NO3frac,     &
                                               pvoc, psoa, soa_yield, tdep, &
                                               hi_nox_yield, lo_nox_yield
   real(dp) :: cno3, cno, cho2, coh, co3
   real(dp) :: total_dVOC_loss, hi_NOx_frac, tt1, tt2, mop
   integer(i4) :: nn, k, i, step, spc, cond
!
!  Declaration of external subroutines
   external msg_toall, timing_start_omp, timing_stop_omp

!===============================================================================
   call msg_toall(CHM_MSG_DEBUG, 'mach_soa_odum [BEGIN]')
   if (chm_timings_L) call timing_start_omp(336, 'mach_soa_odum', 330)

!  Initializations
   dsoa = 0.0

   if (chm_pkg_gas_s(1:5) == 'SAPRC') then
      ! input dvoc in pppmv, return dvoc_soa in ug/m3
      call map_soa_dvoc_svoc_saprc( dvoc, dvoc_soa, tplus, gni, gnk )
   end if

   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
      do i = 1, gni
         if (moi(i, k) <= minconc) cycle
! prepare species conc. convert ppmv to molec/cm3 for rate calculation
         nn = nn + 1
         co3  = dble(trppm(ind_o3, i, k)) * conx1(nn)
         cno3 = dble(trppm(ind_no3, i, k)) * conx1(nn)
         coh  = dble(trppm(ind_oh, i, k)) * conx1(nn)
         cno  = dble(trppm(ind_no, i, k)) * conx1(nn)
         cho2 = dble(trppm(ind_ho2, i, k)) * conx1(nn)

         tt1 = 298.0d0 / dble(tplus(i, k))
         tt2 = (1.d0 / 298.0d0) - (1.0d0 / dble(tplus(i, k)))
!! Calculates oxidation fraction by select pathways (OH, O3, NO3) and estimate
!! the fraction under high/low NOx conditions-- based on R2O2+NO2 or R2O2+HO2 Rx.
         hi_NOx_frac = kNO_r2o2(i, k) * cno / &
                       max(kHO2_r2o2(i, k) * cho2 + kNO_r2o2(i, k) * cno, minconc)
         do spc = 1, nvsoa
            total_dVOC_loss = max(ksoa(i_kOH, spc, i, k) * coh + &
                                  ksoa(i_kO3, spc, i, k) * co3 + &
                                  ksoa(i_kNO3,spc, i, k) * cno3, minconc)
            OHfrac(spc)  = ksoa(i_kOH, spc, i, k) * coh / total_dVOC_loss
            O3frac(spc)  = ksoa(i_kO3, spc, i, k) * co3 / total_dVOC_loss
            NO3frac(spc) = ksoa(i_kNO3, spc, i, k) * cno3 / total_dVOC_loss

            pvoc(spc) = dvoc_soa(spc, i, k) / dble(mstep)
            tdep(spc) = tt1 * exp((h_vap / dble(rgasi)) * tt2)
         end do
         mop = moi(i, k)
         do step = 1, mstep
            do spc = 1, nvsoa
               ! temp dependent C* correction
               cstar(:, :, spc) = cmatx(:, :, spc) * tdep(spc)
               do cond = 1, ncond
                  yield(cond, spc) = alpha(1, cond, spc) * &
                          (1.0d0 / (1.0d0 + ((cstar(1, cond, spc)) / mop))) + &
                                     alpha(2, cond, spc) * &
                          (1.0d0 / (1.0d0 + ((cstar(2, cond, spc)) / mop))) + &
                                     alpha(3, cond, spc) * &
                          (1.0d0 / (1.0d0 + ((cstar(3, cond, spc)) / mop))) + &
                                     alpha(4, cond, spc) * &
                          (1.0d0 / (1.0d0 + ((cstar(4, cond, spc)) / mop))) + &
                                     alpha(5, cond, spc) * &
                          (1.0d0 / (1.0d0 + ((cstar(5, cond, spc)) / mop)))
               end do
! calculate yield for low and high NOx limit and for each oxidant proportion
               hi_nox_yield(spc) = yield(hi_OH, spc) * OHfrac(spc) + &
                                   yield(hi_O3, spc) * O3frac(spc) + &
                                   yield(hi_NO3, spc) * NO3frac(spc)
               lo_nox_yield(spc) = yield(lo_OH, spc) * OHfrac(spc) + &
                                   yield(lo_O3, spc) * O3frac(spc) + &
                                   yield(lo_NO3, spc) * NO3frac(spc)
! final yield of the aerosol from the specific precursor for portion of
!  yield under high NOx and portion under low NOx.
               soa_yield(spc) = hi_nox_yield(spc) * hi_NOx_frac + &
                                lo_nox_yield(spc) * (1.0d0 - hi_NOx_frac)
! calculate the SOA mass contributions from each precursor [unit ug/m3]
               psoa(spc) = soa_yield(spc) * pvoc(spc)
            end do
            mop = mop + sum(psoa(:)) ! increment total condensible organics in mstep
         end do ! end Mop msteps
         ! ensure positivity and return
         dsoa(i, k) = real(max(mop - moi(i, k), 0.0d0))
      end do
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_soa_odum [END]')
   if (chm_timings_L) call timing_stop_omp(336)

   return
end subroutine mach_soa_odum
