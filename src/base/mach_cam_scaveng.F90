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
! Fichier/File   : mach_cam_scaveng.ftn90
! Creation       : S. Gong, W. Gong, A. DASTOOR, C. Stroud, S. Gravel and
!                  B. Pabla for GEM-MACH, June 2008
! Description    : Calculates below cloud scavenge processes of gas and particle
!
! Extra info     : - First version created by S. Gong Jun 05 1994 for CAM
!                  - Vectorized the whole program and add working spaces.
!                    (S. Gong, Jan 19, 1996)
!                  - Separate from incloud processes to deal the below-cloud
!                    scavenging name changed from wetrem to scaveng.
!                    (S. Gong, Dec 20, 1998)
!                  - Added averaging for fractional cloud cover and
!                    evaporatrion modification. also changed the objective
!                    to in and below cloud rain scavenging.
!                    (A. Dastoor, Apr 25, 2001)
!                  - 1. added the subroutine for calculating CC2D
!                       (code supplied by Ashu); and
!                    2. added scavenging of soluble gases and the
!                       accumulation of surface wet fluxes (additional to
!                       the fluxes due to cloud-to-rain conversion).
!                      (W. Gong, May 2001)
!                 - Split aerosol OC into primary and secondary
!                   components !cs>>> (C. Stroud, Aug 2004)
!
! Arguments:  IN
!
!                throw   -> Temperature
!                cf      -> 3d cloud fraction
!                evpfac  -> evap. of strat. precip (consun)
!                pres    -> Local mid-layer pressure [Pa]
!                pdepv   -> gravitational settling velocity
!                qr      -> Liquid and solid precipitation
!                rhsize  -> Unactivated ambient aerosol wet
!                rhop    -> final wet density
!                pdiff   -> diffusion parameter
!                ntr     -> Total number of trace substances (gases and aerosols)
!                thlev   -> Layer thickness [m]
!                roarow  -> Air density (kg/m3)
!                amu     -> Air's dynamic viscosity
!                amfp    -> Mean molecular free path
!
!             IN/OUT
!
!                XROW    -> tracers concentration
!                WETFLX  -> Wet flux
!                FLUX    -> H+ in rain fluxes (including cloud-to-rain conversion)
!
!             OUT
!
!                RTBCLD  -> Rain Scavenging
!
!             LOCAL
!
!                GDREM   -> total removal tendency
!============================================================================
!
!!if_on
subroutine mach_cam_scaveng(throw, xrow, cf, evpfac, pdepv, qr, rtbcld, &
                            rhsize, rhop, pdiff, thlev, roarow, wetflx, &
                            flux, pres, amu, amfp, pni, pnk)
   use mach_cam_utils_mod,      only: isize, nswdep, ntr
!!if_off
   use chm_nml_mod,             only: chm_timings_L, chm_split_snow_rain_l, chm_wdep_scav_coef_s
   use chm_utils_mod,           only: chm_timestep, CHM_MSG_DEBUG
   use chm_consphychm_mod,      only: mwt_air
   use mach_aurams_headers_mod, only: mach_aurams_cldcv2d
   use mach_cam_headers_mod,    only: mach_cam_rain
   use mach_cam_utils_mod,      only: iae1, iae2, icom, tmin
   use mach_cam_utils_mod,      only: igs_SO2, igs_H2O2, igs_ROOH, igs_HNO3,  &
                                      igs_NH3, mwt_igs
   use mach_cam_utils_mod,      only: iae_SU, iae_NI, iae_AM, mwt_aero, ip_wflx
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: throw (pni, pnk)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: cf    (pni, pnk)
   real(kind=4),    intent   (in) :: evpfac(pni, pnk)
   real(kind=4),    intent   (in) :: pres  (pni, pnk)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: qr    (pni, pnk, 2)
   real(kind=4),    intent  (out) :: rtbcld(pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: amfp  (pni, pnk)
   real(kind=4),    intent(inout) :: wetflx(pni, nswdep)
   real(kind=4),    intent(inout) :: flux  (pni, pnk, nswdep)
!!if_off
!
! Local variables
!
   integer(kind=4)  :: nt
   integer(kind=4)  :: l, i, ii, n, np1, np2, np3, np
   real(kind=4)     :: b2, ti, xold, rph, xnew, tend, qr_flux
   real(kind=4)     :: rhso3(pni), rhco3(pni)
   real(kind=4)     :: cc2d(pni, pnk), inv_cc2d(pni, pnk)
   real(kind=4)     :: mwt_hso3, mwt_hco3
   real(kind=4)     :: hso2, hco2, phco3, phso3
   real(kind=4)     :: pxnew(pni, pnk, 2), qr_vel(pni, pnk, 2)
   real(kind=4)     :: gdrem(pni, ntr), gdrem_co2(pni)
   real(kind=4)     :: rtbcld_co2(pni, pnk)
   real(kind=4)     :: mass_per_unit_area(pni, pnk)
   real(kind=4), parameter :: smf1 = 1.0e-15, smf = 1.0e-04
!
!  External subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_cam_scaveng [BEGIN]')
   if (chm_timings_L) call timing_start_omp(375, 'mach_cam_scaveng', 340)

   gdrem     = 0.0
   gdrem_co2 = 0.0

!  get gas molecular weight

   mwt_hso3 = 81.067  !sm(??) % mol_wt   !GEM-MACH doesn't include this yet
   mwt_hco3 = 61.017  !sm() % mol_wt     !GEM-MACH doesn't include this yet

!  calculate CC2D ("2D" cloud cover at a given level)
   call mach_aurams_cldcv2d(cc2d, cf, pni, pnk)

   if (chm_split_snow_rain_l .or. (chm_wdep_scav_coef_s == 'WANG2014')) then
      do n = 1, 2
         do l = 1, pnk
            do i = 1, pni
               inv_cc2d(i, l) = 1.0 / max(cc2d(i, l), smf)
               pxnew(i, l, n) = max(qr(i, l, n), 0.0) * inv_cc2d(i, l)
            end do
         end do
      end do
   else
!  sum up precipitation [stratiform only]
      do l = 1, pnk
         do i = 1, pni
            inv_cc2d(i, l) = 1.0 / max(cc2d(i, l), smf)
            pxnew(i, l, 1) = max(qr(i, l, 1) + qr(i, l, 2), 0.0) * inv_cc2d(i, l)
            pxnew(i, l, 2) = 0.0
         end do
      end do
   end if

   call mach_cam_rain(throw, pxnew, rhsize, pdepv, pres, roarow, pdiff, amu, &
                      amfp, xrow, rtbcld, rhop, cc2d, qr_vel, pni, pnk)

!  qr_vel contains average rain fall velocity field.
!  evaluate H+ in rain fluxes (including cloud-to-rain conversion)
   do i = 1, pni
      rhco3(i) = 0.0
      rhso3(i) = 0.0
   end do
   do l = 1, pnk
      do i = 1, pni
         flux(i, l, 9) = flux(i, l, 9) + 1000.0 * roarow(i, l) * chm_timestep * &
                         (-1.0 * rtbcld(i, l, igs_HNO3) / mwt_igs(igs_hno3)   + &
                         rtbcld(i, l, igs_NH3) / mwt_igs(igs_nh3)) * inv_cc2d(i, l)
      end do
   end do
   do n = 1, isize
      np1 = isize * (iae_SU - 1) + n + (iae1 - 1)
      np2 = isize * (iae_NI - 1) + n + (iae1 - 1)
      np3 = isize * (iae_AM - 1) + n + (iae1 - 1)
      do l = 1, pnk
         do i = 1, pni
            flux(i, l, 9) = flux(i, l, 9) + 1000.0 * roarow(i, l) * chm_timestep *   &
                            (-2.0 * rtbcld(i, l, np1) / mwt_aero(iae_SU)   &
                            - rtbcld(i, l, np2) / mwt_aero(iae_NI)         &
                            + rtbcld(i, l, np3) / mwt_aero(iae_AM)) * inv_cc2d(i, l)
         end do
      end do
   end do

!  Equilibrium scavenging of SO2 and CO2:
   do l = 1, pnk
      do i = 1, pni
!
         mass_per_unit_area(i, l) = thlev(i, l) * roarow(i, l)

!  equivalent Henry's Law coefficients (incorporating dissociation constant)
!  (Same as in mach_incld_upaqr1) for SO2 and CO2
         ti = 1.0 / throw(i, l)
         hso2 = 10.0 ** (-10.65 + 1410.0 * ti) * 10.0 ** (-4.84 + 870.0 * ti)
         hco2 = 10.0 ** (-10.66 + 760.0 * ti + 58000.0 * ti * ti) * 10.0 **   &
                (-14.25 + 5190.0 * ti - 850000.0 * ti * ti)
!  unit conversion factor (MOLAR to ppm)
!  (mwt_air/rho_air/rho_water*RWC*1.E+06, where rho_air in kg/m3 and
!  RWC in kg_water/m3_air and rho_water=1000. (kg/m3))
!  note: assuming QR in kg m-2 s-1
         b2 = mwt_air * pxnew(i, l, 1) / (qr_vel(i, l, 1) + smf) / roarow(i, l) * 1000.0
         qr_flux = qr_vel(i, l, 1) / (pxnew(i, l, 1) + smf1)

!  Equilibrium scavenging of CO2 and SO2:
!  CO2
         xold = 340.0e-06 * 44.01 / mwt_air
         phco3 = rhco3(i) / mass_per_unit_area(i, l) * inv_cc2d(i, l)
         rph = flux(i, l, 9) * qr_flux
         rph = max(2.5e-06, rph)
         xnew = xold - (b2 * hco2 * xold - phco3 * rph) / (rph + b2 * hco2)
         tend = min((xnew - xold) / chm_timestep, 0.0)
         rtbcld_co2(i, l) = tend * cc2d(i, l)
         rhco3(i) = (rhco3(i) - rtbcld_co2(i, l) * chm_timestep         &
                    * mass_per_unit_area(i, l)) * (1.0 - evpfac(i, l))
! SO2
         xold = xrow(i, l, igs_SO2)
         phso3 = rhso3(i) / mass_per_unit_area(i, l) * inv_cc2d(i, l)
! update H+ after scavenging of CO2:
         rph = (flux(i, l, 9) + 1000.0 * roarow(i, l) * chm_timestep * &
                rtbcld_co2(i, l) / mwt_hco3 * inv_cc2d(i, l)) * qr_flux
! temporarily skip CO2 scavenging
         rph = max(2.5e-06, rph)
         xnew = xold - (b2 * hso2 * xold - phso3 * rph) / (rph + b2 * hso2)
         tend = min((xnew - xold) / chm_timestep, 0.0)
         xrow(i, l, igs_SO2) = max(xold + tend * chm_timestep, tmin)          &
                               * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
         rtbcld(i, l, igs_SO2) = tend * cc2d(i, l)
         rhso3(i) = (rhso3(i) - rtbcld(i, l, igs_SO2) * chm_timestep          &
                    * mass_per_unit_area(i, l)) * (1.0 - evpfac(i, l))
      end do
   end do

!  add ground removal (aerosols)
   do np = iae1, iae2
      do l = 1, pnk
         do i = 1, pni
            if (xrow(i, l, np) > tmin) then
!  note gdrem and rtbcld are negative
               gdrem(i, np) = gdrem(i, np) + rtbcld(i, l, np) * mass_per_unit_area(i, l)
               xrow(i, l, np) = xrow(i, l, np) - gdrem(i, np) * evpfac(i, l) &
                                * chm_timestep / mass_per_unit_area(i, l)
               gdrem(i, np) = gdrem(i, np) * (1.0 - evpfac(i, l))
            end if
         end do
      end do
   end do

!  To account for scavenging of soluble gases (Note that the evaporated nitrate and ammonia,
!  originally scavenged from gas phase HNO3 and NH3 are put back into gas phase. It could
!  be argued that they could all, or partly, be in aerosol phase. Also note that gas phase CO2
!  is not updated since its high base value being held at 340 ppm.)
   do l = 1, pnk
      do i = 1, pni
! SO2/HSO3
         if (xrow(i, l, igs_SO2) > tmin) then
            gdrem(i, igs_SO2) = gdrem(i, igs_SO2) + rtbcld(i, l, igs_SO2) * mass_per_unit_area(i, l)
            xrow(i, l, igs_SO2) = xrow(i, l, igs_SO2) - gdrem(i, igs_SO2) * evpfac(i, l) &
                                  * chm_timestep / mass_per_unit_area(i, l) * mwt_igs(igs_so2) / mwt_hso3
            gdrem(i, igs_SO2) = gdrem(i, igs_SO2) * (1.0 - evpfac(i, l))
         end if
! H2O2
         if (xrow(i, l, igs_H2O2) > tmin) then
            gdrem(i, igs_H2O2) = gdrem(i, igs_H2O2) + rtbcld(i, l, igs_H2O2) * mass_per_unit_area(i, l)
            xrow(i, l, igs_H2O2) = xrow(i, l, igs_H2O2) - gdrem(i, igs_H2O2) * evpfac(i, l) &
                                  * chm_timestep / mass_per_unit_area(i, l)
            gdrem(i, igs_H2O2) = gdrem(i, igs_H2O2) * (1.0 - evpfac(i, l))
         end if
! ROOH
         if (xrow(i, l, igs_ROOH) > tmin) then
            gdrem(i, igs_ROOH) = gdrem(i, igs_ROOH) + rtbcld(i, l, igs_ROOH) * mass_per_unit_area(i, l)
            xrow(i, l, igs_ROOH) = xrow(i, l, igs_ROOH) - gdrem(i, igs_ROOH) * evpfac(i, l) &
                                  * chm_timestep / mass_per_unit_area(i, l)
            gdrem(i, igs_ROOH) = gdrem(i, igs_ROOH) * (1.0 - evpfac(i, l))
         end if
! HNO3
         if (xrow(i, l, igs_HNO3) > tmin) then
            gdrem(i, igs_HNO3) = gdrem(i, igs_HNO3) + rtbcld(i, l, igs_HNO3) * mass_per_unit_area(i, l)
            xrow(i, l, igs_HNO3) = xrow(i, l, igs_HNO3) - gdrem(i, igs_HNO3) * evpfac(i, l) &
                                  * chm_timestep / mass_per_unit_area(i, l)
            gdrem(i, igs_HNO3) = gdrem(i, igs_HNO3) * (1.0 - evpfac(i, l))
         end if
! NH3
         if (xrow(i, l, igs_NH3) > tmin) then
            gdrem(i, igs_NH3) = gdrem(i, igs_NH3) + rtbcld(i, l, igs_NH3) * mass_per_unit_area(i, l)
            xrow(i, l, igs_NH3) = xrow(i, l, igs_NH3) - gdrem(i, igs_NH3) * evpfac(i, l) &
                                 * chm_timestep / mass_per_unit_area(i, l)
            gdrem(i, igs_NH3) = gdrem(i, igs_NH3) * (1.0 - evpfac(i, l))
         end if
! HCO3
         gdrem_co2(i) = gdrem_co2(i) + rtbcld_co2(i, l) * mass_per_unit_area(i, l)
         gdrem_co2(i) = gdrem_co2(i) * (1.0 - evpfac(i, l))
      end do
   end do
!  Combine fluxes (surface) with the fluxes from cloud-to-rain; (moles m-2 for chm_timestep, positive)
   do nt = 1, icom
      ii = ip_wflx(nt)
      do n = 1, isize
         np = isize * (nt - 1) + n + (iae1 - 1)
         do i = 1, pni
            wetflx(i, ii) = wetflx(i, ii) - gdrem(i, np) * chm_timestep * 1000.0 / mwt_aero(nt)
         end do
      end do
   end do
   do i = 1, pni
      wetflx(i, 1) = wetflx(i, 1) - gdrem(i, igs_SO2)  * chm_timestep * 1000.0 / mwt_hso3
      wetflx(i, 2) = wetflx(i, 2) - gdrem(i, igs_H2O2) * chm_timestep * 1000.0 / mwt_igs(igs_h2o2)
      wetflx(i, 3) = wetflx(i, 3) - gdrem(i, igs_ROOH) * chm_timestep * 1000.0 / mwt_igs(igs_rooh)
      wetflx(i, 5) = wetflx(i, 5) - gdrem(i, igs_HNO3) * chm_timestep * 1000.0 / mwt_igs(igs_hno3)
      wetflx(i, 6) = wetflx(i, 6) - gdrem(i, igs_NH3)  * chm_timestep * 1000.0 / mwt_igs(igs_nh3)
      wetflx(i, 8) = wetflx(i, 8) - gdrem_co2(i) * chm_timestep * 1000.0 / mwt_hco3
!  wet flux of h+ from charge balence:
      wetflx(i, 9) = wetflx(i, 1) + 2.0 * wetflx(i, 4) + wetflx(i, 5) &
                   - wetflx(i, 6) + wetflx(i, 8) - wetflx(i, 7)
      wetflx(i, 9) = max(wetflx(i, 9), 0.0)
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_cam_scaveng [END]')
   if (chm_timings_L) call timing_stop_omp(375)
   !-----------------------------------------------------------------

   return
end subroutine mach_cam_scaveng
