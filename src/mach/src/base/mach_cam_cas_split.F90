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
! Fichier/File   : mach_cam_cas.ftn90
! Creation       : S. Gong, W. Gong, S. Gravel and B. Pabla for GEM-MACH, June 2008
! Description    : This is a module of calculating the collection efficiency of
!                  the aerosols of type n by collector droplets of radius rcol.
! Extra info     : - First version created by S. Gong Jul 08 1994 for CAM
!                  - Vectorized the whole program and add some working spaces.
!                    (S. Gong, Dec 19, 1996)
!                  - Modified terminal/settling velocity for drops (> 20 um) and
!                    added gas scavenging rate calculations for H2O2, ROOH, HNO3
!                    and NH3. (W. Gong, May 2001)
!
! Arguments:  IN
!               tp     -> Temperature (K)
!               rhop   -> Final wet density
!               roarow -> Air density (kg/m3)
!               rhsize -> Wet radius
!               pres   -> Local mid-layer pressure [Pa]
!               qr     -> rain water/snow
!               pdiff  -> diffusion coefficient
!               pdepv  -> gravitational settling velocity
!
!            OUT
!               qr_vel -> average settling velocity
!               colef  -> Rain scavenging efficiency
!               rscavg -> Rain scavenging rate (by falling rain drops)
!                         for h2o2, rooh, hno3, nh3
!
!============================================================================
!
!!if_on
subroutine mach_cam_cas_split(tp, colef, rscavg, rhop, roarow, rhsize, pres, &
                              qr, qr_vel, pdiff, pdepv, amu, amfp, pni, pnk)
   use mach_cam_utils_mod,   only: isize
!!if_off
   use chm_nml_mod,          only: chm_wdep_scav_coef_s
   use chm_consphychm_mod, only: rgasi, grav, mwt_air, pi, t1s
   use mach_cam_utils_mod, only: mwt_igs, igs_H2O2, igs_HNO3, igs_NH3, igs_ROOH
   implicit none

!!if_on
   integer(kind=4), intent (in) :: pni, pnk
   real(kind=4),    intent (in) :: tp    (pni, pnk)
   real(kind=4),    intent(out) :: colef (pni, pnk, isize, 2)
   real(kind=4),    intent(out) :: rscavg(pni, pnk, 4, 2)
   real(kind=4),    intent (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent (in) :: roarow(pni, pnk)
   real(kind=4),    intent (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent (in) :: qr    (pni, pnk, 2)
   real(kind=4),    intent(out) :: qr_vel(pni, pnk, 2)
   real(kind=4),    intent (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent (in) :: pres  (pni, pnk)
   real(kind=4),    intent (in) :: amu   (pni, pnk)
   real(kind=4),    intent (in) :: amfp  (pni, pnk)
!!if_off
!
!  Local variables
!
   integer(kind=4) :: n, l, i
   real(kind=4)    :: amuw, tl, schm, dg0, dg_h2o2, dg_rooh
   real(kind=4)    :: dg_hno3, dg_nh3, prii, priiv, cfac, re, x, tempc, sigma
   real(kind=4)    :: bond, phys, physp, rephalf, aaa
   real(kind=4)    :: sh_h2o2, sh_rooh, sh_hno3, sh_nh3
   real(kind=4)    :: rn0s, st, rr, vr, sstar, colimp, pin0s, rlambds
   real(kind=4)    :: wd, LOGwd, LOGA, B, LOGC, D, colef1r, colef2r, colef1s, colef2s
   real(kind=4)    :: rrm(pni,2), alpha, gamma, relog

   real(kind=4)    :: mwt_h2o2, mwt_rooh, mwt_hno3, mwt_nh3
   real(kind=4)    :: rel_mwt_h2o2, rel_mwt_rooh, rel_mwt_hno3, rel_mwt_nh3
!  updated based on Seinfeld & Pandis (1998) (Allen & Raabe, 1982)
   real(kind=4), parameter :: aa1 = 1.257, aa2 = 0.4, aa3 = 1.1
!  universal gas constant (in g m2 s-2 mole-1 K-1)
   real(kind=4), parameter :: rhorain = 1000.0
!  parameters for calculating RE for drops (20 um < Dp <= 1 mm)
!  (Jacobson, 1999; Beard, 1976)
   real(kind=4), parameter :: b0 = -3.18657, b1 = 0.992696, b2 = -0.00153193
   real(kind=4), parameter :: b3 = -0.000987059, b4 = -0.000578878, b5 = 0.0000855176
   real(kind=4), parameter :: b6 = -0.00000327815
!  parameters for calculating re for drops (1 mm < dp <= 7 mm)
!  (jacobson, 1999; beard, 1976)
   real(kind=4), parameter :: bb0 = -0.500015e+01, bb1 = 0.523778e+01, bb2 = -0.204914e+01
   real(kind=4), parameter :: bb3 = 0.475294, bb4 = -0.542819e-01, bb5 = 0.238449e-02
!  a small number to avoid overflow
   real(kind=4), parameter :: smf = 1.0e-30

!  parameters for calculating rain aerosol scavenging coefficient:

   real(kind=4), parameter :: ar0 = -6.2609, ar1 = 0.6820, ar2 = 0.8676, ar3 = 0.1282
   real(kind=4), parameter :: br0 = -14.7073, br1 = 51.0427, br2 = -97.3056, br3 = 97.9458, br4 = -53.9227, br5 = 15.3109, br6 = -1.7510
   real(kind=4), parameter :: cr0 = 0.7230, cr1 = 0.0303
   real(kind=4), parameter :: dr0 = -0.6492, dr1 = 9.3483, dr2 = -21.9286, dr3 = 25.3167, dr4 = -15.3951, dr5 = 4.7242, dr6 = -0.5766

!  parameters for calculating snow aerosol scavenging coefficient:

   real(kind=4), parameter :: as0 = -4.4260, as1 = 1.3940, as2 = -1.2020, as3 = -3.2942, as4 = -1.9521, as5 = -0.4904, as6 = -0.0457
   real(kind=4), parameter :: bs0 = -4.3531, bs1 = -0.7828, bs2 = 12.7675, bs3 = -19.8635, bs4 = 13.6176, bs5 = -4.4350, bs6 = 0.5551
   real(kind=4), parameter :: cs0 = 0.5664, cs1 = 0.0085, cs2 = -0.1948, cs3 = -0.6532, cs4 = -0.5462, cs5 = -0.1778, cs6 = -0.0201
   real(kind=4), parameter :: ds0 = 0.5689, ds1 = -0.0923, ds2 = 0.0402, ds3 = 1.4523, ds4 = -2.0780, ds5 = 1.0500, ds6 = -0.18210

   alpha    = 0.0
   amuw     = 1.002e-3  !at 20 c [kg /m/sec]
   cfac     = 0.0
   colef    = 0.0
   colimp   = 0.0
   rlambds  = 0.0
   rn0s     = 0.0
   rrm      = 0.0
   tl       = 0.0
   vr       = 0.0

   rscavg = 0.0
   qr_vel = 0.0

!  molecular weight for h2o2, rooh, hno3, nh3
   mwt_h2o2 = mwt_igs(igs_H2O2)
   mwt_rooh = mwt_igs(igs_ROOH)
   mwt_hno3 = mwt_igs(igs_HNO3)
   mwt_nh3  = mwt_igs(igs_NH3)

   rel_mwt_h2o2 = sqrt((mwt_air + mwt_h2o2) / mwt_h2o2)
   rel_mwt_rooh = sqrt((mwt_air + mwt_rooh) / mwt_rooh)
   rel_mwt_hno3 = sqrt((mwt_air + mwt_hno3) / mwt_hno3)
   rel_mwt_nh3  = sqrt((mwt_air + mwt_nh3)  / mwt_nh3)

   do l = 1, pnk
      do i = 1, pni

         tl = tp(i, l)

!  collector drop size (liquid precipitation)
         if (qr(i, l, 1) > 1.0e-15) then

!  for rain: mass mean raindrop radius (rrm), factor 1.e-3 converts mm into m.
            rrm(i, 1) = 0.35 * (qr(i, l, 1) * 3600.0) ** 0.25 * 1.0e-3  !assuming qr in kg
!  note qr is in m s-1 from input, while the precip. rate
!  needed in the following formula requires qr in mm hr-1!
!  rrm(i, 1) = 0.35 * (qr(i, l, 1) * 3600. * 1000.)**0.25 * 1.e-3
            prii = 2.0 / 9.0 * grav / amu(i, l)
            priiv = prii * (rhorain - roarow(i, l))

!  cunningham slip correction factor settling velocity
            cfac = 1.0 + amfp(i, l) / rrm(i,1) * (aa1 + aa2 * exp(-aa3 * rrm(i,1) / amfp(i, l)))
            qr_vel(i, l, 1) = priiv * rrm(i,1) ** 2 * cfac
            re = roarow(i, l) * qr_vel(i, l, 1) * 2.0 * rrm(i, 1) / amu(i, l)
!  re and qr_vel (prec. velocity) for moderate drops (i.e. 20 um < dp <1 mm),
!  and for large drops (i.e. 1 mm < dp < 7 mm)
!  Beard (1976), Jacobson (1999)
            if (rrm(i,1) > 20.0e-06 .and. rrm(i,1) <= 0.0005) then
               x = log(32.0 * rrm(i,1) ** 3 / 3.0 / amu(i, l) ** 2 *   &
                   (1000.0 - roarow(i, l)) * roarow(i, l) * 9.81)
               re = cfac * exp(b0 + b1 * x + b2 * x * x + b3 * x * x * x + b4 * x * x * x * x +  &
                    b5 * x * x * x * x * x + b6 * x * x * x * x * x * x)
               qr_vel(i, l,1) = re * amu(i, l) / roarow(i, l) / 2.0 / rrm(i,1)
            end if
            if (rrm(i,1) > 0.0005) then
               tempc = max((tl - 273.15),0.)
               sigma = (76.1 - 0.155 * tempc) * 0.001
               bond = 16.0 / 3.0 * (1000.0 - roarow(i, l)) * 9.81 / sigma * rrm(i,1) ** 2
               phys = sigma ** 3 * roarow(i, l) ** 2 / (amu(i, l) ** 4 * (1000.0 - roarow(i, l)) * 9.81)
               physp = phys ** (1.0 / 6.0)
               x = log(bond * physp)
               re = physp * exp(bb0 + bb1 * x + bb2 * x * x + bb3 * x * x * x +  &
                    bb4 * x * x * x * x + bb5 * x * x * x * x * x)
               qr_vel(i, l,1) = re * amu(i, l) / roarow(i, l) / 2.0 / rrm(i,1)
            end if
!
!  gas diffusivity for h2o2, rooh, hno3, nh3,
!  note 3 / (8. * a * d**2) = 3.0751e-06, with a = 6.02213e+23 and d = 4.5e-10 m
            dg0 = 3.0751e-06 / roarow(i, l) * 0.001 * sqrt(rgasi * 1.0E3 * tl * mwt_air / 2.0 / pi)
            dg_h2o2 = dg0 * rel_mwt_h2o2
            dg_rooh = dg0 * rel_mwt_rooh
            dg_hno3 = dg0 * rel_mwt_hno3
            dg_nh3  = dg0 * rel_mwt_nh3
!  Sherwood number
            rephalf = sqrt(re)
            aaa = amu(i, l) / roarow(i, l)
            sh_h2o2 = 2.0 + 0.6 * rephalf * (aaa / dg_h2o2) ** 0.3333
            sh_rooh = 2.0 + 0.6 * rephalf * (aaa / dg_rooh) ** 0.3333
            sh_hno3 = 2.0 + 0.6 * rephalf * (aaa / dg_hno3) ** 0.3333
            sh_nh3  = 2.0 + 0.6 * rephalf * (aaa / dg_nh3)  ** 0.3333
!  Rain scavenging rate for soluble gases (h2o2, rooh, hno3, nh3)
!  Based on gamma = 3/2*qr/ut*sh*dg/dp**2, derived from Seinfeld & Pandis (1998)
            gamma = (3.0 / 2.0) * (qr(i, l, 1) * 1.0e-3 / qr_vel(i, l, 1)) / (rrm(i, 1)**2)
            rscavg(i, l, 1, 1) = sh_h2o2 * dg_h2o2 * gamma
            rscavg(i, l, 2, 1) = sh_rooh * dg_rooh * gamma
            rscavg(i, l, 3, 1) = sh_hno3 * dg_hno3 * gamma
            rscavg(i, l, 4, 1) = sh_nh3 * dg_nh3 * gamma
! particles
            if(chm_wdep_scav_coef_s /= 'WANG2014') then
!  Reynold number
               re = rrm(i, 1) * qr_vel(i, l, 1) * roarow(i, l) / amu(i, l)
               rephalf = sqrt(re)
               relog = log(1.0 + re)
               sstar = (1.2 + (1.0 / 12.0) * relog) / (1.0 + relog)
               do n = 1, isize
!  Ratio of radius of collected particle and collector drop
                  rr   = rhsize(i, l, n) / rrm(i, 1)
!  stokes friction and diffusion coefficients.
                  schm = amu(i, l) / pdiff(i, l, n) / roarow(i, l)
!  Stokes number of collected particels
                  st = 2.0 * pdepv(i, l, n) / grav * (qr_vel(i, l, 1) - pdepv(i, l, n)) / (2.0 * rrm(i,1))
                  if (st > sstar) then
                     colimp = ((st - sstar) / (st - sstar + 2.0 / 3.0)) ** (3.0 / 2.0) * sqrt(1000.0 / rhop(i, l, n))
                  else
                     colimp = 0.0
                  end if
                  vr    = amuw / amu(i, l)
!  Rain scavenging efficiency
                  colef(i, l, n, 1) = 4.0 / (re * schm) * (1.0 + 0.4 * rephalf * schm **  &
                                   (1.0 / 3.0) + 0.16 * sqrt(re * schm)) + 4.0 * rr *  &
                                   (1.0 / vr + (1.0 + 2.0 * rephalf) * rr) + colimp

!  Set the upbound for collection efficiency
                  colef(i, l, n, 1) = max(0.0, min(1.00, colef(i, l, n, 1)))
               end do
            else
!  new aerosol rain scavenging coefficient (s-1) based on Wang et. al (2014)
               do n = 1, isize
!  wet diameter in micrometer
                  wd = min((rhsize(i, l, n) * 2.0 * 1.0e+6), 100.)
                  LOGwd = log10(wd)

!  Rain scavenging efficiency
                  if(wd < 2.2) then
                     LOGA = ar0 + ar1 * LOGwd + ar2 * (LOGwd**2.0) + ar3 * (LOGwd**3.0)
                     B    = cr0 + cr1 * LOGwd
                  else
                     LOGA = br0 + br1 * LOGwd + br2 * (LOGwd**2.0) + br3 * (LOGwd**3.0) + br4 * (LOGwd**4.0) + br5 * (LOGwd**5.0) + br6 * (LOGwd**6.0)
                     B    = dr0 + dr1 * LOGwd + dr2 * (LOGwd**2.0) + dr3 * (LOGwd**3.0) + dr4 * (LOGwd**4.0) + dr5 * (LOGwd**5.0) + dr6 * (LOGwd**6.0)
                  endif

                  colef1r = 10.0 ** LOGA
                  colef2r = (qr(i, l, 1) * 3600.0)**B
                  colef(i, l, n, 1) = colef1r * colef2r

!  Set the upbound for collection efficiency
                  colef(i, l, n, 1) = max(0.0, min(1.00, colef(i, l, n, 1)))
               end do
            end if
         end if

! solid precipitation
         if (qr(i, l, 2) > 1.0e-15) then
!            if (tl <= t1s .and. tl >= t1s - 8.0) then
            if (tl >= t1s - 8.0) then
!  Needle snow scavenging
!  data from slinn(1984) in atmospheic science and power production, ed. darryl randerson.
!  for snow scavenging, the density of snow is set as 1/10 of liquid water.
!  The factor 1.0e-2 in the wetdep calculation takes this into account plus the unit change
!  into m s-1
               qr_vel(i, l, 2) = 50.0e-2    ! average settling velocity
               rrm(i, 2)       = 10.0e-6    ! characteristic capture length
               alpha        = 1.0
!  parameter rn0s used for scavenging of gases (hno3 and nh3)
               rn0s         = 0.05
            end if
            if (tl < t1s - 8.0 .and. tl >= t1s - 25.0) then
!  Steller snow scavenging
               qr_vel(i, l, 2) = 57.0e-2
               rrm(i, 2)       = 100.0e-6
               alpha        = 0.5
               rn0s         = 0.1        ! for gas scavenging
            end if
            if (tl < t1s - 25.0) then
!  Graupel scavenging
               qr_vel(i, l, 2) = 180.0e-2
               rrm(i, 2)       = 1000.0e-6
               alpha        = 2.0 / 3.0
               rn0s         = 1.0        ! for gas scavenging
            end if

!  Snow scavenging rates for gases (hno3 and nh3) same as in adom, where the scavenging rates
!  for hno3 and nh3 are set to be twice the rate for h2so4, which is similar to the rate
!  of collection of supercool water drops by snow/ice.
            pin0s = pi * rn0s
            rlambds = (pin0s * 1.31 / (qr(i, l, 2) / 1000.0 + smf)) ** 0.25 !qr in kg m-2 s-1
            rscavg(i, l, 3, 2) = pin0s / 2.0 * 0.04 * (1.013e5 / pres(i, l)) ** 0.4  &
                                 * (2.22 * 131.01 / rlambds ** 3.11)
            rscavg(i, l, 4, 2) = rscavg(i, l, 3, 2)
            rscavg(i, l, 1, 2) = 0.0
            rscavg(i, l, 2, 2) = 0.0
! particles
            if (chm_wdep_scav_coef_s /= 'WANG2014') then
!  Reynold number

               re = rrm(i, 2) * qr_vel(i, l, 2) * roarow(i, l) / amu(i, l)
               rephalf = sqrt(re)
               relog = log(1.0 + re)
               sstar = (1.2 + (1.0 / 12.0) * relog) / (1.0 + relog)
               do n = 1, isize
!  Ratio of radius of collected particle and collector drop
                  rr   = rhsize(i, l, n) / rrm(i, 2)
!  stokes friction and diffusion coefficients.
                  schm = amu(i, l) / pdiff(i, l, n) / roarow(i, l)
!  Stokes number of collected particels
                  st = 2.0 * pdepv(i, l, n) / grav * (qr_vel(i, l, 2) - pdepv(i, l, n)) / (2.0 * rrm(i, 2))
                  if (st > sstar) then
                     colimp = ((st - sstar) / (st - sstar + 2.0 / 3.0)) ** (3.0 / 2.0) * sqrt(1000.0 / rhop(i, l, n))
                  else
                     colimp = 0.0
                  end if
!  Snow scavenging efficiency
                  colef(i, l, n, 2) = (1.0 / schm) ** alpha +  &
                                      (1.0 - exp(-(1.0 + rephalf) * rr ** 2)) + colimp
!  Set the upbound for collection efficiency
                  colef(i, l, n, 2) = max(0.0, min(1.00, colef(i, l, n, 2)))

               end do
            else
!  new aerosol snow scavenging coefficient (s-1) based on Wang et. al (2014)
               do n = 1, isize
!  wet diameter in micrometer
                  wd = min((rhsize(i, l, n) * 2.0 * 1.0e+6), 100.)
                  LOGwd = log10(wd)
                  if (wd < 1.44) then
                     LOGC = as0 + as1 * LOGwd + as2 * (LOGwd**2.0) + as3 * &
                          (LOGwd**3.0) + as4 * (LOGwd**4.0) + as5 * &
                          (LOGwd**5.0) + as6 * (LOGwd**6.0)
                     D    = cs0 + cs1 * LOGwd + cs2 * (LOGwd**2.0) + cs3 * &
                          (LOGwd**3.0) + cs4 * (LOGwd**4.0) + cs5 * &
                          (LOGwd**5.0) + cs6 * (LOGwd**6.0)
                  else
                     LOGC = bs0 + bs1 * LOGwd + bs2 * (LOGwd**2.0) + bs3 * &
                          (LOGwd**3.0) + bs4 * (LOGwd**4.0) + bs5 * &
                          (LOGwd**5.0) + bs6 * (LOGwd**6.0)
                     D    = ds0 + ds1 * LOGwd + ds2 * (LOGwd**2.0) + ds3 * &
                          (LOGwd**3.0) + ds4 * (LOGwd**4.0) + ds5 * &
                          (LOGwd**5.0) + ds6 * (LOGwd**6.0)
                  end if
                  colef1s = 10.0 ** LOGC
                  colef2s = (qr(i, l, 2) * 3600.0)**D
                  colef(i, l, n, 2) = colef1s * colef2s
!  Set the upbound for scavenging coefficient
                  colef(i, l, n, 2) = max(0.0, min(1.00, colef(i, l, n, 2)))

               end do
            end if
         end if
      end do
   end do

   return
end subroutine mach_cam_cas_split



