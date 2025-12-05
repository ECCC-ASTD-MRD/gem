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
! Fichier/File   : mach_cam_rain.ftn90
! Creation       : S. Gong, W. Gong, S. Gravel and B. Pabla for GEM-MACH, June 2008
! Description    : Scavenging calculation by rain and snow
!
! Extra info     : - First version created by S. Gong Jul 08 1994 for CAM
!                    Using scheme of Slinn (1977) - Water, Air, and Soil Pollution 7(1977)
!                  - Vectorized the whole program and add someworking spaces. (S. Gong, Dec 19, 1996)
!                  - Adding scavenging of gases (W. Gong, May 2001)
!
! Arguments:  IN
!                 tp       -> Temperature (K)
!                 qr       -> rain water/snow
!                 rhsize   -> Wet radius
!                 pdepv    -> gravitational settling velocity
!                 pres     -> Local mid-layer pressure [Pa]
!                 roarow   -> Air density (kg/m3)
!                 pdiff    -> diffusion parameter
!                 rhop     -> Final wet density
!                 cc2d     -> 2d cloud cover above
!                 amu      -> Air's dynamic viscosity
!                 amfp     -> Mean molecular free path
!
!             IN/OUT
!                 xrow     -> tracers concentration
!
!             OUT
!                 rtbcld   -> Rain Scavenging
!                 qr_vel   -> Average precipitation settling velocity
!
!============================================================================
!
!!if_on
subroutine mach_cam_rain(tp, qr, rhsize, pdepv, pres, roarow, pdiff, &
                         amu, amfp, xrow, rtbcld, rhop, cc2d, qr_vel, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
!!if_off
   use chm_utils_mod,        only: chm_timestep
   use mach_cam_headers_mod, only: mach_cam_cas, mach_cam_cas_split
   use chm_consphychm_mod,   only: t1s
   use chm_nml_mod,          only: chm_split_snow_rain_l, chm_wdep_scav_coef_s
   use mach_cam_utils_mod,   only: icom, iae1, igs_H2O2, igs_ROOH, igs_HNO3, igs_NH3
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: tp    (pni, pnk)
   real(kind=4),    intent   (in) :: qr    (pni, pnk, 2)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pres  (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent  (out) :: rtbcld(pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: cc2d  (pni, pnk)
   real(kind=4),    intent  (out) :: qr_vel(pni, pnk, 2)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: amfp  (pni, pnk)
!!if_off
!
! Local variables
!
   integer(kind=4)         :: l, i, n, nn, np
   real(kind=4)            :: tl, rrm, xold, xnew, tmin, dm, tend
   real(kind=4)            :: rscavg(pni, pnk, 4, 2)
   real(kind=4), parameter :: bcrain = 0.5, bcsnow = 0.15
   real(kind=4)            :: colef (pni, pnk, isize, 2)
   real(kind=4)            :: wetdep(pni, pnk, isize, 2)

   rtbcld   = 0.0
!  qr_vel = 0.0

!
!  Rain scavenging
!
!  call to compute collection efficientcy coefficients as well as the scavenging
!  rates for soluble gases
!  Note that in the case of WANG2014 semi-empirical scavenging parameterization
!  it is the scavenging rate for particles being computed in mach_cam_cas_split.

!  minimum tendency for below cloud scavenging:
   tmin = 0.0
!  options for split snow-rain scavenging or combined scavenging (orig)
   if (chm_split_snow_rain_l .or. (chm_wdep_scav_coef_s == 'WANG2014')) then
      call mach_cam_cas_split(TP, COLEF, rscavg, RHOP, ROAROW, RHSIZE, PRES, &
                              QR, qr_vel, PDIFF, PDEPV, amu, amfp, pni, pnk)
!  Note qr_vel now contains fall velocity at a given grid,
!  and passed back to SCAVENG to be used for equilibrium scavenging
!  of SO2 and CO2. (wg/June/2001)
      wetdep = 0.0
      if(chm_wdep_scav_coef_s == 'WANG2014') then
         wetdep = colef
      else
         do l = 1, pnk
            do i = 1, pni
               tl = tp(i, l) - t1s
!  rain scavenging rate
               if (qr(i, l, 1) > 1.0e-15) then
!  the unit of qr is mm s-1 (or, kg-of-water m-2 s-1)
!  the 1.0e-3 converts rrm (average drop size) into m
                  rrm = 0.35 * (qr(i, l, 1) * 3600.0) ** 0.25 * 1.0e-3
                  do n = 1, isize
                     wetdep(i, l, n, 1) = bcrain * qr(i, l, 1) * 1.0e-3 * &
                                          colef(i, l, n, 1) / rrm
                  enddo
               end if
!  for snow scavenging, the dnesity of snow is set as 1/10 of liquid water.
!  the factor 1.0e-2 in the wetdepcalculation takes this into account plus the unit change
!  into m s-1
               if (qr(i, l, 2) > 1.0e-15) then
!  Needle snow scavenging
                  if (tl >= -8.0) dm = 3.8e-5                 !characteristic length
!  steller snow scavenging
                  if (tl < -8.0 .and. tl >= -25.0) dm = 2.7e-5
!  graupel scavenging
                  if (tl < -25.0) dm = 1.4e-4
                  do n = 1, isize
                     wetdep(i, l, n, 2) = bcsnow * qr(i, l, 2) * 1.0e-2 * &
                                          colef(i, l, n, 2) / dm
                  end do
               end if
            end do
         end do
      end if


!  add the rain scavenging tendency
      do nn = 1, icom
         do n = 1, isize
            np = isize * (nn - 1) + n + (iae1 - 1)
            do l = 1, pnk
               do i = 1, pni
                  xold = xrow(i, l, np)
                  xnew = xold * exp(-chm_timestep * (wetdep(i, l, n, 1) + &
                                                     wetdep(i, l, n, 2)))
                  tend = min((xnew - xold) / chm_timestep, 0.0)
                  rtbcld(i, l, np) = tend * cc2d(i, l)
                  xrow(i, l, np) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
               end do
            end do
         end do
      end do

!  scavenging of soluble gas (H2O2, ROOH, HNO3, NH3)
      do l = 1, pnk
         do i = 1, pni
! h2o2
            xold = xrow(i, l, igs_H2O2)
            xnew = xold * exp(-chm_timestep * (rscavg(i, l, 1, 1) + &
                                               rscavg(i, l, 1, 2)))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_H2O2) = tend * cc2d(i, l)
            xrow(i, l, igs_H2O2) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
! rooh
            xold = xrow(i, l, igs_ROOH)
            xnew = xold * exp(-chm_timestep * (rscavg(i, l, 2, 1) + &
                                               rscavg(i, l, 2, 2)))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_ROOH) = tend * cc2d(i, l)
            xrow(i, l, igs_ROOH) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
! hno3
            xold = xrow(i, l, igs_HNO3)
            xnew = xold * exp(-chm_timestep * (rscavg(i, l, 3, 1) + &
                                               rscavg(i, l, 3, 2)))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_HNO3) = tend * cc2d(i, l)
            xrow(i, l, igs_HNO3) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
! nh3
            xold = xrow(i, l, igs_NH3)
            xnew = xold * exp(-chm_timestep * (rscavg(i, l, 4, 1) + &
                                               rscavg(i, l, 4, 2)))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_NH3) = tend * cc2d(i, l)
            xnew = max(xold + tend * chm_timestep, tmin)
            xrow(i, l, igs_NH3) = max(xold + tend * chm_timestep, tmin)      &
                                  * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
         end do
      end do

   else

      call mach_cam_cas(tp, colef, rscavg, rhop, roarow, rhsize, pres, &
                        qr, qr_vel, pdiff, pdepv, amu, amfp, pni, pnk)

!  Note qr_vel now contains fall velocity at a given grid,
!  and passed back to SCAVENG to be used for equilibrium scavenging
!  of SO2 and CO2. (wg/June/2001)
      wetdep = 0.0
      do l = 1, pnk
         do i = 1, pni
            if (qr(i, l, 1) > 1.0e-15) then
               tl = tp(i, l) - t1s
               if (tl > 0.0) then
!  rain scavenging rate
!  the unit of qr is mm s-1 (or, kg-of-water m-2 s-1)
!  the 1.0e-3 converts rrm (average drop size) into m
                  rrm = 0.35 * (qr(i, l, 1) * 3600.0) ** 0.25 * 1.0e-3
                  do n = 1, isize
                     wetdep(i, l, n, 1) = bcrain * qr(i, l, 1) * 1.0e-3 * &
                                          colef(i, l, n, 1) / rrm
                  end do
               else
!  for snow scavenging, the dnesity of snow is set as 1/10 of liquid water.
!  the factor 1.0e-2 in the wetdep calculation takes this into account
!  plus the unit change into m s-1
               ! characteristic length
!  Needle snow scavenging
                  if (tl <= 0.0 .and. tl >= -8.0) dm = 3.8e-5
!  steller snow scavenging
                  if (tl < -8.0 .and. tl >= -25.0) dm = 2.7e-5
!  graupel scavenging
                  if (tl < -25.0) dm = 1.4e-4
!
                  do n = 1, isize
                     wetdep(i, l, n, 1) = bcsnow * qr(i, l, 1) * 1.0e-2 * &
                                          colef(i, l, n, 1) / dm
                  end do
               end if
            end if
         end do
      end do

!  add the rain scavenging tendency
      do nn = 1, icom
         do n = 1, isize
            np = isize * (nn - 1) + n + (iae1 - 1)
            do l = 1, pnk
               do i = 1, pni
                  xold = xrow(i, l, np)
                  xnew = xold * exp(-chm_timestep * wetdep(i, l, n, 1))
                  tend = min((xnew - xold) / chm_timestep, 0.0)
                  rtbcld(i, l, np) = tend * cc2d(i, l)
                  xrow(i, l, np) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
               end do
            end do
         end do
      end do

!  scavenging of soluble gas (H2O2, ROOH, HNO3, NH3)

      do l = 1, pnk
         do i = 1, pni
! h2o2
            xold = xrow(i, l, igs_H2O2)
            xnew = xold * exp(-chm_timestep * rscavg(i, l, 1, 1))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_H2O2) = tend * cc2d(i, l)
            xrow(i, l, igs_H2O2) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
! rooh
            xold = xrow(i, l, igs_ROOH)
            xnew = xold * exp(-chm_timestep * rscavg(i, l, 2, 1))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_ROOH) = tend * cc2d(i, l)
            xrow(i, l, igs_ROOH) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
! hno3
            xold = xrow(i, l, igs_HNO3)
            xnew = xold * exp(-chm_timestep * rscavg(i, l, 3, 1))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_HNO3) = tend * cc2d(i, l)
            xrow(i, l, igs_HNO3) = max(xold + tend * chm_timestep, tmin)     &
                                   * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
! nh3
            xold = xrow(i, l, igs_NH3)
            xnew = xold * exp(-chm_timestep * rscavg(i, l, 4, 1))
            tend = min((xnew - xold) / chm_timestep, 0.0)
            rtbcld(i, l, igs_NH3) = tend * cc2d(i, l)
            xnew = max(xold + tend * chm_timestep, tmin)
            xrow(i, l, igs_NH3) = max(xold + tend * chm_timestep, tmin)      &
                                  * cc2d(i, l) + xold * (1.0 - cc2d(i, l))
         end do
      end do
   end if  ! chm_wdep_split_snow_rain options

   return
end subroutine mach_cam_rain
