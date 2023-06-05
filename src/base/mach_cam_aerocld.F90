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
! Fichier/File   : mach_cam_aerocld.ftn90
! Creation       : S. Gong, W. Gong, S. Menard, V. Bouchet,  P. Huang, C. Stroud,
!                  S. Gravel and B. Pabla for GEM-MACH, June 2008
! Description    : Aerosol-cloud interaction module:
!                  (1) Sub-grid vertical velocity distribution.
!                  (2) Aerosol activation and cloud spectrum calculation.
!                  (3) Cloud chemistry of sulphate oxidation
!
! Extra info     : - First version created by S. Gong Jan 11 1998 for CAM
!                  - The in-cloud removal due to rain-out and production due
!                    to evporation of aerosols and gases are parameterized
!                    as the same fraction of cloud water is removed or
!                    produced in a cloud module. (S. Gong Dec 02 1998)
!                  - Implementation of the vectorize version of the ADOM
!                    stratus aqueous phase chemistry solver. This solver
!                    is vectorized along the 3 directions (X-Y-Z). This
!                    version of the solver is running one order of magnitude
!                    faster than the original ADOM scalar code. Because CAM
!                    is run by vertical slices (X-Z) along the Y coordinate,
!                    it is not possible at the moment to use the full potential
!                    of this solver. The latitudinal loop (Y) will be
!                    a scalar loop. (W. Gong, and S. Menard, Jun 01, 2000)
!                  - Adding mass redistribution and number diagnosis after
!                    the cloud chemistry process. (W. Gong, Nov 2000)
!                  - Implementing wet flux (including cloud-to-rain transfer
!                     and part of evaporation) (W. Gong, Apr 2001)
!                  - Implementing changes due to heterogeneous chemistry
!                    (V. Bouchet, Aug 2001)
!                  - Split aerosol OC into primary and secondary components
!                    (C. Stroud, Jul 2004) components !cs>>>
!                  - Added EC and CM for evaporation (W. Gong, Jul 2005)
!                  - Corrected molecular weight assignment for EC and CR
!                    (M. Moran, Sept 2005)
!
!
!**********************************************************************
!Author Sylvain Menard (AES/contractor)          August, 2000.
!       Wanmin Gong    (ARQI/AES)
!
!
!Revision v1.0 -  Sylvain Menard - August, 2000.
!       QI Model Application Team (MAT)
!
!Language
!       Fortran
!
!Object
!
!Arguments and more:
!
!____________________________________________________________________________
!          |                                      | T |           |    |
!  NAME    |          DESCRIPTION                 | Y |DIMENSIONS |IN/ |
!          |                                      | P |           |OUT |
!          |                                      | E |           |    |
!---------------------------------------------------------------------------
! throw    | Temperature (K)                      | R |pni, pnk   | I  |
! xrow     | Tracer array (kg/kg) with moon level | R |pni, pnk,  |    |
!          |                                      |   |ntr        |    |
! pnk      | No of vertical levels                | I | scalar    |    |
! pni      | number of longitude grid points      | I | scalar    |    |
! rhsize   | Unactivated ambient aerosol wet      | R | pni, pnk, |    |
!          | Radius  [m] (cf. CLSIZE)             |   | isize     |    |
! aeronum  | Number-Concentration of Aerosol      | R | pni, pnk, |    |
!          |                                      |   | isize     |    |
! ntr      | Number of tracers                    | I | scalar    | I  |
! icom     | Number of aerosol species            | I | scalar    |    |
! isize    | Number of size bins                  | I | scalar    |    |
! thlev    | Layer thickness [m]                  | R | pni, pnk  | I  |
! roarow   | Air Density  [kg/m3]                 | R | pni, pnk  | I  |
! pres     | pressure [Pa]                        | R | pni, pnk  | I  |
! zmlwc    | Cloud liquid water [kg/kg]           | R | pni, pnk  | I  |
! ccn      | cloud droplet number density [1/m3]  | R | pni, pnk  | I  |
! jlat     | slice numer (1-71)                   | I | scalar    | I  |
! rcrits   | Bin Number (+ Fraction Un-Activated) | R |pni, pnk,2 | I  |
! clsize   | Activated ambient aerosol wet        | R |pni, pnk,  |    |
!          | radius  [m] (cf. RHSIZE)             |   | isize     |    |
! q_bin    | CWC per activated size bin (kg/m3)   | R |pni, pnk,  |    |
!          |                                      |   |isize      |    |
! igs_H2O2 | pointer location of h2o2             | I | scalar    |    |
! igs_ROOH | pointer location of h2o2             | I | scalar    |    |
! igs_HNO3 | pointer location of hno3             | I | scalar    |    |
! igs_NH3  | pointer location of nh3              | I | scalar    |    |
! igs_O3   | pointer location of o3               | I | scalar    |    |
! igs_SO2  | pointer location of so2              | I | scalar    | I  |
!---------------------------------------------------------------------------
!
!!if_on
subroutine mach_cam_aerocld(throw, xrow, rhsize, aeronum, thlev, roarow, &
                            pres, zmlwc, jlat, rcrits, tcldcv, flux,     &
                            wetflx, fctr, frevp, rhrow, ccn, kount, pni, pnk)
   use mach_cam_utils_mod,      only: isize, ntr, nswdep
!!if_off
   use mach_cam_headers_mod,    only: mach_cam_intrsec1_outer, mach_cam_aeroact
   use mach_incld_headers_mod,  only: mach_incld_main
   use mach_hetv_headers_mod,   only: mach_hetv_hetchem
   use mach_cam_utils_mod,      only: maxnsg, maxns
   use mach_cam_utils_mod,      only: igs_SO2, igs_H2O2, igs_ROOH, mwt_igs
   use mach_cam_utils_mod,      only: rhop0, pvol, icom, iae1, ip_wflx, &
                                      iae_SU, iae_NI, iae_AM, mwt_aero
   use chm_nml_mod,             only: chm_hetchem_s, chm_aqueous_s
   use chm_utils_mod,           only: chm_error_l
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: jlat
   integer(kind=4), intent   (in) :: kount
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: throw  (pni, pnk)
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhsize (pni, pnk, isize)
   real(kind=4),    intent(inout) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev  (pni, pnk)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: zmlwc  (pni, pnk)
   real(kind=4),    intent  (out) :: rcrits (pni, pnk)
   real(kind=4),    intent   (in) :: tcldcv (pni, pnk)
   real(kind=4),    intent  (out) :: flux   (pni, pnk, nswdep)
   real(kind=4),    intent  (out) :: wetflx (pni, nswdep)
   real(kind=4),    intent   (in) :: fctr   (pni, pnk)
   real(kind=4),    intent   (in) :: frevp  (pni, pnk)
   real(kind=4),    intent   (in) :: rhrow  (pni, pnk)
   real(kind=4),    intent   (in) :: ccn    (pni, pnk)
!!if_off
!
! local variables
!
   integer(kind=4)  :: ii, ll, kk
   integer(kind=4)  :: ibulk
   integer(kind=4)  :: mm, nt, nn
   integer(kind=4)  :: nn1, nn2, nn3
   real(kind=4)     :: qbin_min, tramass
   real(kind=4)     :: clsize   (pni, pnk, isize)
   real(kind=4)     :: psacw    (pni, pnk)
   real(kind=4)     :: gaz      (pni, pnk, maxns)
   real(kind=4)     :: aerocon  (pni, pnk, icom, isize)
   real(kind=4)     :: aerocon0 (pni, pnk, icom, isize)
   real(kind=4)     :: q_bin    (pni, pnk, isize)
   real(kind=4)     :: q_bin12  (pni, pnk, isize, 2)
   real(kind=4)     :: daqchm   (pni, pnk, isize)
   real(kind=4)     :: rmass    (pni, pnk, isize)
   real(kind=4)     :: totmass  (pni, pnk, isize)
   real(kind=4)     :: rhopd    (pni, pnk, isize)
   real(kind=4)     :: zmlwc_new(pni, pnk)
   real(kind=4)     :: tmass    (pni, pnk, icom)
   real(kind=4)     :: evap     (pni, pnk, nswdep)
   real(kind=4)     :: fbin_cld (pni, pnk, icom, isize)
   real(kind=4)     :: evapp
   real(kind=4)     :: net_vol
!
! Code Begins 
!
!  switch for bulk chemistry or bin resolved
!  switch =1 for bulk chemistry
!  switch =0 for bin resolved
   ibulk = 1

!  Adding correction for cloud fraction (using total cloud fraction)
!  ( an arbitary number 1.e-10 is added to avoid overflow). WG, Feb, 2001
   zmlwc_new = 0.0
   where (tcldcv > 0.0)
      zmlwc_new = zmlwc / tcldcv
   end where

!  Aerosol activation and cloud spectrum [ghan scheme]; and
!  calculation of cwc per activated size bin
   call mach_cam_aeroact(q_bin, rhsize, rcrits, aeronum, zmlwc_new, &
                         roarow, clsize, ccn, ibulk, pni, pnk)

!  Aqueous phase chemistry is not executed if chm_aqueous_s == 'NIL'. So the following
!  aqueous phase chemistry code is executed only if chm_aqueous_s == 'GONG'. See also
!  mach_cam_main.ftn90
!
!----------------------------------------------------------------------------------------
!  Start of aqueous phase chemistry
!----------------------------------------------------------------------------------------

   if (chm_aqueous_s == 'GONG') then

      psacw    = 0.0
      qbin_min = 1e-07
!
!  Copy the gaseous species concentrations to gaz (for maxns species)
      do ll = 1, maxns
         do kk = 1, pnk
            do ii = 1, pni
               gaz(ii, kk, ll) = xrow(ii, kk, ll)
            end do
         end do
      end do

!  Similarly, copy the aerosol species into aerocon (isize bins)
      do mm = 1, isize
         do nt = 1, icom
            nn = (iae1 - 1) + (nt - 1) * isize + mm
            do kk = 1, pnk
               do ii = 1, pni
                  aerocon(ii, kk, nt, mm) = xrow(ii, kk, nn)
                  aerocon0(ii, kk, nt, mm) = xrow(ii, kk, nn)
               end do
            end do
         end do
      end do
!
!  Calculate dry density before cloud chemistry
      totmass = 0.0
      do mm = 1, isize
         do kk = 1, pnk
            do ii = 1, pni
               net_vol = 0.0
               do nt = 1, icom
                  tramass = max(1.0e-33, aerocon(ii, kk, nt, mm))
                  totmass(ii, kk, mm) = totmass(ii, kk, mm) + tramass
                  net_vol = net_vol + tramass / rhop0(nt)
               end do
               rhopd(ii, kk, mm) = totmass(ii, kk, mm) / net_vol
            end do
         end do
      end do
!
!  transfer q_bin into q_bin12
      do mm = 1, isize
         do kk = 1, pnk
            do ii = 1, pni
               q_bin12(ii, kk, mm, 1) = q_bin(ii, kk, mm) ! water only
               q_bin12(ii, kk, mm, 2) = 0.0               ! ice/snow not used
            end do
         end do
      end do

!  Assign q_bin=0 if q_bin <1e-07 kg/m3.

      where (q_bin12 <= qbin_min)
          q_bin12 = 0.0
      end where
!      
!  Section 3  new aqueous phase chemistry solver to be run by slices.
!  Sylvain M. july 2000.
      call mach_incld_main(gaz, aerocon, q_bin12, throw, psacw, clsize, rcrits, &
                           roarow, ibulk, flux, fctr, aeronum, pni, pnk)
      if (chm_error_l) return

!  Return updated gaz back into xrow
!  adjusted to take into account of cloud fraction - WG, Feb, 2001
      do ll = 1, maxns
         do kk = 1, pnk
            do ii = 1, pni
               xrow(ii, kk, ll) = (1.0 - tcldcv(ii, kk)) * xrow(ii, kk, ll) +  &
                                       tcldcv(ii, kk) * gaz(ii, kk, ll)
            end do
         end do
      end do

!  Calculate daqchm, the net change in aerosol mass needed for rebinning
!  making use of xrow before it is updated.
!  And copy the updated "aerocon" back into "xrow"
      daqchm  = 0.0
      do mm = 1, isize
         do nt = 1, icom
            nn = (iae1 - 1) + (nt - 1) * isize + mm
            do kk = 1, pnk
               do ii = 1, pni
                  daqchm(ii, kk, mm) = daqchm(ii, kk, mm) + &
                                   (aerocon(ii, kk, nt, mm) - xrow(ii, kk, nn))
                  xrow(ii, kk, nn) = aerocon(ii, kk, nt, mm)
               end do
            end do
         end do
      end do


!  mass redistribution and number diagnosis
      call mach_cam_intrsec1_outer(xrow, rhopd, daqchm, aeronum, q_bin, &
                                   rcrits, 1, pni, pnk)

!  In-cloud size distribution ratios (used to distribute
!  evaporated bulk aerosol mass to bins).  WG, May 2001
      tmass = 0.0
      do nt = 1, icom
         do mm = 1, isize
            nn = (iae1 - 1) + (nt - 1) * isize + mm
            do kk = 1, pnk
               do ii = 1, pni
                  tramass = max(1.0e-33, xrow(ii, kk, nn))
                  tmass(ii, kk, nt) = tmass(ii, kk, nt) + tramass
               end do
            end do
         end do
         do mm = 1, isize
            nn = (iae1 - 1) + (nt - 1) * isize + mm
            do kk = 1, pnk
               do ii = 1, pni
                  fbin_cld(ii, kk, nt, mm) = max(0.0, xrow(ii, kk, nn) &
                                             / tmass(ii, kk, nt))
!
!  Adjust for cloud fraction - WG, Feb. 2001
                  xrow(ii, kk, nn) = tcldcv(ii, kk) * xrow(ii, kk, nn) +  &
                                     (1.0 - tcldcv(ii, kk)) * &
                                     aerocon0(ii, kk, nt, mm)
               end do
            end do
         end do
      end do

!  wet flux from cloud-to-rain (moles per m2)
!  add evaporation (note k=1 at model top) - WG, May 2001
      wetflx = 0.0

      do ll = 1, nswdep
         do kk = 1, pnk
            do ii = 1, pni
               wetflx(ii, ll) = wetflx(ii, ll) + flux(ii, kk, ll) *  &
                                tcldcv(ii, kk) * thlev(ii, kk)
!  evap in moles per m3:
               evap(ii, kk, ll) = wetflx(ii, ll) * frevp(ii, kk) / thlev(ii, kk)
               wetflx(ii, ll) = wetflx(ii, ll) * (1.0 - frevp(ii, kk))
            end do
         end do
      end do

!  evaporation of particles and size distribution
!  update gas (so2, h2o2, rooh) and particles (all in kg/kg) after evaporation
      do kk = 1, pnk
         do ii = 1, pni
            xrow(ii, kk, igs_SO2)  = xrow(ii, kk, igs_SO2) + evap(ii, kk, 1) * &
                                     mwt_igs(igs_SO2) * 1.0e-3 / roarow(ii, kk)
            xrow(ii, kk, igs_H2O2) = xrow(ii, kk, igs_H2O2) + evap(ii, kk, 2)* &
                                     mwt_igs(igs_h2o2) * 1.0e-3 / roarow(ii, kk)
            xrow(ii, kk, igs_ROOH) = xrow(ii, kk, igs_ROOH) + evap(ii, kk, 3)* &
                                     mwt_igs(igs_ROOH) * 1.0e-3 / roarow(ii, kk)
         end do
      end do
      
      do nt = 1, icom
         ll = ip_wflx(nt)
         do mm = 1, isize
            nn = (iae1 - 1) + (nt - 1) * isize + mm
            do kk = 1, pnk
               do ii = 1, pni
                  evapp = evap(ii, kk, ll) * fbin_cld(ii, kk, nt, mm)
                  xrow(ii, kk, nn) = xrow(ii, kk, nn) + evapp * &
                                     mwt_aero(nt) * 1.0e-3 / roarow(ii, kk)
               end do
            end do
         end do
      end do
!
   else
      flux   = 0.0
      wetflx = 0.0
   end if ! end of aqueous phase chemistry
!-------------------------------------------------------------------------------
!  end of aqueous phase chemistry
!-------------------------------------------------------------------------------
!
!  Call heterogeneous chemistry (apply to cases when aq & het. are run 
!  independantly of each other or when it is run for the cloud free portion of
!  the grid - then the cloudy het. chem is done in IN_cloud)
!  heterogeneous chemistry called independently of aqueous chemistry for now
   if (chm_hetchem_s /= 'NIL') then
!-------------------------------------------------------------------------------
!  start  of heterogeneous chemistry
!-------------------------------------------------------------------------------
!  mass balance checks for hetchem done inside mach_hetv_hetchem

!  update RHOPD and AERONUM
      totmass = 0.0
      do mm = 1, isize
         do kk = 1, pnk
            do ii = 1, pni
               net_vol = 0.0
               do nt = 1, icom
                  nn = (iae1 - 1) + (nt - 1) * isize + mm
                  tramass = max(1.0e-33, xrow(ii, kk, nn))
                  totmass(ii, kk, mm) = totmass(ii, kk, mm) + tramass
                  net_vol             = net_vol + tramass / rhop0(nt)
               end do
               rhopd(ii, kk, mm) = totmass(ii, kk, mm) / net_vol
               rmass(ii, kk, mm) = pvol(mm) * rhopd(ii, kk, mm)
               aeronum(ii, kk, mm) = totmass(ii, kk, mm) / rmass(ii, kk, mm)
            end do
         end do
      end do
!
!  xrow is transfered into gaz (for maxns gazeous species)
      do ll = 1, maxns
         do kk = 1, pnk
            do ii = 1, pni
               gaz(ii, kk, ll) = xrow(ii, kk, ll)
            end do
         end do
      end do

!  xrow is transfer into aerocon for aerosol species
      do mm = 1, isize
         nn1 = (iae1 - 1) + (iae_SU - 1) * isize + mm
         nn2 = (iae1 - 1) + (iae_NI - 1) * isize + mm
         nn3 = (iae1 - 1) + (iae_AM - 1) * isize + mm
         do kk = 1, pnk
            do ii = 1, pni
               aerocon(ii, kk, iae_SU, mm) = xrow(ii, kk, nn1) !so4
               aerocon(ii, kk, iae_NI, mm) = xrow(ii, kk, nn2) !no
               aerocon(ii, kk, iae_AM, mm) = xrow(ii, kk, nn3) !nh
            end do
         end do
      end do

      call mach_hetv_hetchem(gaz, aerocon, throw, pres, aeronum, &
                             rhrow, roarow, ibulk, jlat, kount, pni, pnk)
      if (chm_error_l) return

!  transfer gaz back to xrow for gaseous species 
      do ll = 1, maxns   ! transfer first maxns species
         do kk = 1, pnk ! or maxcnz equivalent
            do ii = 1, pni
               xrow(ii, kk, ll) = gaz(ii, kk, ll)
            end do
         end do
      end do
!
!daqchm:
!  calculating the net change in aerosol mass needed for rebinning
!  making use of xrow before it is updated; and
!
!  transfer "aerocon" back to "xrow"
      daqchm = 0.0
      do mm = 1, isize
         nn1 = (iae1 - 1) + (iae_SU - 1) * isize + mm
         nn2 = (iae1 - 1) + (iae_NI - 1) * isize + mm
         nn3 = (iae1 - 1) + (iae_AM - 1) * isize + mm
         do kk = 1, pnk
            do ii = 1, pni
               daqchm(ii, kk, mm) = &
                           aerocon(ii, kk, iae_SU, mm) - xrow(ii, kk, nn1) + &
                           aerocon(ii, kk, iae_NI, mm) - xrow(ii, kk, nn2) + &
                           aerocon(ii, kk, iae_AM, mm) - xrow(ii, kk, nn3)
!
!  transfer "aerocon" back to "xrow"
               xrow(ii, kk, nn1) = aerocon(ii, kk, iae_SU, mm) !so4
               xrow(ii, kk, nn2) = aerocon(ii, kk, iae_NI, mm) !no3
               xrow(ii, kk, nn3) = aerocon(ii, kk, iae_AM, mm) !nh4
            end do
         end do
      end do

!  mass redistribution
      call mach_cam_intrsec1_outer(xrow, rhopd, daqchm, aeronum, q_bin, &
                                   rcrits, 0, pni, pnk)

!-------------------------------------------------------------------------------
!  end  of heterogeneous chemistry
!-------------------------------------------------------------------------------
   end if

   return
end subroutine mach_cam_aerocld
