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
! Fichier/File   : mach_mie_opt.ftn90
! Creation       : P.A. Makar, W. Gong, and D. Akingunola
!
! Description    : Estimate aerosol optical properties optical depth, asymmetry
!                  factor, and the single-scattering albedo
!
! Arguments:  IN
!               chem_tr -> Chemical tracers concentations (ug/kg)
!               trwtrow -> Aerosol liquid water content for each bin
!               rhoa    -> Air density (kg/m3)
!               thlev   -> Layer thickness [m]
!
!             IN/OUT
!               busper  -> Permanent bus for chemistry
!               busvol  -> Volatile bus for chemistry
!
!============================================================================
!
!!if_on
subroutine mach_mie_opt(busper, busvol, aero_tr, trwtrow, thlev, rhoa, pni, pnk,&
                        nmod, kmod)
   use chm_ptopo_grid_mod,   only: chm_nk
   use mach_cam_utils_mod,   only: isize, ntr
!!if_off
   use chm_utils_mod,        only: ik
   use chm_nml_mod,          only: nk_start_pm
   use chm_ptopo_grid_mod,   only: chm_ni
   use chm_species_idx_mod,  only: sp_OPTD, sp_ASYM, sp_SSCA
   use chm_species_info_mod, only: species_master
   use mach_cam_headers_mod, only: mach_cam_main, mach_cam_sfss
   use mach_cam_utils_mod,   only: icom, iae1, rhop0, pvol, binrange
   use mach_mie_data_mod,    only: mach_mie_lookup, nwl_aod
   implicit none

!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   integer(kind=4), intent   (in) :: nmod   (pni)
   integer(kind=4), intent   (in) :: kmod   (pni, chm_nk)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: aero_tr(pni, pnk, ntr)
   real(kind=4),    intent   (in) :: trwtrow(pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev  (pni, pnk)
   real(kind=4),    intent   (in) :: rhoa   (pni, pnk)
!!if_off
!
!  Local variables
!
   integer(kind=4)            :: ii, ix, kk, kz, mm, iae, nn, nt
   integer(kind=4)            :: opt_id, ssc_id, asy_id, this_ik
   real(kind=4)               :: rd, aerotmp, wattmp, rwmax
   real(kind=4)               :: aerovol(pni, pnk, isize)
   real(kind=4)               :: aerowat_vol(pni, pnk, isize)
   real(kind=4)               :: numconc(pni, pnk, isize)
! WFRAC   - water mass fraction for determining refractive
!           index (interpolating between dry and pure water particles)
   real(kind=4)               :: wfrac  (pni, pnk, isize)
! RW      - wet radius for use in Mie code
   real(kind=4)               :: rw     (pni, pnk, isize)

   real(kind=4), dimension(nwl_aod, pni, pnk) :: bext, ssca, asym
   real(kind=4), dimension(pni, pnk, nwl_aod) :: totopt, sscalb, asymbn
   real(kind=4), dimension(pni, nwl_aod)      :: zoptd
!  mfac is an upper limit on the volume ratio of wet to dry particles.
!  The value of mfac = 102000 roughly corresponds to bext(wet)/bext(dry) = 26,
!  which is the value assumed at 99% RH for the IMPROVE formula, and
!  corresponds to a wet to dry particle mass ratio of about 9.7E5.
   real(kind=4),   parameter :: mfac = 102000.0
!  for unit convertion
   real(kind=4),   parameter :: ug2kg = 1.E-9
!

   do mm = 1, isize
      do kk = 1, pnk
         do ii = 1, pni
! Evaluate the volume of aerosol water / m3 of air from the aerosol-bound water
            aerowat_vol(ii, kk, mm) = trwtrow(ii, kk, mm) * ug2kg * &
                                      rhoa(ii, kk) * 1.0E-3
!
!  Evaluate the total dry particle volume for the given size from the chemistry
            aerovol(ii, kk, mm) = 0.0
            do iae = 1, icom
               nt = mm + isize * (iae - 1) + (iae1 - 1)
               aerovol(ii, kk, mm) = aerovol(ii, kk, mm) + &
                                     (aero_tr(ii, kk, nt) * rhoa(ii, kk) / &
                                      rhop0(iae))
            end do
!
!Number concentration from dry particle volume / m3 air / single particle volume
            numconc(ii, kk, mm) = aerovol(ii, kk, mm) / pvol(mm)
         end do
      end do
   end do
!
! Creates water fraction, number concentration and wet radius arrays for
! 12 bin distribution for later use in bhmie.
!
! For each bin, determine the mass and wet particle size. The latter is
! truncated if the aerosols are close to deliquescence, since the Mie code
! will have trouble converging for very large particle sizes i.e. very large
! values of the Mie parameter.
   do mm = 1, isize
!
!  Dry particle radius:
      rd = (binrange(1, mm) + binrange(2, mm)) * 0.5

!  Determine the wet particle size and volume:
      do kk = 1, pnk
         do ii = 1, pni
            aerotmp = aerovol(ii, kk, mm)
            wattmp  = aerowat_vol(ii, kk, mm)
! Wet particle radius default is dry particle radius.
            rw(ii, kk, mm) = rd
! A limiting value of 100 is used here for consistency with the lower limit
! in the Mie lookup table application, and to avoid going below the low number
! limit in adding water to the smaller sized aerosols.
            wfrac(ii, kk, mm) = 0.0
            if (numconc(ii, kk, mm) > 100.0) then
!
! Maximum allowed wet radius is from (volume of single particle * mfac /(4/3 pi ) )**1/3
               rwmax = (mfac * aerotmp / (numconc(ii, kk, mm) * 4.189))**0.33333333
!
!  Calculate wet particle radius:  add volume of aerosol water / m3 air to volume of
!  dry aerosol / m3 air, divide by number of particles / m3 air, divide by
!  (4/3 pi ), take the cube root.
               rw(ii, kk, mm) = max(rd, ((aerotmp + wattmp) / &
                                     (numconc(ii, kk, mm) * 4.189))**0.3333333)
!
               if (rw(ii, kk, mm) > rwmax) then
! add 4/3 pi * number of particles/m3 air * (rwmax^3 -rd^3) to wet particle volume
                  wattmp = 4.189 * numconc(ii, kk, mm) * &
                                max((rwmax * rwmax * rwmax - rd * rd * rd), 0.0)
                  rw(ii, kk, mm) = rwmax
               end if
!
!  add kg H2O/kg * m3 H2O/kg H2O * kg air/m3 air (= m3 H2O / m3 air) to wfrac
               wfrac(ii, kk, mm) = wfrac(ii, kk, mm) + wattmp
!
!  add ditto to total particle volume:
!
               aerotmp = aerotmp + wattmp
!
!  wfrac(i,l,is) now holds the volume of water for bin size "is", and aerotmp
!  holds the total particle volume. Determine the water volume fraction:
               wfrac(ii, kk, mm) = wfrac(ii, kk, mm) / aerotmp
            end if

         end do
      end do
   end do
!
! Evaluate the aerosol optical properties at specified wavelengths.
! NOTE: These 4 wavelengths are mid-band values required for the GEM radiative
! transfer feedbacks.  The remaining values (if any) are for diagnostic
! output calculations.  Units of wavelength here are metres.
! lam_for_cccmarad = (/0.4448E-06, 0.9401E-06, 0.17857E-05, 0.31905E-05/)
!
   call mach_mie_lookup(wfrac, rw, numconc, bext, ssca, asym, pni, pnk)

   do nn = 1, nwl_aod
      do kk = 1, pnk
         do ii = 1, pni
! note that BEXT from the code is in m**-1
! The limits imposed below are to prevent division by zero errors in s/r raddriv.
!
            totopt(ii, kk, nn) = bext(nn, ii, kk) * thlev(ii, kk)
!
            asymbn(ii, kk, nn) = asym(nn, ii, kk)
            if (totopt(ii, kk, nn) > 1.0e-10) then
               sscalb(ii, kk, nn) = max(min(1.0, ssca(nn, ii, kk)), 1.0e-5)
!               asymb(ii, kk, nn) = max(min(1.0, asym(nn, ii, kk)), 1.0e-5)
            else
               totopt(ii, kk, nn) = 1.0e-10
               sscalb(ii, kk, nn) = 1.0
!               asymb(ii, kk, nn) = 1.0
            end if
         end do
      end do
   end do
!
! *  Place the (3-D) aerosol optical depth, single scattering albedo and
! *  asymmetry factor into the permanent bus for storage and output:
   do nn = 1, nwl_aod
      opt_id = nn + sp_OPTD - 1
      ssc_id = nn + sp_SSCA - 1
      asy_id = nn + sp_ASYM - 1
      do kz = nk_start_pm, chm_nk
         do ii = 1, pni
            kk = kmod(ii, kz) - nk_start_pm + 1
            ix = nmod(ii)
            this_ik = ik(ix, kz, chm_ni)
            busper(species_master(opt_id) % per_offset + this_ik) = totopt(ii, kk, nn)
            busper(species_master(ssc_id) % per_offset + this_ik) = sscalb(ii, kk, nn)
            busper(species_master(asy_id) % per_offset + this_ik) = asymbn(ii, kk, nn)
         end do
      end do
!
! Save column (2-D) total aerosol optical depth to output:
      if (species_master(opt_id) % out_offset > 0) then
         do ii = 1, pni
            ix = nmod(ii) - 1
            zoptd(ii, nn) = 0.0
            do kk = 1, pnk
               zoptd(ii, nn) = totopt(ii, kk, nn) + zoptd(ii, nn)
            end do
            busvol(species_master(opt_id) % out_offset + ix) = zoptd(ii, nn)
         end do
      end if
   end do

   return
end subroutine mach_mie_opt
