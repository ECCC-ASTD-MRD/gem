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
! Projet / Project : GEM-MACH version 2
! Fichier / File   : mach_drydep_gas_solver.ftn90
! Creation         : S. Gravel - July 2014, from A. Robichaud - Dec 2002, UPGRADE DEC 2008 (version 2)
! Description      : Calculates aerodynamic resistance for all landuse
!
! Arguments:  IN
!
!    iseason  -> Seasonal categories
!                  1  midsummer with lush vegetation
!                  2  autumn with cropland before harvest
!                  3  later autumn after frost, no snow
!                  4  winter, snow on ground and subfreezing
!                  5  transitional spring with partially green short annuals
!
!   lfu ->       Land Form Use
!
!   metvar2d  -> meteorological variable bloc
!
! Arguments:     OUT
!
! AERO_RESIST -> Aerodynamic resistance for species "species" (s/m)
!
!============================================================================!
!
!!if_on
subroutine mach_gas_drydep_ra2(aero_resist, iseason, lfu, metvar2d)
   use chm_metvar_mod,     only: SIZE_MV2D
   use mach_drydep_mod,    only: lucprm
   use chm_ptopo_grid_mod, only: chm_ni
!!if_off
   use chm_metvar_mod,     only: MV2D_TDIAG, MV2D_WSDIAG, MV2D_DLAT, MV2D_ILMO, MV2D_SNODP
   use chm_utils_mod,      only: chm_lun_out, global_debug
   use chm_consphychm_mod, only: karman, pi, tcdk
   use mach_drydep_mod,    only: zz0

   implicit none
!!if_on
   integer(kind=4), intent (in) :: iseason    (chm_ni)
   real(kind=4),    intent (in) :: metvar2d   (chm_ni, SIZE_MV2D)
!  Land use and land use-dependant aerodynamic resistance - implicit conversion from the 1D vector
   real(kind=4),    intent (in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent(out) :: aero_resist(chm_ni, lucprm)
!!if_off
!
!  Local variables
!
   real(kind=4)            :: ilmob, ylat, rad_to_deg
   real(kind=4)            :: psic, psim, tsurf, href1, vvent
   real(kind=4)            :: psic1, psic2, psim1, psim2
   real(kind=4)            :: ustar, roughlen, ref_hgt, zeta_ref, zeta_zero
   integer(kind=4)         :: i, ii
   integer(kind=4)         :: isea, nlus
   logical(kind=4)         :: local_dbg
! Local constants
   real(kind=4), parameter :: max_ilmob = 100.

!  Code for dry deposition of gaseous species begins here

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   if (local_dbg) then
      write(chm_lun_out, *) ' Begin mach_gas_drydep_ra2'
   end if

   rad_to_deg = 180.0 / pi

   do i = 1, chm_ni
      tsurf = amax1 (metvar2d(i, MV2D_TDIAG) - tcdk, -60.)
      ylat  = metvar2d(i, MV2D_DLAT) * rad_to_deg
!
      isea  = iseason(i)
      if (metvar2d(i, MV2D_SNODP) < .001 .and. isea == 4 .and. tsurf > 5.) then
         isea  = 5    ! transition season
      end if
      if (metvar2d(i, MV2D_SNODP) > 0.005) then
         isea  = 4
      end if
      if (ylat >= 60.0 .and. ylat <= 75.0) then
         if (isea == 4 .and. tsurf >= 5.0) isea = 5
      end if

      ilmob = sign(min(abs(metvar2d(i, MV2D_ILMO)),max_ilmob), metvar2d(i, MV2D_ILMO))

      do nlus = 1, lucprm
         if (lfu(i, nlus) < 0.01) then
            cycle
         end if

         ref_hgt=10.*zz0(nlus, isea)
         if (ref_hgt <= 1.5) ref_hgt=1.5
         if (nlus == 13 .or. nlus == 14) ref_hgt=10.0

! test debug 12 janvier 2012: ref_hgt=10
!        ref_hgt=10.0;

         roughlen = zz0(nlus, isea)

         psic = 0.0
         psim = 0.0

         href1 = amin1((ref_hgt - roughlen),0.)

         if ((ref_hgt * ilmob) > 0.01) then
            psic = -5.0 * (href1 * ilmob)
            psim = psic
         else if ((ref_hgt * ilmob) < -0.01) then
            zeta_ref = max(ref_hgt*ilmob, -1.0)
            zeta_zero = max(max(roughlen*ilmob, -roughlen/ref_hgt), -1.0)
            psic1 = exp(0.598 + 0.39 * log(-zeta_ref) - 0.09 * (log(-zeta_ref))**2)
            psic2 = exp(0.598 + 0.39 * log(-zeta_zero) - 0.09 * (log(-zeta_zero))**2)
            psim1 = exp(0.032 + 0.448 * log(-zeta_ref) - 0.132 * (log(-zeta_ref))**2)
            psim2 = exp(0.032 + 0.448 * log(-zeta_zero) - 0.132 * (log(-zeta_zero))**2)
            psic  = psic1 - psic2
            psim  = psim1 - psim2
         end if

         vvent= metvar2d(i, MV2D_WSDIAG)
         ustar=amax1(karman*vvent/(log(ref_hgt/roughlen) - psim),0.001)

         ! if over ocean or lake, define roughlen (depends on wind over water)
         ! and do 10 iterations on ustar

         if (nlus == 13 .or. nlus == 14) then
            do ii = 1,10
               if (nlus == 13) roughlen = 0.11 * 1.5E-5 / ustar
               if (nlus == 14) roughlen = 0.0144*ustar**2/9.81 + 1.463e-5/(9.1*ustar)
               if ((ref_hgt * ilmob) < -0.01) then
                  zeta_ref = max(ref_hgt*ilmob, -1.0)
                  zeta_zero = max(max(roughlen*ilmob, -roughlen/ref_hgt), -1.0)
                  psim1 = exp(0.032 + 0.448 * log(-zeta_ref) - 0.132 * (log(-zeta_ref))**2)
                  psim2 = exp(0.032 + 0.448 * log(-zeta_zero) - 0.132 * (log(-zeta_zero))**2)
                  psim  = psim1 - psim2
               endif
               ustar=karman * vvent/(log(ref_hgt/roughlen) - psim)
               if ((ref_hgt * ilmob) < -0.01) then
                  zeta_ref = max(ref_hgt*ilmob, -1.0)
                  zeta_zero = max(max(roughlen*ilmob, -roughlen/ref_hgt), -1.0)
                  psic1 = exp(0.598 + 0.39 * log(-zeta_ref) - 0.09 * (log(-zeta_ref))**2)
                  psic2 = exp(0.598 + 0.39 * log(-zeta_zero) - 0.09 * (log(-zeta_zero))**2)
                  psic  = psic1 - psic2
               end if
            end do
         end if

         aero_resist(i, nlus) = (log((ref_hgt + roughlen) / roughlen) - psic) / (karman * (ustar + 0.001))

         if (aero_resist(i, nlus) < 10.0) then
            aero_resist(i, nlus) = 10.0
         elseif (aero_resist(i, nlus) > 9999.0) then
            aero_resist(i, nlus)  =  9999.0
         end if

         enddo  ! loop over all land use categories

      enddo  ! loop over grid point row

   return

end subroutine mach_gas_drydep_ra2
