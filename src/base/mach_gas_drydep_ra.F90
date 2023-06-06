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
! Projet / Project : GEM-MACH
! Fichier / File   : mach_gas_drydep_ra.ftn90
! Creation         : S. Gravel - Apr 2013 from A. Robichaud - Dec 2002
! Description      : Calculates aerodynamic resistance for all landuse
!
!
! Arguments:  IN
!
!    season   -> Seasonal categories
!                  1  midsummer with lush vegetation
!                  2  autumn with cropland before harvest
!                  3  later autumn after frost, no snow
!                  4  winter, snow on ground and subfreezing
!                  5  transitional spring with partially green short annuals
!
!   lfu       -> Land Form Use
!   metvar2d  -> meteorological variable bloc
!
! Arguments:     OUT
!
! AERO_RESIST -> Aerodynamic resistance for species "species" (s/m)
!
!
!============================================================================!
!
!!if_on
subroutine mach_gas_drydep_ra(aero_resist, iseason, lfu, metvar2d)
   use mach_drydep_mod,    only: lucprm
   use chm_metvar_mod,     only: SIZE_MV2D
   use chm_ptopo_grid_mod, only: chm_ni
!!if_off
   use chm_metvar_mod,     only: MV2D_UE, MV2D_ILMO, MV2D_SNODP
   use chm_utils_mod,      only: chm_lun_out, global_debug
   use chm_consphychm_mod, only: karman
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
   real(kind=4)            :: ilmob
   real(kind=4)            :: psic, href1
   real(kind=4)            :: psic1, psic2
   real(kind=4)            :: ustar, roughlen, zeta_ref, zeta_zero
   integer(kind=4)         :: i
   integer(kind=4)         :: isea, nlus
! Local constants
   real(kind=4), parameter :: ref_hgt = 10.0
   real(kind=4), parameter :: max_ilmob = 100.
! Debug logical switch
   logical(kind=4)         :: local_dbg

!  Code for aerodynamic resistance of various landuse categories

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   if (local_dbg) then
      write(chm_lun_out, *) ' Begin mach_gas_drydep_ra:'
   end if

  !  Initialize
   aero_resist = 0.0

! Begin main loop over Domain

   do i = 1, chm_ni

      ustar = max(metvar2d(i, MV2D_UE), 0.001)

      ! put a filter on Inverse Monin-Obukhov length from GEM
      ilmob = sign(min(abs(metvar2d(i, MV2D_ILMO)),max_ilmob), metvar2d(i, MV2D_ILMO))

      isea     = iseason(i)

      ! Check for snow depth at particular grid point, and
      ! override Season value if there is.
      if (metvar2d(i, MV2D_SNODP)  > 0.005) then
         isea  = 4
      end if

      ! loop over all LUC (Land Use Categories)

      do nlus = 1, lucprm

         if (lfu(i, nlus) < 0.01) then
            cycle
         end if

         ! Take roughness length from look-up table
         roughlen = zz0(nlus, isea)

         psic = 0.0

         ! following discussion between A. ROBICHAUD & Y. DELAGE re: GEM (fall 2002), psic = psic'(href) - psic'(z0)
         href1 = ref_hgt - roughlen

         if ((ref_hgt * ilmob) > 0.01) then
            psic = -5.0 * (href1 * ilmob)
         else if ((ref_hgt * ilmob) < -0.01) then
            zeta_ref = max(ref_hgt*ilmob, -1.0)
            zeta_zero = max(max(roughlen*ilmob, -roughlen/ref_hgt), -1.0)
            psic1 = exp(0.598 + 0.39 * log(-zeta_ref) - 0.09 * (log(-zeta_ref))**2)
            psic2 = exp(0.598 + 0.39 * log(-zeta_zero) - 0.09 * (log(-zeta_zero))**2)
            psic  = psic1 - psic2
         end if

         ! if over ocean or lake, define roughlen (depends on wind over water)
         if ((nlus == 13) .or. (nlus == 14)) then
            roughlen = 0.11 * 1.5E-5 / ustar
         end if

         aero_resist(i, nlus) = (log((ref_hgt + roughlen) / roughlen) - psic) / (karman * (ustar + 0.001))

         if (aero_resist(i, nlus) < 0.0) then
            aero_resist(i, nlus) = 0.0
         elseif (aero_resist(i, nlus) > 9999.0) then
            aero_resist(i, nlus)  =  9999.0
         end if


      end do  ! loop over all land use categories

   end do  ! loop over grid point row

   if (local_dbg) then
      write(chm_lun_out, *) 'mach_gas_drydep_ra EXIT'
   end if

   return

end subroutine mach_gas_drydep_ra
