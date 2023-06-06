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
! Fichier/File   : mach_landuse.ftn90
! Creation       : S. Menard, A. Kallaur  - August 2007
! Description    : Convert into 15 landuse categories the original 26 landuse of GEM.
!                  This is required to do the dry deposition of chemical species.
!
! Extra info     : 1) This mapping is also used in the GRAHM model (Ref: D. Davignon, L.P. Crevier / S. Menard)
!                  2) A bug was found in the mapping and a correction was made in Feb 2012 by A. Robichaud. (see Mantis 2125)
!                  3) Mantis 2680
!                  4) Modified to account for ICE/snow presence on vegetation surfaces
!                     (S. Beagley, S. gravel and D. Akingunola, 2019)
!
!                  26 landuse (original GEM)    --> 15 landuse (GEM-MACH dry deposition)
!                  =====================================================================
!                  1   water                        14  Ocean
!                  2   ice                          12  Ice caps and glaciers
!                  3   inland lake                  13  Inland water
!                  4   evergreen neddleleaf trees    1  Evergreen needleleaf forest
!                  5   evergreen broadleaf trees     2  Evergreen broadleaf forest
!                  6   deciduous neddleleaf trees    2  Evergreen broadleaf forest
!                  7   deciduous broadleaf trees     3  Deciduous needleleaf forest
!                  8   tropical broadleaf trees      4  Deciduous broadleaf forest
!                  9   drought deciduous trees       4  Deciduous broadleaf forest
!                  10  evergreen broadleaf shrub    10  Dwarf trees, shrubs with ground cover (tundra)
!                  11  deciduous shrubs             10  Dwarf trees, shrubs with ground cover (tundra)
!                  12  thorn shrubs                 10  Dwarf trees, shrubs with ground cover (tundra)
!                  13  short grass and forbs        10  Dwarf trees, shrubs with ground cover (tundra)
!                  14  long grass                    6  Grassland
!                  15  crops                         7  Crops, mixed farming
!                  16  rice                          7  Crops, mixed farming
!                  17  sugar                         7  Crops, mixed farming
!                  18  maize                         7  Crops, mixed farming
!                  19  cotton                        7  Crops, mixed farming
!                  20  irrigated crops               7  Crops, mixed farming
!                  21  urban                        15  Urban
!                  22  tundra                        9  Tundra
!                  23  swamp                        11  Wet land with plants
!                  24  desert                        8  Desert
!                  25  mixed wood forests            5  Mixed forest
!                  26  mixed shrubs                 10  Dwarf trees, shrubs with ground cover (tundra)
!
! Arguments:
!            IN
!
!            OUT
!              landuse_out  --> 15 landuse for dry deposition
!
!==============================================================================
!
!!if_on
subroutine mach_landuse(busper, metvar2d, landuse_out)
   use chm_ptopo_grid_mod,  only: chm_ni
   use chm_metvar_mod,      only: SIZE_MV2D
   use mach_drydep_mod,     only: lucprm
!!if_off
   use chm_utils_mod,       only: chm_lun_out, global_debug, ik, chm_error_l
   use chm_metvar_mod,      only: MV2D_GLSEA, MV2D_SNODP, MV2D_SNOF
   use chm_phyvar_mod,      only: vegf, urban
! From rpnphy
   use sfc_options,         only: schmurb
   implicit none
!!if_on
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: metvar2d   (chm_ni, SIZE_MV2D)
   real(kind=4),    intent  (out) :: landuse_out(chm_ni, lucprm)
!!if_off
!
! Local variables
!
   real(kind=4), dimension(chm_ni, 26) :: metvar_vegf
   real(kind=4), dimension(chm_ni) :: lucsum
   real(kind=4)                    :: glsea, snof
   real(kind=4)                    :: qchange, nwluf13, nwluf14, lufsum
   real(kind=4), parameter         :: tolerance = 0.01
   integer(kind=4)                 :: i, j
   logical(kind=4)                 :: local_dbg

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   if (local_dbg) then
      write (chm_lun_out, *) 'Remapping of 26 landuse types to 15 landuse categories'
   end if
!
! Load vegetation fraction from permanent bus
!
   do j = 1, 26
      do i = 1, chm_ni
         metvar_vegf(i, j) = busper(vegf + ik(i, j, chm_ni))
      end do
   end do

  ! Update the URBAN fraction if running with 'TEB'
   if (schmurb == 'TEB') then
      do i = 1, chm_ni
         metvar_vegf(i, 21) = busper(urban + i - 1)
      end do
   end if

   do i = 1, chm_ni
      landuse_out(i,  1) = metvar_vegf(i,  4)
      landuse_out(i,  2) = metvar_vegf(i,  5) + metvar_vegf(i,  8)
      landuse_out(i,  3) = metvar_vegf(i,  6)
      landuse_out(i,  4) = metvar_vegf(i,  7) + metvar_vegf(i,  9)
      landuse_out(i,  5) = metvar_vegf(i, 25)
      landuse_out(i,  6) = metvar_vegf(i, 14)
      landuse_out(i,  7) = metvar_vegf(i, 15) + metvar_vegf(i, 16) + &
                           metvar_vegf(i, 17) + metvar_vegf(i, 18) + &
                           metvar_vegf(i, 19) + metvar_vegf(i, 20)
      landuse_out(i,  8) = metvar_vegf(i, 24)
      landuse_out(i,  9) = metvar_vegf(i, 22)
      landuse_out(i, 10) = metvar_vegf(i, 10) + metvar_vegf(i, 11) + &
                           metvar_vegf(i, 12) + metvar_vegf(i, 13) + &
                           metvar_vegf(i, 26)
      landuse_out(i, 11) = metvar_vegf(i, 23)
      landuse_out(i, 12) = metvar_vegf(i,  2)
      landuse_out(i, 13) = metvar_vegf(i,  3)
      landuse_out(i, 14) = metvar_vegf(i,  1)
      landuse_out(i, 15) = metvar_vegf(i, 21)
      lucsum(i) = 0.0
   end do
!
! Land use corrections
! NOTE: The variable *glsea* represents the fraction of the **water**
!       (inland water or sea) that is covered by ice. It does not represent
!       the fraction of the **grid box** that is covered by ice.
   do i = 1, chm_ni
! PART 1:  If sea-ice present we reassign ice covered part to "ice cap"
      glsea = min(1.0, max(metvar2d(i, MV2D_GLSEA), 0.0))
      if (glsea > 0.0) then ! some ice covered water
         ! reduce sea and inland water fraction by ice fraction
         nwluf13 = landuse_out(i, 13) * (1.0 - glsea)
         nwluf14 = landuse_out(i, 14) * (1.0 - glsea)
         qchange = (landuse_out(i, 13) - nwluf13) + &
                   (landuse_out(i, 14) - nwluf14)
         landuse_out(i, 12) = min(landuse_out(i, 12) + qchange, 1.0)
         landuse_out(i, 13) = nwluf13
         landuse_out(i, 14) = nwluf14
      end if
!
! PART 2: If snow present we ignore real land use category and impose ICE conditions
      snof = metvar2d(i, MV2D_SNOF)
      if (metvar2d(i, MV2D_SNODP) > 0.005 .and. snof > 0.0) then
! Due to snow presence impose ICE deposition conditions
! Ice/snow fre water surfaces dealt with above, adjust other classes except
! water and where ice which is already present.
         lufsum = 0.0
         do j = 1, lucprm
            if (j /= 12 .and. j /= 13 .and. j /= 14) then
! open water remains unchanged, ice adds original below.
               landuse_out(i, 12) = landuse_out(i, 12) + &
                                    landuse_out(i, j) * snof
               ! reduce landuse by snow-covered fraction
               landuse_out(i, j) = landuse_out(i, j) * (1.0 - snof)
            end if
         end do
!
      end if

   end do
!
! Look if the summation of the 15 landuse_out types equal 1.
   do j = 1, lucprm
      do i = 1, chm_ni
         lucsum(i) = lucsum(i) + landuse_out(i, j)
      end do
   end do
   if (local_dbg) then
      write (chm_lun_out, *) 'Land use sum', lucsum
   end if

   do i = 1, chm_ni
      if (abs(1.0 - lucsum(i)) > tolerance) then
         write(0, *) '### Error in mach_landuse ###'
         write(0, *) '# Summation of vegetation fractions over all landuse ', &
                    & 'is not equal to 1: lucsum(', i, ') = ', lucsum(i)
         write(0, *) '###         ABORT         ###'
         chm_error_l = .true.
         return
      end if
   end do

   return

end subroutine mach_landuse
