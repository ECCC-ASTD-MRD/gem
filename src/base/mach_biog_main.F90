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

!
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
! Fichier/File   : mach_biog_main.ftn90
! Creation       : H. Landry, S. Menard, A. Kallaur - mai 2006
! Description    : This is the main biogenic subroutine.  It drives the
!                  the computation and the output of the biogenic emissions
!
! Extra info     : From here, there should not be any mention of the source of the
!                  data: we don't want to carry all the bus interface with us.
!                  This function is detached from GEM and CHEM.
!                  It is also assumed that all arrays are 2D of size ni
!                  PAR here stands for Photosynthetically Active Radiation.
! Extra Info     : Added VOC emission into canopy shaded levels
!                  (P. Makar and D. Akingunola, Fall 2020)
!
! Arguments:  IN
!                metvar2d       --> local storage (chemistry lib only, see chm_exe) for met. fields
!                metvar3d           copied from Physics buses (see chm_load_metvar).
!                iseasons       --> seasons field (5 seasons defined)
!                emisbio_can    --> 3D Biogenic emissions within canopy columns
!
!             IN/OUT
!	         busper         --> permanent bus
!
!==============================================================================
!
!!if_on
subroutine mach_biog_main(busper, metvar2d, metvar3d, iseasons,  metvar3dcan, &
                          emisbio_can, kcan, ni_can)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
   use mach_pkg_gas_mod,     only: num_be_sp
!!if_off
   use chm_metvar_mod,       only: MV2D_CANG, MV3D_TPLUS, MV3D_ZMOM
   use chm_utils_mod,        only: chm_error_l, CHM_MSG_DEBUG
   use chm_nml_mod,          only: chm_timings_L, chm_canopy_shading_l, &
                                   chm_sat_seasons_l, chm_biog_s
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_be_std, sp_LAI, sp_STSE, sp_FRT,  &
                                   sp_CRL, sp_FRL, sp_HC
   use mach_pkg_gas_mod,     only: num_be_std, be_species
   use mach_headers_mod,     only: mach_biog_parcalc, mach_biog_beis
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: ni_can
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: metvar2d   (chm_ni,SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d   (chm_ni,chm_nk,SIZE_MV3D)
   real(kind=4),    intent   (in) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   integer(kind=4), intent   (in) :: iseasons   (chm_ni)
   integer(kind=4), intent   (in) :: kcan       (ni_can, nkc)
   real(kind=4),    intent  (out) :: emisbio_can(ni_can, nkt, num_be_sp)
!!if_off
!
!   Local variables
!
!  * Loop variable
   integer(kind=4) :: i, ic, kc, nn, nsp, nlev
!
!   *  Constants
   real(kind=4), parameter :: seasons_coef(5) = (/ 1.0, 0.66, 0.25, 0.0, 0.5 /)
!
   real(kind=4), dimension(num_be_sp, chm_ni)  :: emisbio
!     *  Standard emission (without correction)
   real(kind=4), dimension(num_be_std, chm_ni) :: summer_std_p, winter_std_p
!     * Leaf Area Index
   real(kind=4) :: lai_p
!     * meteorological variables needed for some calculation
   real(kind=4), dimension(chm_ni) :: cosine_solar_p   ! cosine of solar angle
   real(kind=4), dimension(nkc + 1, chm_ni) :: air_temp  ! Air temperature
   real(kind=4), dimension(ni_can, nkt)     :: zmomcan   ! Momentum height
!
   real(kind=4), dimension(chm_ni) :: sesn_coef
   real(kind=4) :: kbe, canparscat, canpardif, parshade, fracshade
   real(kind=4) :: par, parsun, laisun, fracsun, tmpsun, tmpshade
   real(kind=4), dimension(nkc + 1, chm_ni) :: light_correction
!  * Photosynthetically Active Radiation variables
   real(kind=4), dimension (chm_ni) :: pardb, pardif  ! PAR direct beam, PAR diffuse
!
!  * Canopy shading  variables
   real(kind=4),    dimension(nkc + 1)         :: lai_acc
   real(kind=4),    dimension(nkc + 1, chm_ni) :: laifrac
   real(kind=4),    dimension(ni_can)          :: canopy_hgt
   logical(kind=4), dimension(chm_ni)          :: canopy_column
!
!  External subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp
!
   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_biog [BEGIN]')
   if (chm_timings_L) call timing_start_omp(320, 'mach_biog', 480)
!
   do i = 1, chm_ni
      nn = sp_be_std
      do nsp = 1, num_be_std
         summer_std_p(nsp, i) = busper(sm(nn)     % per_offset + i - 1)
         winter_std_p(nsp, i) = busper(sm(nn + 1) % per_offset + i - 1)
         nn = nn + 2
      end do
!  * Ensure that there're no negative values in  the cosine of solar zenith (cang)
      cosine_solar_p(i) = max(metvar2d(i, MV2D_CANG), 0.0)
   end do
!
!  * Computation of the Photosynthetically Active Radiation (PARcalc)
   call mach_biog_parcalc(cosine_solar_p, pardb, pardif, metvar2d)
   if (chm_error_l) return
!
!  * Seasonal coefficients
   if (chm_sat_seasons_l) then
!  sesn_coef based on proximity to maximum monthly LAI:
      do i = 1, chm_ni
         sesn_coef(i) = busper(sm(sp_STSE) % per_offset + i - 1)
      end do
   else
!  sesn_coef based on spatial mask inherited from ADOM model (?)
      do i = 1, chm_ni
         sesn_coef(i) = seasons_coef(iseasons(i))
      end do
   end if

   if (chm_canopy_shading_l) then
      do i = 1, chm_ni
         canopy_column(i) = busper(sm(sp_FRT) % per_offset + i - 1) > 0.0
      end do
   else
      canopy_column = .false.
      canopy_hgt    = 0.0
   end if
!
!* -----------------------------------------------------------------------------
!*    Light correction
!*       The formula used is (ref. Geron et al):
!*
!*          cl(i, j, k) = (a*cl1*L(k)/(1+a*a*L(k)*L(k))**0.5
!*       where:
!*             L(k) : PAR value for leaf area index type "k"
!* -----------------------------------------------------------------------------
!
   light_correction = 0.0    ! Initialize
   ic = 0
   do i = 1, chm_ni
      lai_p = busper(sm(sp_LAI) % per_offset + i - 1)
!
      if (canopy_column(i)) then
! Fraction of LAI at heights of 1.0, 0.5 and 0.2 canopy height ,
! counting down from top:
         laifrac(1, i) = busper(sm(sp_FRL) % per_offset + i - 1)
         laifrac(2, i) = busper(sm(sp_FRL + 1) % per_offset + i - 1)
         laifrac(3, i) = busper(sm(sp_FRL + 2) % per_offset + i - 1)
         laifrac(4, i) = 1.0
!
! Fraction of cumulative LAI above each canopy layer interface,
! counting down from top: sm(sp_CRL), sm(sp_CRL + 1), sm(sp_CRL + 2)

! Accumulated LAI for determining amount of LAI reaching leaves: from values
! at 0.75 hc, 0.5 hc, 0.35 hc, and 0.20 hc calculated earlier
!  value of lai at 0.75 hc used to represent accumulated LAI within top canopy
!  layer for isoprene PAR correction:
         lai_acc(1) = busper(sm(sp_CRL) % per_offset + i - 1) * lai_p
!  value of lai at 0.5 hc is used to represent accumulated LAI within middle
!  canopy layer for isoprene PAR correction:
         lai_acc(2) = busper(sm(sp_CRL + 1) % per_offset + i - 1) * lai_p
!  value of lai at 0.2 hc is used to represent accumulated LAI within bottom
!  layer for isoprene PAR correction:
         lai_acc(3) = busper(sm(sp_CRL + 2) % per_offset + i - 1) * lai_p
!
         ic = ic + 1
         do kc = 1, nkc
            air_temp(kc, i) = metvar3dcan(ic, kcan(ic, kc), MV3D_TPLUS)
         end do
         canopy_hgt(ic) = busper(sm(sp_HC) % per_offset + i - 1)
         zmomcan(ic, :) = metvar3dcan(ic, :, MV3D_ZMOM)
!
         nlev = nkc + 1
      else
         laifrac(:, i) = 1.0
         nlev = 1
      end if
!
!  Generate the light correction at the middle of each layer:
!  lai_p is the total lai within the layer
      lai_acc(nlev) = lai_p
      if (lai_p > 0.0) then
         par = pardb(i) + pardif(i)
         if ((par > 0.01) .and. (lai_p <= 0.1) .and. &
            (cosine_solar_p(i) > 0.02079483)) then

            light_correction(1:nlev, i) = get_light_correction(min(par, 2550.0))

         else
            kbe = 0.5 * sqrt(1. + (tan(acos(cosine_solar_p(i))))**2)
!  Light correction is the value for the accumulated canopy down
!  to the middle of each canopy layer in this version.
            do kc = 1, nlev
               canparscat = 0.5 * pardb(i) * (exp(-0.894 * kbe * lai_acc(kc)) - &
                           exp(-1.0 * kbe * lai_acc(kc)))
               canpardif = pardif(i) * (1.0 - exp(-0.6 * lai_acc(kc))) / &
                           (0.61 * lai_acc(kc))
               parshade  = canpardif + canparscat
               parsun    = kbe * par + parshade
               laisun    = (1.0 - exp(-kbe * lai_acc(kc))) / kbe
               fracsun   = laisun / lai_acc(kc)
               fracshade = 1.0 - fracsun

               if (parsun <= 0.01) then
                  tmpsun = 0.0
               else
                  tmpsun = fracsun * get_light_correction(parsun)
               end if
               if (parshade <= 0.01) then
                  tmpshade = 0.0
               else
                  tmpshade = fracshade * get_light_correction(parshade)
               end if
               light_correction(kc, i) = tmpsun + tmpshade
            end do
         end if
      end if
!
      air_temp(nlev, i) = metvar3d(i, chm_nk, MV3D_TPLUS)
   end do
!
! Compute the biogenic emissions
!
   if (chm_biog_s(1:4) == 'BEIS') then
      call mach_biog_beis(emisbio, emisbio_can, summer_std_p, winter_std_p, &
                          sesn_coef, air_temp, zmomcan, light_correction,   &
                          laifrac, canopy_hgt, canopy_column, ni_can)
      if (chm_error_l) return
   end if
!
! Return the biogenic emissions to the bus
   do i = 1, chm_ni
      do nsp = 1, num_be_sp
         nn = be_species(nsp)
         busper(sm(nn) % be_offset + i - 1) = emisbio(nsp, i)
      end do
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_biog [END]')
   if (chm_timings_L) call timing_stop_omp(320)
   !-----------------------------------------------------------------

   return

   contains
!
!  A simple function declaration for a repetitive calculation
!
   real function get_light_correction(PAR_flux)
      real PAR_flux
      get_light_correction = (0.0028782 * PAR_flux) / sqrt (1.0 + 7.29e-06 * (PAR_flux**2))
   end function get_light_correction

end subroutine mach_biog_main
