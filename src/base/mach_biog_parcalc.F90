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
! Fichier/File   : mach_biog_parcalc.ftn90
! Creation       : H. Landry, S. Menard, A. Kallaur - mai 2006
! Description    : Compute the Photosynthetically Active Radiation (PAR)
!
! Extra info     : (from the original file in AURAMS v1.3.1 modele/pre/parcalc.ftn90)
!                  Computes the radiation terms needed for biogenic calculations.
!                  two terms are computed: total solar radiation (langleys/min)
!                  and par (ue/m**2-s).  the methodology for this calculation
!                  was taken from iqbal, m., 1983. an introduction to solar
!                  radiation. academic press, new york, pp. 202-210.
!
! Arguments: IN
!              cosine_solar_p -> cosine of solar angle
!              metvar2d       -> local storage (chemistry lib only, see chm_exe) for met. fields 
!                                copied from Physics buses (see chm_load_metvar).
!            OUT
!              pardb  -> PAR direct beam
!              pardif -> PAR diffuse
!
!==============================================================================
!
!!if_on
subroutine mach_biog_parcalc(cosine_solar_p, PARDB, PARDIF, metvar2d)
   use chm_ptopo_grid_mod, only: chm_ni
   use chm_metvar_mod,     only: SIZE_MV2D
!!if_off
   use chm_metvar_mod,     only: MV2D_PPLUS, MV2D_FLUSOLIS
   use chm_utils_mod,      only: chm_error_l
   implicit none
!!if_on
   real(kind=4), intent (in) :: cosine_solar_p(chm_ni)           ! cosine of solar angle
   real(kind=4), intent(out) :: pardb         (chm_ni)           ! PAR direct beam
   real(kind=4), intent(out) :: pardif        (chm_ni)           ! PAR diffuse!
   real(kind=4), intent (in) :: metvar2d      (chm_ni,SIZE_MV2D) ! Met. variables
!!if_off
!
!  Local variables
!
   real(kind=4)            :: ot, wa, ratio
   real(kind=4)            :: rad_direct_vis, rad_diffuse_vis, rad_direct_ir, rad_diffuse_ir
   real(kind=4)            :: rad_visible_total, rad_ir_total
   real(kind=4)            :: fract_rad_visible, fract_direct_vis, fract_diffuse_vis
!  * Loop variable
   integer(kind=4)         :: i
!  * some constants
   real(kind=4), parameter :: min_cosine   = 0.0208   ! Minimum cosine of solar angle allowed for calculation ( about cos(88.8 deg) )
   real(kind=4), parameter :: min_pressure = 10000.0  ! Minimum surface pressure allowed for calculation ( in pascal )
   real(kind=4), parameter :: std_pressure = 101325.0 ! Standard pressure at sea level ( in pascal )
   real(kind=4), parameter :: watt_2_umol  = 4.6      ! conversion factor from watt to umol

   pardb  = 0.0
   pardif = 0.0

!    Loop on the slab
   do i = 1, chm_ni
!     First do some value check
      if (metvar2d(i, MV2D_PPLUS) < min_pressure) then
         write(0, *) 'ERROR, surface pressure is too low at metvar2d(i, MV2D_PPLUS)(i, MV2D_PPLUS) = ', metvar2d(i, MV2D_PPLUS)
         write(0, *) '*** ABORTING ***'
         chm_error_l = .true.
         return
      endif

      if (cosine_solar_p(i) >= min_cosine .AND. metvar2d(i, MV2D_FLUSOLIS) > 0.0) then
!        Compute clear sky radiation terms (maximum)
         ot = metvar2d(i, MV2D_PPLUS) / std_pressure / cosine_solar_p(i)
         rad_direct_vis  = 600.0 * exp(-0.185 * ot)         * cosine_solar_p(i)
         rad_diffuse_vis = 0.42  * (600.0 - rad_direct_vis) * cosine_solar_p(i)
         wa = 1320.0 * 0.077 * (2.0 * ot)**0.3

         rad_direct_ir  =        (720.0 * exp(-0.06 * ot)- wa) * cosine_solar_p(i)
         rad_diffuse_ir = 0.65 * (720.0 - wa - rad_direct_ir)  * cosine_solar_p(i)

         rad_visible_total = rad_direct_vis + rad_diffuse_vis
         rad_ir_total      = rad_direct_ir  + rad_diffuse_ir
         fract_rad_visible = rad_visible_total / (rad_ir_total + rad_visible_total)

         ratio = metvar2d(i, MV2D_FLUSOLIS) / (rad_ir_total + rad_visible_total)

!        Compute fraction of visible that is direct beam

         if (ratio >= 0.89) then
            fract_direct_vis = rad_direct_vis / rad_visible_total * 0.941124
         else if (ratio < 0.21) then
            fract_direct_vis = rad_direct_vis / rad_visible_total * 9.55E-3
         else
            fract_direct_vis = rad_direct_vis / rad_visible_total * (1.0 - ((0.9 - ratio) / 0.7)**0.666667)
         endif

         fract_diffuse_vis = 1.0 - fract_direct_vis

!     Compute PAR, direct beam and diffuse, i umol/m2.sec

         pardb(i)  = metvar2d(i, MV2D_FLUSOLIS) * fract_rad_visible * fract_direct_vis  * watt_2_umol
         pardif(i) = metvar2d(i, MV2D_FLUSOLIS) * fract_rad_visible * fract_diffuse_vis * watt_2_umol

      endif

   enddo  ! end of slab loop

end subroutine mach_biog_parcalc
