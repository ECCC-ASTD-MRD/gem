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
! Fichier/File   : mach_cam_drydep_main.ftn90
! Creation       : S. Gong, P. Makar, L. Zhang, S. Gravel and B. Pabla for GEM-MACH, June 2008
!
! Description    : This is an entry subroutine of the dry deposition processes
!
! Extra info     : - First version created by S. Gong Aug 04 1994 for CAM
!                    Method:
!                       parameterization of dry deposition scheme by Filippo Giorgi
!                       [JGR, 91(D9), 1986] was used.
!                       Two layer model.
!                  - Vectorized the whole program and add working spaces.
!                    fveg, cd012, cvocdl, aexp (S. Gong, Dec 19, 1996)
!                  - Refined parameterization for particles compared with empirical
!                    model resutls. (L. Zhang, Jan 15, 1998)
!                  - Pass all control variables for 15 land-use from the
!                    control file. (S. Gong, Aug 01, 1999)
!
! Arguments:  IN
!               iseasn  -> Assigned season descriptors
!               Ra      -> Aerodynamic resistance for species "species" (s/m)
!               usi     -> friction velocity
!               thlev   -> Layer thickness [m]
!             roarow    -> Air density (kg/m3)
!              rhsize   -> Wet radius
!              fland    -> Landuse
!               pdiff   -> diffusion coefficient
!               amu     -> Air's dynamic viscosity
!
!             OUT
!              DRYFLX   -> Dry deposition flux
!
!             IN/OUT
!              XROW     -> tracer concentration
!              PDEPV    -> gravitational settling velocity
!
!============================================================================
!
!!if_on
 subroutine mach_cam_drydep_main(iseasn, ra, usi, thlev, roarow, rhsize, fland, &
                                 xrow, amu, pdepv, pdiff, dryflx, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr, icom
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_utils_mod,        only: CHM_MSG_DEBUG
   use chm_nml_mod,          only: chm_timings_L, chm_pm_drydep_s
   use mach_cam_headers_mod, only: mach_cam_drypar, mach_cam_drydep1, &
                                   mach_cam_drydep2

   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   integer(kind=4), intent   (in) :: iseasn (pni)
   real(kind=4),    intent   (in) :: usi   (pni)
   real(kind=4),    intent   (in) :: ra    (pni, lucprm)
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: fland (pni, lucprm)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent(inout) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent  (out) :: dryflx(pni, icom)
!!if_off
   external msg_toall, timing_start_omp, timing_stop_omp

!
   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_cam_drydepo [BEGIN]')
   if (chm_timings_L) call timing_start_omp(377, 'mach_cam_drydepo', 340)

   call mach_cam_drypar(thlev, roarow, pdiff, rhsize, amu, fland, &
                        xrow, pdepv, iseasn, ra, usi, pni, pnk)
!
!  vertical mass balance.
!
   if (chm_pm_drydep_s == 'ZHANG_MAKAR') then
      call mach_cam_drydep2(xrow, dryflx, pdepv, thlev, roarow, pni, pnk)
   else !(chm_pm_drydep_s == 'ZHANG')
      call mach_cam_drydep1(xrow, dryflx, pdepv, thlev, roarow, pni, pnk)
   end if

   call msg_toall(CHM_MSG_DEBUG, 'mach_cam_drydepo [END]')
   if (chm_timings_L) call timing_stop_omp(377)
   !-----------------------------------------------------------------

   return
end subroutine mach_cam_drydep_main
