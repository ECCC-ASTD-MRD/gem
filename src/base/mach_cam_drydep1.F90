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
! Fichier/File   : mach_cam_drydep1.ftn90
! Creation       : S. Gong, B. Pabla and S. Gravel for GEM-MACH, June 2008
! Description    : Mass balance calculation for each layer.
!
! Extra info     : - First version created by S. Gong Jul 08 1994 for CAM
!                    Method:
!                     1. For each level, compute the falling distance in the
!                        time step far(I)
!                     2. Check the layer thichness below the level and determine
!                        which level it will land
!                     3. Add the mass to the level
!
!                  - Vectorized the whole program and add working spaces.
!                    (S. Gong, Jan 19, 1996)
!
! Arguments:  IN
!
!                pdepv   -> gravitational settling velocity
!                thlev   -> Layer thickness [m]
!                rho     -> Layer air density [kg/m3]
!
!             IN/OUT
!
!                XROW    -> species concentrations
!
!             OUT
!
!                DRYFLX  -> Dry deposition flux
!
!============================================================================
!
!!if_on
subroutine mach_cam_drydep1(xrow, dryflx, pdepv, thlev, rho, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr, icom
!!if_off
   use chm_utils_mod,        only: chm_timestep
   use mach_cam_utils_mod,   only: iae1, iae2, mwt_aero, tmin
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent  (out) :: dryflx(pni, icom)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: rho   (pni, pnk)
!!if_off
!
! Local variables
!
   integer(kind=4) :: nt, n, np, l, i, ll1, lar
   real(kind=4)    :: dhl, far, mass_per_area
   real(kind=4)    :: pnew(pni, pnk, ntr)
   real(kind=4)    :: rtdry(pni, pnk, ntr)

   pnew  = 0.0
   dryflx = 0.0

   do nt = 1, icom
      do n = 1, isize
         np = isize * (nt - 1) + n + (iae1 - 1)
         do i = 1, pni
            do l = 1, pnk - 1
               far = 0.0
               lar = 1
               dhl = pdepv(i, l, n) * chm_timestep
               mass_per_area = thlev(i, l) * rho(i, l) 
!  create deposition tendencies for each layer
               rtdry(i, l, np) = xrow(i, l, np) * &
                                 (exp(-dhl / thlev(i, l)) - 1.0) / chm_timestep
               do ll1 = l + 1, pnk
                  far = far + thlev(i, ll1)
                  if (dhl <= far .and. lar == 1) then
!  not reaching the ground, lar = 0
!  assign the mass to the layer it will be in during the time step
                     pnew(i, ll1, np) = -rtdry(i, l, np) * mass_per_area / &
                                         (thlev(i, ll1) * rho(i, ll1))
                     lar = 0
                  end if
               end do
            end do

            rtdry(i, pnk, np) = xrow(i, pnk, np) * &
                                  (exp(-pdepv(i, pnk, n) * chm_timestep / &
                                       thlev(i, pnk)) - 1.0) / chm_timestep
!  evaluate the dry deposition flux
            dryflx(i, nt) = dryflx(i, nt) + &
                            rtdry(i, pnk, np) * chm_timestep * thlev(i, pnk)
         end do
      end do
   end do

!  update and return the aerosol concentration at each layer
   do np = iae1, iae2
      do l = 1, pnk
         do i = 1, pni
            xrow(i, l, np) = max(xrow(i, l, np) + (rtdry(i, l, np) + &
                                 pnew(i, l, np)) * chm_timestep, tmin)
         end do
      end do
   end do

   ! convert dryflx to mol/m2
   do nt = 1, icom
      do i = 1, pni
         dryflx(i, nt) = dryflx(i, nt) * 1.e3 * rho(i, pnk) /  mwt_aero(nt)
      end do
   end do

   return
end subroutine mach_cam_drydep1
