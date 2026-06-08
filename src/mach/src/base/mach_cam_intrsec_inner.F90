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
!  Modification:  This version expects only particle mass as input,
!                 and is called from an outer routine that subdivides
!                 the particle bins to achieve higher accuracy within
!                 this stage.  P.A. Makar, W. Gong, S. Gong, Feb 2008
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_cam_intrsec.ftn90
! Creation       : S. Gong, S. Gravel and B. Pabla for GEM-MACH, June 2008
!
! Description    : This module computes intersectional transport of aerosols
!                  due to comdensation or cloud processes
!
! Extra info     : First version created by S. Gong Aug 11 1997 for CAM
!
! Arguments:  IN
!               rtcond   -> Mass transfer rate onto each particle size bin
!               totmas   -> Total mass of aerosol in each bin
!               nn       -> Loop index on aerosol types
!              rtcond    -> Mass transfer rate onto each particle size bin
!               aeronum  -> Number concentration of aerosols
!                  mae   -> 0
!
!             IN/OUT
!              XROW      -> Particle tracer concentration in each bin before/after intersection tranport
!
!             MODULE: CAM_UTILS
!               nbnd     -> Number of (sub-divided)size bins
!============================================================================
!
!!if_on
 subroutine mach_cam_intrsec_inner(xrow, nn, rtcond, aeronum, mae, pni, pnk)
   use mach_cam_utils_mod, only: nbnd, nsb
!!if_off
   use chm_utils_mod,      only: chm_timestep
   use mach_cam_utils_mod, only: rhop0, icom, sub_pvol
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nn
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, nsb)
   real(kind=4),    intent   (in) :: rtcond (pni, pnk, nbnd)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, nbnd)
!!if_off
!
!  local variables
!
   integer(kind=4)   :: k, l, n, i, nt, no, nk, k0
   real(kind=4)      :: voij, vok, vokm1, vokp1
   real(kind=4)      :: rth  (nbnd)
   real(kind=4)      :: rgrid(pni, pnk, nsb)
!

   rgrid = xrow
   
   do l = 1 + mae, pnk
      do i = 1, pni
         do n = 1, nbnd

            if (aeronum(i, l, n) <= 0.0) cycle

!  new dry volume of size bin n
            voij = sub_pvol(n) + &
                   (rtcond(i, l, n) * chm_timestep / (aeronum(i, l, n) * rhop0(nn)))

            rth = 0.0
            do k = n, nbnd
               vok = sub_pvol(k)
               if (k == nbnd .and. voij >= vok) then
                  rth(k) = 1.0
               end if
               if (k < nbnd) then
                  k0 = min(k+1, nbnd)
                  vokp1 = sub_pvol(k0)
                  if (voij >= vok .and. voij < vokp1) then
                     rth(k) = vok / voij * (vokp1 - voij) / (vokp1 - vok)
                  endif
               end if
               if (k > 1) then
                  k0 = max(1, k-1)
                  vokm1 = sub_pvol(k0)
                  if (voij > vokm1 .and. voij < vok) then
                     rth(k) = 1.0 - rth(k0)
                  endif
               end if
               do nt = 1, icom
                  no = n + nbnd * (nt - 1)
                  nk = k + nbnd * (nt - 1)
                  if (rth(k) > 0.0) then
!  zero bin n for re-distribution
                     if (n == k) then
                        xrow(i, l, no) = max(0.0, xrow(i, l, no) - rgrid(i, l, no))
                     end if
!  distributed into bin k
                     xrow(i, l, nk) = xrow(i, l, nk) + rgrid(i, l, no) * rth(k)
                     if (nt == nn) then
                        xrow(i, l, nk) = xrow(i, l, nk) + rth(k) * rtcond(i, l, n) * chm_timestep
                     end if
                  end if
               end do
            end do

         end do
      end do
   end do

   return
end subroutine mach_cam_intrsec_inner
