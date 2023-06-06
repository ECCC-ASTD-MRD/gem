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
! Fichier/File   : mach_cam_condsoa.ftn90
! Creation       : S. Gong, S. Gravel and B. Pabla for GEM-MACH, June 2008
!
! Description    : SOA Condensation routine
!
! Extra info     : First version created by S. Gong Aug 29 1999 for CAM
!                  Use the same fraction of suphate condensation into
!                  each bin for the total soa.
!
! Arguments:  IN
!               throw   -> Temp
!               aeronum -> Number concentration of aerosols
!               ntr     -> Total number of trace substances (gases and aerosols)
!               roarow  -> Air density (kg/m3)
!               pcond   -> fractional condensation to each bin
!               soa     -> Secondary organics aerosols
!
!             OUT
!              RTCOND   -> Mass transfer rate onto each particle size bin
!
!             INOUT
!              XROW     -> Tracers concentration
!
!============================================================================
!
!!if_on
subroutine mach_cam_condsoa(aeronum, xrow, roarow, rtcond, pcond, soa, pres, &
                            tp, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
!!if_off
   use mach_cam_utils_mod,   only: iae_OC
   use mach_cam_headers_mod, only: mach_cam_intrsec_outer
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
   real(kind=4),    intent  (out) :: rtcond (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pcond  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: soa    (pni, pnk)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: tp     (pni, pnk)
!!if_off
!
!  local variables
!
   integer(kind=4)  :: n, l, i

!  condensation rate of soa for each bin
   do n = 1, isize
      do l = 1, pnk
         do i = 1, pni
            if (soa(i, l) > 0.0 .and. (aeronum(i, l, n) * roarow(i, l)) >  1.0e6) then
!  rate of soa condensed to each bin
               rtcond(i, l, n) = pcond(i, l, n) * soa(i, l)
            else
               rtcond(i, l, n) = 0.0
            end if
         end do
      end do
   end do

!  call to compute the intersectional transport due to condensation process.
   call mach_cam_intrsec_outer(xrow, iae_OC, rtcond, aeronum, 0, pres, &
                               tp, roarow, pni, pnk)
   return
end subroutine mach_cam_condsoa
