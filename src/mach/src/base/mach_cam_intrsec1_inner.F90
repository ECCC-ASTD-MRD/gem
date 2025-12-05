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
!  Modification:  subdivision of the bin distribution is now done inside
!  mach_cam_intrsec1_outer.ftn90.   Consequently, the size of xrow
!  has decreased here to include only the particle bins.
!  P.A. Makar, W. Gong, S. Gong, February 2009.
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_cam_intrsec1.ftn90
! Creation       : S. Gong, W. Gong, V. Bouchet, S. Gravel and
!                  B. Pabla for GEM-MACH, June 2008
! Description    : This module computes intersectional transport of aerosols
!                  due to comdensation or cloud processes
!
! Extra info     : - First version created by S. Gong Aug 11 1997 for CAM
!                  - Modified for mass redistribution after cloud chemistry
!                    process (W. Gong, Oct 2000)
!                    1. change from volume fraction to mass fraction;
!                    2. extend to moving mass in either directions;
!                    3. deal with net mass change in a single aerosol
!                       as opposed to one component at a time;
!                    4. correcting an error to insure mass conservation
!                       in the case of skipping bins.
!                  - Correction to the reallocation loop (V. Bouchet, Jul 2000)
!
! Arguments:  IN
!               daqchm  -> Combined mass change in each particle size bin
!               aeronum -> Number concentration of aerosols
!               rmass   -> Single aerosol mass in each size bin before cloud chemistry
!
!             OUT
!
!             IN/OUT
!               XROW    -> Particle tracer concentration in each bin before/after intersection tranport
!
!============================================================================
!
!!if_on
subroutine mach_cam_intrsec1_inner(xrow, daqchm, aeronum, rmass, pni, pnk)
   use mach_cam_utils_mod, only: nbnd, nsb
!!if_off
   use mach_cam_utils_mod, only: icom
   implicit none
!!if_on
   integer(kind=4),  intent   (in) :: pni, pnk
   real(kind=4),     intent(inout) :: xrow   (pni, pnk, nsb)
   real(kind=4),     intent   (in) :: daqchm (pni, pnk, nbnd)
   real(kind=4),     intent   (in) :: aeronum(pni, pnk, nbnd)
   real(kind=4),     intent   (in) :: rmass  (pni, pnk, nbnd)
!!if_off
!
!  local variables
!
   integer(kind=4) :: k0
   integer(kind=4) :: k, pos, l, il, n, i, nt, no, nk
   real(kind=4)    :: rcond, rmoij, rmok, rmokp1, rmokm1
   real(kind=4)    :: rgrid(pni, pnk, nsb)
   real(kind=4)    :: rth(pni, pnk, nbnd)
   real(kind=4), parameter    :: smf = 0.001, smf1 = 1.0e-35

   rmokm1   = 0.0
   rmokp1   = 0.0
   rth      = 0.0

   do pos = 1, nsb
      do l = 1, pnk
         do il = 1, pni
            rgrid(il, l, pos) = xrow(il, l, pos)
         end do
      end do
   end do

!  condensation rate
   do n = 1, nbnd
      rth = 0.0

!  growth
      do k = n, nbnd
         do l = 1, pnk
            do i = 1, pni
               rcond = daqchm(i, l, n)

!  new dry mass (single aerosol) of size bin n
               rmoij = max(0.0, rmass(i, l, n) + rcond / (aeronum(i, l, n) + smf))
               rmok = rmass(i, l, k)
               if (k == nbnd .and. rmoij >= rmok) then
                  rth(i, l, k) = 1.0
               end if
               if (k < nbnd) then
                  k0 = min ( k+1, nbnd )
                  rmokp1 = rmass(i, l, k0)
                  if (rmoij >= rmok .and. rmoij < rmokp1) then
                     rth(i, l, k) = rmok / (rmoij + smf1) * (rmokp1 - rmoij) / (rmokp1 - rmok)
                  end if
               end if
               if (k > 1) then
                  k0 = max ( 1, k-1 )
                  rmokm1 = rmass(i, l, k0)
                  if( rmoij > rmokm1 .and. rmoij < rmok) then
                     rth(i, l, k) = 1.0 - rth(i, l, k0)
                  end if
               end if
            end do
         end do
      end do

!  reduction
      do k = n, 1, -1
         do l = 1, pnk
            do i = 1, pni
               rcond = daqchm(i, l, n)
!  new dry mass (single aerosol) of size bin n
               rmoij = max(0.0, rmass(i, l, n) + rcond / (aeronum(i, l, n) + smf))
               rmok = rmass(i, l, k)
               if (k == 1 .and. rmoij <= rmok) then
                  rth(i, l, k) = 1.0
               end if
               if (k > 1) then
                  rmokm1 = rmass(i, l, k - 1)
                  if (rmoij >= rmokm1 .and. rmoij < rmok) then
                     rth(i, l, k) = rmok / (rmoij + smf1) * (rmoij - rmokm1) / (rmok - rmokm1)
                  end if
               end if
               if (k < nbnd) then
                  rmokp1 = rmass(i, l, k + 1)
                  if (rmoij > rmok .and. rmoij < rmokp1) then
                     rth(i, l, k) = 1.0 - rth(i, l, k + 1)
                  end if
               end if
            end do
         end do
      end do

      do k = 1, nbnd
         do nt = 1, icom
            no = n + nbnd * (nt - 1) 
            nk = k + nbnd * (nt - 1) 
            do l = 1, pnk
               do i = 1, pni
!  zero bin n for re-distribution
                  if (n == k) then
                     xrow(i, l, no) = max(0.0, xrow(i, l, no) - rgrid(i, l, no))
                  end if
!  distributed into bin k
                  xrow(i, l, nk) = max(0.0, xrow(i, l, nk) + rgrid(i, l, no) * rth(i, l, k))
               end do
            end do
         end do
      end do
   end do

   return
end subroutine mach_cam_intrsec1_inner
