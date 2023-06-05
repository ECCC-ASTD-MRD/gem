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
! Fichier/File   : mach_cam_aeroact.ftn90
! Creation       : S. Gong, W. Gong, V. Bouchet, S. Gravel, and
!                  B. Pabla for GEM-MACH, June 2008
!
! Description    : Calculation of activated aerosol
!
! Extra info     : - First version created by S. Gong Jan 12 1998 for CAM
!                    A simplified version of ghan 1993 activation
!                    parameterization scheme with constant d parameter.
!                  - The cloud number density is now computed from jones et al
!                    1994 parameterization which relates the total aerosol
!                    number to cloud droplet number. (S. Gong, Dec 02, 1998)
!
! Arguments:  IN
!               pnk     -> no. z-direct. vertical levels
!               pni     -> number of longitude grid points
!               rhsize  -> Unactivated ambient aerosol wet
!               aeronum -> Number conconcenration
!               zmlwc   -> CWC content (bulk) (kg/kg)
!               roarow  -> air density (kg/m3)
!               ncp     -> cloud droplet number density from microphysics
!
!             OUT
!               Q_BIN   -> CWC per activated size bin (kg/m3)
!               RCRITS  -> Bin Number (+ Fraction Un-Activated)
!               CLSIZE  -> Activated ambient aerosol wet
!
!             LOCAL
!               CCN     -> cloud droplet number density [1/m3]
!               TOTNUM  -> Accumulated number of aerosol (1/m3)
!============================================================================
!
!!if_on
subroutine mach_cam_aeroact(q_bin, rhsize, rcrit, aeronum, zmlwc, &
                            roarow, clsize, ncp, ibulk, pni, pnk)
   use mach_cam_utils_mod, only: isize
!!if_off
   use chm_nml_mod,        only: chm_indirect_l
   
   implicit none
!!if_on
   integer(kind=4), intent (in) :: ibulk
   integer(kind=4), intent (in) :: pni, pnk
   real(kind=4),    intent(out) :: q_bin  (pni, pnk, isize)
   real(kind=4),    intent (in) :: rhsize (pni, pnk, isize)
   real(kind=4),    intent(out) :: rcrit  (pni, pnk)
   real(kind=4),    intent (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent (in) :: zmlwc  (pni, pnk)
   real(kind=4),    intent (in) :: roarow (pni, pnk)
   real(kind=4),    intent (in) :: ncp    (pni, pnk)
   real(kind=4),    intent(out) :: clsize (pni, pnk, isize)
!!if_off
!
!  local variables
!
   integer(kind=4) :: n, n1, l, i, isize_0
   real(kind=4)    :: totnum(pni, pnk, isize)
   real(kind=4)    :: h2o(pni, pnk), ccn(pni, pnk)
   real(kind=4)    :: rleft, rwi
   real(kind=8)    :: tona
   real, parameter :: cub = 1.0 / 3.0, rhopw_inv = 1.0e-3

!  Initializations
   clsize = 0.0
   totnum = 0.0
   q_bin  = 0.0
   rcrit  = real(isize + 1)

!  accumulate number of aerosols[/m^3] for each grid
   do n = isize, 1, -1
      do n1 = n, isize
         do l = 1, pnk
            do i = 1, pni
               totnum(i, l, n) = totnum(i, l, n) + aeronum(i, l, n1)  * roarow(i, l)
            end do
         end do
      end do
   end do

!  Find the activated bin number
!  define: rcrit = bin number + % un-actived
!  e.g.    rcrit = 4.3 means bin 4 and over are activated with 30% un-activated in bin 4.

   if (chm_indirect_l) then
!  cloud droplet number density from the microphysics scheme
      ccn = ncp
   else
!  cloud droplet number density [1/m3] from jones et al 1994
      do l = 1, pnk
         do i = 1, pni
            tona = dble(totnum(i, l, 1)) * 1.0d-6
            ccn(i, l) = real(375.0d6 * (1.0d0 - dexp(-2.5d-3 * tona)))
         end do
      end do
   end if

   do n = isize, 2, -1
      do l = 1, pnk
         do i = 1, pni
            if (ccn(i, l) > totnum(i, l, n) .and. ccn(i, l) <= totnum(i, l, n - 1)) then
               rleft = (totnum(i, l, n - 1) - ccn(i, l)) / (aeronum(i, l, n - 1) * roarow(i, l))
               if (zmlwc(i, l) > 0) then
                  rcrit(i, l) = real(n - 1) + min(1.0, rleft)
               end if
            end if
            if (ccn(i, l) > totnum(i, l, 1)) then
               rcrit(i, l) = 1.0
            end if
         end do
      end do
   end do

   where(ccn > 0.0)
     h2o = zmlwc * roarow / ccn
   elsewhere
     h2o = 0.0
   end where

!  Distribute bulk cloud liquid water [kg/kg] into activated bins according
!  to the cloud number density in each activated bin and estimate cloud
!  droplet radius.
   if (ibulk == 1) then
      do l = 1, pnk
         do i = 1, pni
            clsize(i, l, 1) = ((h2o(i, l) / 4.189) * rhopw_inv) ** cub
         end do
      end do
!
   else
!!!Cloud droplet radius at the interrupted bin where aerosol and cloud doplet 
!!!co-exist is saved in rcoex(i, l). For fully actviated bins, the rhsize is 
!!!updated.   
      do l = 1, pnk
         do i = 1, pni
            if (ccn(i, l) > 1.0e-7) then
               isize_0 = int(rcrit(i, l))
               do n = 1, isize
                  if (n >= isize_0) then
                     rwi = rhsize(i, l, n)
!                    if (n == isize_0) rcoex(i, l) = rhsize(i, l, n)

!  CLSIZE(I, L, N)=RHSIZE(I, L, N) we replace that line !!!!!!
!
!  explanation:

!  In mach_cam_AEROACT
!  CLSIZE in this bin is not calculated, instead it is given the aerosol dry
!  size. The wet size and dry size can be 100 times different. With CLSIZE
!  much smaller the droplet size for the cloud drops in this bin, the mass
!  transfer rate will be greatly enhanced (proportional to 1/r^3 - 1/r^2 at
!  the moment in UPAQR.f). What I did in my box model is to calculate the
!  CLSIZE in this bin the same way as for the other activated bins. Could you
!  implement this change and see how much it affects the threshold you set
!  for q_bin? Wanmin G .12 June 2000.

                     clsize(i, l, n) = ((rwi * rwi * rwi +  &
                                       h2o(i, l) / 4.189) * rhopw_inv) ** cub
                  end if
               end do
            end if
         end do
      end do
   end if

! Calculate q_bin:  CWC per activated size bin (kg/m3)
   do l = 1, pnk
      do i = 1, pni
         if (h2o(i, l) > 0.0) then
!            isize_0 = int(rcrit(i, l))
!
!   updated according to Wanmin's email Jul22, 2008 (see Mantis842)
!
            isize_0 = min(isize, int(rcrit(i, l)))
            if (isize_0 < isize) then
               do n = isize, isize_0, -1
                  q_bin(i, l, n) = h2o(i, l) * aeronum(i, l, n) * roarow(i, l)
               end do
               q_bin(i, l, isize_0) = h2o(i, l) * &
                                      (ccn(i, l) - totnum(i, l, isize_0 + 1))
            else
               q_bin(i, l, isize) = h2o(i, l) * ccn(i, l)
            end if
         end if
      end do
   end do
      
   return
end subroutine mach_cam_aeroact