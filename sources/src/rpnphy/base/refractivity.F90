!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

subroutine refractivity2(zdct_bh, zdct_count, zdct_lvl, zdct_lvlmax, &
        zdct_lvlmin, zdct_sndmax, zdct_sndmin, zdct_str, zdct_thick, &
        zdct_trdmax, zdct_trdmin, zdct_moref, &
        zpmoins, zgz_plus, zsigm, ztplus, zhuplus, ni, nk)
   use tdpack_const
   use phy_options
   implicit none
!!!#include <arch_specific.hf>

   !@Arguments
   ! ni       horizontal running length
   ! nk       vertical dimension
   integer, intent(in) :: ni, nk
   real, dimension(ni), intent(in) :: zpmoins
   real, dimension(ni,nk), intent(in) :: zgz_plus, zsigm, ztplus, zhuplus
   real, dimension(ni), intent(inout) ::  zdct_bh, zdct_count, zdct_lvl, &
        zdct_lvlmax, zdct_lvlmin, zdct_sndmax, zdct_sndmin, zdct_str, &
        zdct_thick, zdct_trdmax, zdct_trdmin
   real, dimension(ni,nk), intent(inout) ::  zdct_moref

   !@Author Anna Glazer   (Nov 2008), S. Gaudreault (Mar 2012)
   !@Object IR/RF Refractivity
   !  M = 77.6/T{P+4810e/T} + z/R*10^6
   !  e and P in hPa, z and R in m, T in deg K
   !
   !  ref: Haack & Burk JAM 2001, vol 40 pp 673-687

   integer :: k, nk1, sk, sk1, sk2, i
   real :: a1, a2, a3, a4, ee(ni,nk), ppp(ni,nk), pp(ni,nk), geop(ni,nk)
   real, parameter :: safeguard = 1.e-5

   a1 = 77.6 * 0.01
   a2 = 4810.
   a3 = 1e+6 / (grav * rayt)
   a4 = 1 / grav
   geop = zgz_plus * grav

   !  pre-calculate water vapour pressure (ee)

   do k=1,nk
      do i=1,ni
         ppp(i,k) = zsigm(i,k)*zpmoins(i)
      end do
   end do

   call mfoefq3(ee, zhuplus, ppp, ni, nk, ni)

   !  Both ee (from mfoefq3) and pplus are in pa and should be in hpa
   !  cst a1 (a1 --> a1*0.01) has been modified to include this

   do k=1,nk
      do i=1,ni
         zdct_moref(i,k) = a1/ztplus(i,k)           &
              *(zsigm(i,k)*zpmoins(i) &
              + a2 * ee(i,k)/ztplus(i,k))   &
              + a3 * geop(i,k)
      end do
   end do

   !  compute above topography height (as used in physics)

   do k=1,nk
      do i=1,ni
         pp(i,k) = zdct_moref(i,k)
      end do
   end do

   !  Eliminate too weak ducts (strength lower than 1 M-unit)
   !  (Start from k=nk-1 (instead of k=nk))

   !  detecte first local max
   nk1=nk-1
   do i=1,ni
      do k=nk1,2,-1
         if (pp(i,k-1) < pp(i,k)) then
            zdct_lvlmax(i) = k
            go to 20
         endif
      end do
20 end do
   !
   !   the level of local max is dct_lvlmax, just above dct_moref has lower value,
   !   go up and look for the level when dct_moref starts to increase again
   !   the level of local min is dct_lvlmin
   !
   do i=1,ni
      sk = zdct_lvlmax(i)
      if (sk .gt. 0) then
         sk1 = sk-1
         do k=sk1,2,-1
            if (pp(i,k-1) > pp(i,k)) then
               zdct_lvlmin(i) = k
               go to 25
            endif
         end do
      endif
25 end do

   !  calcule duct strength

   do i=1,ni
      sk  = zdct_lvlmax(i)
      sk1 = zdct_lvlmin(i)
      if (sk .gt. 0) then
         zdct_str(i) = zdct_moref(i,sk) - zdct_moref(i,sk1)
      endif
   end do

   !  eliminate too weak ducts

   do i=1,ni
      sk = zdct_str(i)
      if (sk .lt. 1.) then
         zdct_str   (i) = 0.0
         zdct_lvlmax(i) = 0.0
         zdct_lvlmin(i) = 0.0
      endif
   end do

   !  for sfc based ducts :     calcule duct base height and thickness
   !  for non sfc based ducts : go down and find the first level (dct_lvl)
   !  such that the value of dct_moref at this level is smaler than at dct_lvlmin

   do i=1,ni
      sk  = zdct_lvlmax(i)
      sk1 = zdct_lvlmin(i)
      if (sk .gt. 0) then
         if (zdct_moref(i,sk1) .le. zdct_moref(i,nk)) then
            zdct_bh(i) = -100
            zdct_thick(i) = a4 * geop(i,sk1)
         else
            do k=sk,nk
               if (zdct_moref(i,sk1) .ge. zdct_moref(i,k)) then
                  zdct_lvl(i) = k
                  go to 30
               endif
            end do
         endif
      endif
30 end do

   !  for non sfc based ducts : calcule duct base height and thickness
   !  linear interpolation between levels : dct_lvl and dct_lvl-1
   !  dct_lvl > 0 means there is a non surface based duct

   do i=1,ni
      sk  = zdct_lvl(i)
      sk1 = sk-1
      sk2 = zdct_lvlmin(i)
      if (sk .gt. 0) then
         zdct_bh(i) = a4 * &
              (geop(i,sk) + (geop(i,sk1) - geop(i,sk)) / &
              (zdct_moref(i,sk1) -zdct_moref(i,sk) + safeguard) * &
              (zdct_moref(i,sk2)-zdct_moref(i,sk)))
         zdct_thick(i) = a4 * geop(i,sk2) - zdct_bh(i)
      endif
   end do

   !  find second local max

   do i=1,ni
      sk = zdct_lvlmin(i)
      if (sk .gt. 0) then
         do k=sk,2,-1
            if (pp(i,k-1) .lt. pp(i,k)) then
               zdct_sndmax(i) = k
               go to 40
            endif
         end do
      endif
40 end do

   !  find second local min

   do i=1,ni
      sk = zdct_sndmax(i)
      if (sk .gt. 0) then
         sk1 = sk-1
         do k=sk1,2,-1
            if (pp(i,k-1) .gt. pp(i,k)) then
               zdct_sndmin(i) = k
               go to 45
            endif
         end do
      endif
45 end do

   !  find third local max

   do i=1,ni
      sk=zdct_sndmin(i)
      if (sk .gt. 0) then
         do k=sk,2,-1
            if (pp(i,k-1) .lt. pp(i,k)) then
               zdct_trdmax(i)=k
               go to 60
            endif
         end do
      endif
60 end do

   !  find third local min

   do i=1,ni
      sk = zdct_trdmax(i)
      if (sk .gt. 0) then
         sk1 = sk-1
         do k=sk1,2,-1
            if (pp(i,k-1) .gt. pp(i,k)) then
               zdct_trdmin(i) = k
               go to 65
            endif
         end do
      endif
65 end do

   !  the number of ducts will be put in dct_count (DNBR)
   !  if dct_lvlmax, dct_sndmax, dct_trdmax have not zero values 
   !  these values indicate k's of local max that are followed higher by min

   do i=1,ni
      sk  = zdct_lvlmax(i)
      sk1 = zdct_sndmax(i)
      sk2 = zdct_trdmax(i)
      if (sk .gt. 0)  then 
         zdct_count(i) = zdct_count(i) +1
      endif
      if (sk1 .gt. 0) then 
         zdct_count(i) = zdct_count(i) +1
      endif
      if (sk2 .gt. 0) then
         zdct_count(i) = zdct_count(i) +1
      endif
   end do
   return
end subroutine refractivity2
