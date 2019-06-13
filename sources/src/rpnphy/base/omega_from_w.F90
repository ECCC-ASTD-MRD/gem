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
!-------------------------------------- LICENCE END --------------------------


subroutine omega_from_w(dotp,ww,tt,qq,ps,sigma,ni,nk)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>

   integer                  :: ni,nk
   real, dimension(ni)      :: ps
   real, dimension(ni,nk)   :: dotp,ww,tt,qq,sigma

   !Author
   !          L.Spacek, April 2012
   !
   !Revisions
   ! 001
   !
   !Object
   !          Interface to convection/condensation
   !
   !Arguments
   !
   !          - Input -
   ! ww       vertical velocity in z-coordinate [m/s]
   ! tt       temperature
   ! qq       specific humidity
   ! ps       surface pressure
   ! sigma    vertical coordinate
   ! ni       horizontal dimension
   ! nk       vertical dimension
   !
   !          - Output -
   ! dotp     vertical velocity in pressure coordinates [Pa/s]

   real, dimension(ni,nk) :: tv, press

   call mfotvt(tv,tt,qq,ni,nk,ni)
   press(:,:) = spread( ps(:), dim=2, ncopies=nk )
   press(:,:) = sigma(:,1:nk)*press(:,:)
   dotp=-grav*ww*press/(rgasd*tv)

   return
end subroutine omega_from_w
