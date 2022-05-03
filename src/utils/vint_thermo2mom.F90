!-------------------------------------- LICENCE BEGIN ------------------------
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

subroutine vint_thermo2mom(fldout,fldin,vcoef,ni,nk)
   use phygridmap, only: phydim_nk
   implicit none
!!!#include <arch_specific.hf>
   !@Arguments
   integer, intent(in)  :: ni, nk         !# dimensions
   real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
   real,    intent(in)  :: fldin(ni,nk)   !# var on thermo levels
   real,    intent(out) :: fldout(ni,nk)  !# var on momentum levels
   !@Author L. Spacek (Dec 2007)
   !@Object
   ! Vertical interpolation of T/Q to momentum or thermo levels
   !*@/
   integer :: k, km1

   do k=nk,1,-1
      km1 = max(1,k-1)
      fldout(:,k) = fldin(:,k) - vcoef(:,k,1) * (fldin(:,k) - fldin(:,km1))
   enddo

   return
end subroutine vint_thermo2mom


!/@*
subroutine vint_mom2thermo(fldout,fldin,vcoef,ni,nk)
   use phygridmap, only: phydim_nk
   implicit none
!!!#include <arch_specific.hf>
   !@Arguments
   integer, intent(in)  :: ni, nk         !#  dimensions
   real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
   real,    intent(in)  :: fldin(ni,nk)   !# var on momentum levels
   real,    intent(out) :: fldout(ni,nk)  !# var on thermo levels
   !@Author L. Spacek (Dec 2007)
   !@Object
   ! Vertical interpolation of T/Q to momentum or thermo levels
   !*@/
   integer :: k, kp1

   do k=1,nk
      kp1 = min(nk,k+1)
      fldout(:,k) = fldin(:,k) + vcoef(:,k,2) * (fldin(:,kp1) - fldin(:,k))
   enddo

   return
end subroutine vint_mom2thermo
