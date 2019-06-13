!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------

subroutine rigrad3(shr2,tve,qe,zn,ps,se,z,rig,gama,gamaq,ribkg,n,nk)
   use tdpack
   implicit none
!!!#include <arch_specific.hf>

   ! Input arguments
   integer, intent(in) :: n,nk                      !slab dimensions
   real, dimension(n,nk), intent(in) :: shr2        !square of the wind shear (s-2)
   real, dimension(n,nk), intent(in) :: tve         !virtual temperature on energy levels (K)
   real, dimension(n,nk), intent(in) :: qe          !specific humidity on energy levels (kg/kg)
   real, dimension(n,nk), intent(in) :: zn          !mixing length (m)
   real, dimension(n), intent(in) :: ps             !surface pressure (Pa)
   real, dimension(n,nk), intent(in) :: se          !sigma value for energy levels
   real, dimension(n,nk), intent(in) :: z           !height of energy levels (m AGL)
   logical, intent(in) :: ribkg                     !compute a mixing length-based background Ri

   ! Output arguments
   real, dimension(n,nk), intent(out) :: rig        !gradient Richardson number
   real, dimension(n,nk), intent(out) :: gama,gamaq !temperature,moisture profile adjusmtent

   ! Internal variables
   integer :: i,k,down,up
   real :: te
   real, dimension(n) :: buoy
   real, dimension(n,nk) :: thetav,rig_bkg

   ! Compute gradient Richardson number profile
   do k=1,nk
      thetav(:,k) = tve(:,k) * (1e5/(ps(:)*se(:,k)))**CAPPA
   enddo
   do k=2,nk-1
      buoy = (thetav(:,k+1)-thetav(:,k-1))/(z(:,k+1)-z(:,k-1))
      rig(:,k) = (GRAV / thetav(:,k)) * buoy / shr2(:,k)
   enddo
   rig(:,1) = (GRAV / thetav(:,1)) * ((thetav(:,2)-thetav(:,1))/(z(:,2)-z(:,1))) / shr2(:,1)
   rig(:,nk) = (GRAV / thetav(:,nk)) * ((thetav(:,nk)-thetav(:,nk-1))/(z(:,nk)-z(:,nk-1))) / shr2(:,nk)

   ! Create mixing length-based Ri values
   BACKGROUND_PROFILE: if (ribkg) then
      do i=1,n
         do k=1,nk     
            up = k
            do while (up > 1 .and. z(i,up)-z(i,k) < zn(i,k))
               up = up - 1
            enddo
            up = min(up+1,k)
            down = k
            do while (down < nk .and. z(i,k)-z(i,down) < zn(i,k))
               down = down + 1
            enddo
            down = max(down-1,k)
            rig_bkg(i,k) = sum(rig(i,up:down))/real(down-up+1)
         enddo
      enddo
      rig = rig_bkg 
   endif BACKGROUND_PROFILE

   ! Compute (gamma) profile adjustment terms
   do k=1,nk
      do i=1,n
         te = fottv(tve(i,k),qe(i,k))
         gama(i,k) = -GRAV / (CPD * (tve(i,k)/te))
         gamaq(i,k) = 0.
      enddo
   enddo

end subroutine rigrad3

