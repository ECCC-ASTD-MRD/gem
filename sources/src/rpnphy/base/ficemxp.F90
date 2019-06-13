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

subroutine ficemxp2(fice,tf,df,t,n,nk)

   implicit none
!!!#include <arch_specific.hf>

   ! Arguments
   integer, intent(in) :: n                             !horizontal dimension
   integer, intent(in) :: nk                            !vertical dimension
   real, dimension(n,nk), intent(in) :: t               !dry air temperature (K)
   real, dimension(n,nk), intent(out) :: fice           !fraction of ice
   real, dimension(n,nk), intent(out) :: tf             !threshold for saturation wrt ice or liquid
   real, dimension(n,nk), intent(out) :: df             !derivative of fraction wrt T for ice or liquid

   !Author
   !          J. Mailhot (Jan 2000)
   
   !Object
   !          Calculate the fraction of ice and set the values of threshold
   !          and derivative w/r to T for computation of saturation values
   !          in the presence of mixed phases.

   !Notes
   !          Based on the definition in:
   !          - Burk et al. 1997, JGR 102, 16529-16544
   !          and the observations of:
   !          - Curry et al. 1990, Int. J. Climatol. 10, 749-764.
   !          - Curry et al. 1997, JGR 102, 13851-13860.
   !
   !          For F (fraction of ice), linear variation between Tmin and Tmax
   !          For TF and DF, values are set such that saturation is w/r to liquid for T > Tmin
   !                 "         "             "        saturation is w/r to ice    for T < Tmin

   ! Local parameter declarations
   real, parameter :: TMIN=248.16,TMAX=258.16

   ! Local variable declarations
   integer :: j,k
   real :: dt
   
   ! Compute ice fraction and saturation measures
   dt = 1.0/(tmax-tmin)
   do k=1,nk
      do j=1,n
         fice(j,k) = max(0.0,min(1.0,(tmax-t(j,k))*dt))
         tf(j,k) = fice(j,k)
         df(j,k) = -dt
         if( t(j,k).lt.tmin .or. t(j,k).gt.tmax) df(j,k) = 0.0
      end do
   end do

   return
end subroutine ficemxp2
