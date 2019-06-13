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

subroutine blweight2(w,s,ps,n,nk)
   implicit none
!!!#include <arch_specific.hf>

   ! Arguments
   integer, intent(in) :: n                             !horizontal dimension
   integer, intent(in) :: nk                            !vertical dimension
   real, dimension(n,nk), intent(in) :: s               !sigma levels
   real, dimension(n), intent(in) :: ps                 !surface pressure (Pa)
   real, dimension(n,nk), intent(out) :: w              !weighting function [0,1]
   
   !Author
   !          J. Mailhot and B. Bilodeau (Dec 2001)
   
   !Revision
   !001       A-M. Leduc and S. Belair (Jun 2003) - ps as argument
   !                   blweight ---> blweight2. Change weighting
   !                   profile from sigma to pressure dependent.
   
   !Object
   !          Compute a weighting profile to be used with the moist
   !          turbulence scheme.
   
   !Notes
   !          The profile is set to:
   !            1 in the lower part of the atmosphere (S .ge. SMAX) if pres .ge. pmax
   !            0 in the upper part of the atmosphere (S .le. SMIN) if pres .le. pmin
   !            (with a linear interpolation in-between)
   
   ! Local parameters
   real, parameter :: PMIN=45000,PMAX=55000

   ! Local variables
   integer :: j,k
   real :: pres

   ! Compute profile of weighting functions
   do k=1,nk
      do j=1,n
         pres = s(j,k)*ps(j)
         w(j,k) = (1.-(pmax-pres)/(pmax-pmin))
         w(j,k) = min(max(w(j,k),0.),1.)
      end do
   end do

   return
end subroutine blweight2
