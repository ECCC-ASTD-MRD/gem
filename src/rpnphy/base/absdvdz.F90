!-------------------------------------- LICENCE BEGIN -------------------------
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

subroutine ABSDVDZ3(DUDZ2,U,V,TVE,SE,DSGDZ,SIGMA,N,M,NK)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>
   integer N,M,NK
   real DUDZ2(N,NK),U(M,NK),V(M,NK),TVE(N,NK)
   real SE(n,NK),DSGDZ(n,NK),SIGMA(N,NK)
   !
   !Author
   !          J. Cote (RPN 1983)
   !
   !Revision
   ! 001      J. Cote RPN(Nov 1984)SEF version documentation
   ! 002      M. Lepine  -  RFE model code revision project (Feb 87)
   ! 003      R. Benoit  -  'E' Levels intervals option (Mar 89)
   !                      -  To pass in parameter the workspace
   !                         WORK to DVRTEF
   ! 004      N. Brunet  (May90)
   !                Standardization of thermodynamic functions
   ! 005      C. Girard  (Mar93) Cleanup
   ! 006      B. Bilodeau (May 94) - New physics interface
   ! 007      C. Girard  (Jan96) - added DUDZ2 in computation
   !
   !Object
   !          to calculate the vertical shear of the wind
   !          DUDZ2 = module (Dwind/DZ) **2
   !
   !Arguments
   !
   !          - Output -
   ! DUDZ2    vertical shear of the wind squared (N,NK)
   !
   !          - Input -
   ! U        east-west component of the wind
   ! V        north-south component of the wind
   ! TVE      virtual temperature on 'E' level
   ! SIGMA    sigma levels
   ! DSGDZ    distance between the sigma levels
   ! W        W
   ! N        horizontal dimension
   ! M        1st dimension of U,V,TVE
   ! NK       vertical dimension
   !
   !Notes
   !   U and V are wind images
   !   S and U can share the same space if U is not required
   !   W and V can share the same space if V is not required

#include "clefcon.cdk"

   real DUDSG, DVDSG
   integer j,k

   do k=1,NK-1
      do j=1,N

         DUDSG = (U(j,k+1)-U(j,k))/(SIGMA(j,k+1)-SIGMA(j,k))
         DVDSG = (V(j,k+1)-V(j,k))/(SIGMA(j,k+1)-SIGMA(j,k))

         DSGDZ(j,k) = (GRAV/RGASD)*SE(j,k)/TVE(j,k)

         DUDZ2(j,k) = (DUDSG**2 + DVDSG**2) * DSGDZ(j,k) ** 2

      end do
   end do

   do j=1,N
      DUDSG = U(j,NK)/(SIGMA(j,NK)-1.)
      DVDSG = V(j,NK)/(SIGMA(j,NK)-1.)
      DSGDZ(j,NK) = (GRAV/RGASD)*SE(j,NK)/TVE(j,NK)
      DUDZ2(j,NK) = (DUDSG**2 + DVDSG**2) * DSGDZ(j,NK) ** 2
   end do

   return
end subroutine ABSDVDZ3
