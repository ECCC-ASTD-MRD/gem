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

subroutine riflux2(R,PRI,N,NK)
   implicit none
!!!#include <arch_specific.hf>
   integer N,NK
   real R(N,NK),PRI(N,NK)

   !@Author J. Cote (RPN 1983)
   !@Revision
   ! 001      J. Cote RPN(Nov 1984)SEF version documentation
   ! 002      J. Mailhot RPN(Feb 1985)'Clipping' of GAMA (CG)
   ! 003      J. Mailhot RPN(Mar 1985)Derive by logarithm
   ! 004      M. Lepine - RFE model code revision project (Feb 87)
   ! 005      R. Benoit - Option of 'E' levels interval (Mar 89)
   ! 006      Y. Delage  (May89) Revision of the vertical diffusion
   ! 007      J. Mailhot RPN(Feb 1990)Compatible with *RIGRAD*
   ! 008      C. Girard (March 1993) Clean-up
   !@Object calculate the flux Richardson number
   !@Arguments
   !          - Input/Output -
   ! R        gradient Richardson number as input (if STAGE, R(*,NK)=0)
   !          flux Richardson number as output (if STAGE, R(*,NK)=0)
   !          - Input -
   ! PRI      inverse of generalized Prandtl number (KT/KM)
   ! N        horizontal dimension
   ! NK       vertical dimension

#include "clefcon.cdk"

   integer J,K

   do K=1,NK
      do J=1,N
         R(J,K)=PRI(J,K)*R(J,K)
      enddo
   enddo

   return
end subroutine riflux2
