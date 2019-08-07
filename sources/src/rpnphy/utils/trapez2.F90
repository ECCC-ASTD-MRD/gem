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

real function TRAPEZ2(DEL,F,N,NM)
   implicit none
!!!#include <arch_specific.hf>
   integer N, NM, I
   real F(N),DEL(NM)

   !@Author L.Garand (1989)

   !@Revision
   ! 001      G.Pellerin(Mar90)Standard documentation

   !@Object calculate integrals

   !@Arguments
   !          - Input -
   ! DEL      thickness of layers
   ! F        variable to be integrated
   ! N        dimension of F
   ! NM       dimension of DEL

   TRAPEZ2=0.
   do I=1,NM
      TRAPEZ2=TRAPEZ2+(F(I)+F(I+1))/2.*DEL(I)
   enddo
   return
end function TRAPEZ2
