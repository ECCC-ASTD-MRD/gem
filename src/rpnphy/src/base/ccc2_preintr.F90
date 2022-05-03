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
!-------------------------------------- LICENCE END ---------------------------

subroutine ccc2_preintr3(inpr, dir, qq, co2, rhc, il1, il2, ilg, lay)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, il1, il2
   real dir(ilg,lay), qq(ilg,lay), co2(ilg,lay), rhc(ilg,lay)
   integer inpr(ilg,lay)

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Revisions
   ! 001 L. Solheim (Apr 21,2008)  -  previous version printr2 for gcm15g:
   !                                - cosmetic change to add threadprivate
   !                                  for common block "trace", in support
   !                                  of "radforce" model option.
   ! 002  J. Li (Feb 09,2009)      -  new version for gcm15h:
   !                               - co2 passed directly, thus no need
   !                                 for "trace" common block.
   ! 003  P. vaillancourt, A-M Leduc  (Apr 2010) - Integrating the Li/Lazare 2009 version. (preintr3)

   !@Object
   !       This subroutine determines the interpolation points for the ratio
   !       of h2o and co2.

   !@Arguments
   !              - Output -
   ! inpr   number of the ratio level for the standard 5 ratios
   ! dir    interpolation factor for mass ratio of h2o / co2
   !        between two neighboring standard input ratios
   ! rhc    the ratio of the h2o mass mixing to co2 mass mixing
   !            - Input -
   ! qq           water vapor mass mixing ratio
   ! co2          co2 mass mixing ratio
   ! il1          1
   ! il2          horizontal dimension
   ! ilg          horizontal dimension
   ! lay          number of model levels
   !----------------------------------------------------------------------

#include "ccc_tracegases.cdk"

   integer jendr, k, i, j, l, lp1
   real standr(5)

   data standr / .06,  .24,  1., 4., 16. /

   jendr  =  5
   DO_K: do k = 1, lay
      do i = il1, il2
         inpr(i,k)   =  0
         rhc(i,k)    =  qq(i,k) / max(co2(i,k), 1.0E-10)

         do j = 1, jendr
            if (rhc(i,k) .gt. standr(j)) then
               inpr(i,k) =  inpr(i,k) + 1
            endif
         enddo

         l   =  inpr(i,k)
         lp1 =  l + 1
         if (l .ge. 1 .and. l .lt. 5) then
            dir(i,k)  = (rhc(i,k) - standr(l)) / &
                 (standr(lp1) - standr(l))
         else
            !---------------------------------------------------------------
            !     dir is not used with values of {0,5} in tlinehc, but we
            !     initialize here to avoid problems with nan when used
            !     in multitasking mode.
            !---------------------------------------------------------------
            dir(i,k)  =  0.0
         endif

      enddo
   enddo DO_K

   return
end subroutine ccc2_preintr3
