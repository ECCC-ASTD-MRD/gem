
subroutine ccc1_preintr (inpr, dir, qq, rhc, il1, il2, ilg, lay)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, il1, il2, jendr, k, i, j, l, lp1
   real rrmco2
   real dir(ilg,lay), qq(ilg,lay), rhc(ilg,lay), standr(5)
   integer inpr(ilg,lay)

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Object
   !       This subroutine determines the interpretion points for the ratio
   !       of h2o and co2.

   !@Arguments
   ! inpr   number of the ratio level for the standard 5 ratios
   ! dir    interpolation factor for mass ratio of h2o / co2
   !        between two neighboring standard input ratios
   ! q:     water vapor mass mixing ratio
   ! rhc    the ratio of the h2o mass mixing to co2 mass mixing
   !----------------------------------------------------------------------

#include "ccc_tracegases.cdk"

   data standr / .06,  .24,  1., 4., 16. /

   rrmco2 =  1.0 / (rmco2 + 1.e-10)
   jendr  =  5
   DO_K: do k = 1, lay
      do i = il1, il2
         inpr(i,k)   =  0
         rhc(i,k)    =  qq(i,k) * rrmco2

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

            !--------------------------------------------------------------
            !     dir is not used with values of {0,5} in tlinehc, but we
            !     initialize here to avoid problems with nan when used
            !     in multitasking mode.
            !--------------------------------------------------------------

            dir(i,k)  =  0.0
         endif

      enddo
   enddo DO_K

   return
end subroutine ccc1_preintr
