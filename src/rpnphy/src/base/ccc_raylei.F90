
subroutine ccc_raylei(taur, ib, dp, il1, il2, ilg, lay)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, ib, il1, il2, ibm1, k, i
   real ri(3), taur(ilg,lay), dp(ilg,lay)

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Object Rayleigh scattering for bands2-bands4, near infrared region

   !@Arguments
   ! taur  rayleigh optical depth
   ! dp    air mass path for a model layer (explained in raddriv).

   !----------------------------------------------------------------------

   data ri / .16305e-04, .17997e-05, .13586e-06 /

   ibm1 = ib - 1
   do  k = 1, lay
      do  i = il1, il2
         taur(i,k) =  ri(ibm1) * dp(i,k)
      enddo
   enddo
   return
end subroutine ccc_raylei
