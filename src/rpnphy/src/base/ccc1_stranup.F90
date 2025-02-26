!**S/P STRANUP - CALCULATION OF THE UPWARD SOLAR
!
      subroutine ccc1_stranup (refl, dp, o3, ib, ig, lev1, &
                          il1, il2, ilg, lay, lev)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lev, ib, ig, lev1, il1, il2, lev1m1, k, kp1, i
      real tau, dtr
      real refl(ilg,2,lev), dp(ilg,lay), o3(ilg,lay)
!
!Authors
!
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
!
! 001
!
!Object
!
!        Calculation of the upward solar flux above 1 mb, no scattering
!        effect is considered
!
!Arguments
!
! refl   reflectivity
! tau    optical depth
! dtr    direct transmission function
!        1.6487213 diffusivity factor (Li, 2000 JAS p753)
!
!Implicites
!
#include "ccc_tracegases.cdk"
#include "ccc1_bands.cdk"
!
!*
      lev1m1 =  lev1 - 1
!
      if (ib .eq. 1)                                                then
!
        if (ig .eq. 1)                                              then
          do 105 k = lev1m1, 1, - 1
            kp1 = k + 1
            do 100 i = il1, il2
              tau         = (cs1o3(ig) * o3(i,k) + cs1o21 * rmo2) * &
                             dp(i,k)
              dtr         =  exp( - 1.6487213 * tau)
              refl(i,1,k) =  refl(i,1,kp1) * dtr
              refl(i,2,k) =  refl(i,2,kp1) * dtr
  100       continue
  105     continue
!
        else
          do 115 k = lev1m1, 1, - 1
            kp1 = k + 1
            do 110 i = il1, il2
              tau         =  cs1o3(ig) * o3(i,k) * dp(i,k)
              dtr         =  exp( - 1.6487213 * tau)
              refl(i,1,k) =  refl(i,1,kp1) * dtr
              refl(i,2,k) =  refl(i,2,kp1) * dtr
  110       continue
  115     continue
!
        endif
!
      else
!
        do 250 k = lev1m1, 1, - 1
          kp1 = k + 1
          do 200 i = il1, il2
            refl(i,1,k)   =  refl(i,1,kp1)
            refl(i,2,k)   =  refl(i,2,kp1)
  200     continue
  250   continue
!
      endif
!
      return
      end
