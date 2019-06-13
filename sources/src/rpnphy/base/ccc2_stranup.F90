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
!**S/P STRANUP - CALCULATION OF THE UPWARD SOLAR
!
      subroutine ccc2_stranup3 (refl, dp, dt, o3, o2, ib, ig, lev1, &
                          il1, il2, ilg, lay, lev)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lev, ib, ig, lev1, il1, il2
      real refl(ilg,2,lev), dp(ilg,lay), dt(ilg,lay), o3(ilg,lay), &
           o2(ilg,lay)
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
! 001  L. Solheim, J. Li (Apr 21,2008) - previous version stranup2 for gcm15g:
!                                - cosmetic change to add threadprivate
!                                  for common block "trace", in support
!                                  of "radforce" model option.
!                                - more accuate treatment of o3 .
! 002  J. Li (Feb 09,2009) -     new version for gcm15h:
!                               - o2 passed directly, thus no need
!                                 for "trace" common block.
! 003  P. vaillancourt, A-M Leduc  (Apr 2010) - Integrating the  Li/Lazare 2009 version. (stranup3)
!                                   - add calculation of dto3
!
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
#include "ccc2_bands.cdk"
!
!*
      integer lev1m1, k, kp1, i
      real tau, dtr, dto3

      lev1m1 =  lev1 - 1
!
      if (ib .eq. 1)                                                then
!
        if (ig .eq. 1)                                              then
          do 105 k = lev1m1, 1, - 1
            kp1 = k + 1
            do 100 i = il1, il2
              dto3        =  dt(i,k) + 23.13
              tau         = ((cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                             dto3 * cs1o3(3,ig))) * o3(i,k) + &
                             cs1o21 * o2(i,k)) * dp(i,k)
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
              dto3        =  dt(i,k) + 23.13
              tau         = (cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                             dto3 * cs1o3(3,ig))) * o3(i,k) * dp(i,k)
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
