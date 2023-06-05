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

subroutine ccc2_strandn3(tran, attn, attntop, rmu, dp, dt, o3, o2, &
     rmu3, ib, ig, lev1, il1, il2, ilg, lay, lev)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, lev, ib, ig, lev1, il1, il2
   real tran(ilg,2,lev), attn(ilg), attntop(ilg), rmu(ilg), &
        dp(ilg,lay), dt(ilg,lay), o3(ilg,lay), o2(ilg,lay), &
        rmu3(ilg)

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Revisions
   ! 001 L. Solheim , J. Li  (Apr 21,2008) -  previous version strandn2 for gcm15g:
   !                                       - cosmetic change to add threadprivate
   !                                         for common block "trace", in support
   !                                         of "radforce" model option.
   !                                       - more accuate treatment of o3
   ! 002  J. Li (Feb 09,2009) - new version for gcm15h:
   !                          - o2 passed directly, thus no need
   !                           for "trace" common block.
   ! 003  P. vaillancourt, A-M Leduc  (Apr 2010) - Integrating the  Li/Lazare 2009 version. (strandn3)
   !                                  add dto3 calculation

   !@Object
   !        Calculation of the downward flux from top level to 1 mb, no
   !        scattering effect is considered

   !@Arguments
   ! tran     transmisivity
   ! attn     attenuation factor for reducing solar flux from model
   !          top level to 1 mb
   ! attntop  attenuation factor for reducing solar flux from toa to
   !          model top level
   ! rmu      cosine of solar zenith angle
   ! dp       air mass path for a model layer (explained in raddriv).
   ! o3       o3 mass mixing ratio
   ! o2       o2 mass mixing ratio

#include "ccc_tracegases.cdk"
#include "ccc2_bands.cdk"


   integer i, k, kp1, lev1m1
   real tau, xx
   real dto3


   lev1m1 =  lev1 - 1

   if (ib .eq. 1) then

      !----------------------------------------------------------------------
      !     band (14500 - 50000 cm^-1), nongray gaseous absorption of o3
      !     ans o2.
      !----------------------------------------------------------------------

      do i = il1, il2
         tran(i,1,1)       =  attntop(i)
         tran(i,2,1)       =  attntop(i)
      enddo

      if (ig .eq. 1) then
         do k = 1, lev1m1
            kp1 = k + 1
            do i = il1, il2
               dto3          = dt(i,k) + 23.13
               xx            = (cs1o21 - 0.881e-05 * rmu3(i)) * o2(i,k)
               tau           = ((cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                    dto3 * cs1o3(3,ig))) * o3(i,k) + xx) * &
                    dp(i,k)
               tran(i,1,kp1) =  tran(i,1,k) * exp(- tau / rmu(i))
               tran(i,2,kp1) =  tran(i,1,kp1)
            enddo
         enddo

      else
         do k = 1, lev1m1
            kp1 = k + 1
            do i = il1, il2
               dto3          =  dt(i,k) + 23.13
               tau           =  (cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                    dto3 * cs1o3(3,ig))) * o3(i,k) * dp(i,k)
               tran(i,1,kp1) =  tran(i,1,k) * exp(- tau / rmu(i))
               tran(i,2,kp1) =  tran(i,1,kp1)
            enddo
         enddo
      endif

      !----------------------------------------------------------------------
      !     flux adjustment for region below 1 mb
      !----------------------------------------------------------------------

      do i = il1, il2
         attn(i)           =  tran(i,1,lev1)
      enddo

   else

      do k = 1, lev1
         do i = il1, il2
            tran(i,1,k)       =  1.0
            tran(i,2,k)       =  1.0
         enddo
      enddo

      do i = il1, il2
         attn(i)           =  1.0
      enddo
   endif

   return
end subroutine ccc2_strandn3
