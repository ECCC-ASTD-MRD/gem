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
!**S/P  SATTENU - CALCULATION OF SOLAR ATTENUATION
!
      subroutine ccc1_sattenu (atten, ib, ig, rmu, o3, &
                          qq, dp, dip, dt, dt0, &
                          inpt, gh, il1, il2, ilg, s1)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, ib, ig, il1, il2, i, ng, im
      real tau
      real atten(ilg), rmu(ilg), o3(ilg), qq(ilg), dp(ilg), dip(ilg), &
           dt(ilg), dt0(ilg), s1(ilg)
      integer inpt(ilg)
      logical gh
      integer isl
!
!Authors
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
! 001    M.Lazarre,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers
!
!Object
!
!        Calculation of solar attenuation above the model top level. for
!        band1 only o3 and o2 are considered, the contribution of other
!        gases is small. for band 3 and 4, co2 is considered for gh
!
!Arguments
!
! atten  attenuation factor for downward flux from toa to the
!        model top level
! o3     o3 mass mixing ratio
! qq     water vapor mass mixing ratio
! dp     here dp is only the pressure difference, different from
!        that defined in raddriv. so there is a factor 1.02
! dip    interpolation factor for pressure
! dt     layer temperature - 250 k
! dt0    temperature in moon layer - 250 k
!
!Implicites
!
#include "ccc_tracegases.cdk"
#include "ccc1_bands.cdk"
#include "ccc1_bandsh.cdk"
!
!*
      if (ib .eq. 1)                                                then
        if (gh)                                                     then
          if (ig .eq. 3)                                            then
            do 100 i = il1, il2
              tau      =  1.02 * (cs1o3gh(ig) * o3(i) + &
                          cs1o2gh3 * rmo2) * dp(i)
              atten(i) =  exp(- tau / rmu(i))
  100       continue
          else
            do 110 i = il1, il2
              tau      =  1.02 * cs1o3gh(ig) * o3(i) * dp(i)
              atten(i) =  exp(- tau / rmu(i))
  110       continue
          endif
        else
!
          if (ig .eq. 1)                                            then
            do 120 i = il1, il2
              tau      =  1.02 * (cs1o3(ig) * o3(i) + cs1o21 * rmo2) * &
                          dp(i)
              atten(i) =  exp(- tau / rmu(i))
  120       continue
          else
            do 130 i = il1, il2
              tau      =  1.02 * cs1o3(ig) * o3(i) * dp(i)
              atten(i) =  exp(- tau / rmu(i))
  130       continue
          endif
        endif
!
      else if (ib .eq. 2)                                           then
        if (ig .eq. 1)                                              then
          do 200 i = il1, il2
            atten(i)   =  1.0
  200     continue
        else
          ng = 6
          im = ig - 1
          isl=1
          call ccc1_attenue (atten, cs2o2gh(1,1,im), o3, qq, dp, dip, dt,dt0, &
                        rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
        endif
!
      else if (ib .eq. 3)                                           then
        ng = 3
        isl=1
        call ccc1_attenue (atten, cs3co2gh(1,1,ig), o3, qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
!
      else if (ib .eq. 4)                                           then
        ng = 3
        if (ig .ne. 4 .and. ig .ne. 6 .and. ig .ne. 8)              then
          if (ig .le. 3)  im =  ig
          if (ig .eq. 5)  im =  ig - 1
          if (ig .eq. 7)  im =  ig - 2
          if (ig .eq. 9)  im =  ig - 3
          isl=1
          call ccc1_attenue (atten, cs4co2gh(1,1,im), o3, qq, dp, dip,dt,dt0, &
                        rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
        else
          do 400 i = il1, il2
            atten(i)   =  1.0
  400     continue
        endif
      endif
!
      return
      end
