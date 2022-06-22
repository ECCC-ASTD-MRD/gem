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
      subroutine ccc2_sattenu4 (atten, ib, ig, rmu, o3, co2, &
                          ch4, o2, dp, dip, dt, dt0, &
                          inpt, gh, il1, il2, ilg)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, ib, ig, il1, il2
      real tau,dto3
      real atten(ilg), rmu(ilg), o3(ilg), co2(ilg), ch4(ilg),o2(ilg), &
           dp(ilg), dip(ilg),dt(ilg), dt0(ilg)
      integer inpt(ilg)
      logical gh
!
!Authors
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
! 001    M.Lazare  (may 05,2006)- previous version sattenu2 for gcm15e/f:
!                         - pass integer variables "init" and
!                           "nit" instead of actual integer
!                           values, to "attenue" routines.
! 002    M.Lazare,J.Li,L.Solheim (apr 21,2008) - previous version sattenu3 for gcm15g:
!                         - cosmetic change to add threadprivate
!                            for common block "trace", in support
!                            of "radforce" model option.
!                          - calls attenue3 instead of attenue2.
!                          - update o3 and add ch4 effect.
! 003    M.Lazarre,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers
!
! 004    J.Li (feb 09,2009) - new version for gcm15h:
!                           - 3d ghg implemented, thus no need
!                             for "trace" common block or
!                             temporary work arrays to hold
!                             mixing ratios of ghg depending on
!                             a passed, specified option.
!                           - calls attenue4 instead of attenue3.
! 005    P.vaillancourt, A-M.Leduc (apr 2010)- Integrating the Li/Lazare 2009 version. (sattenu4).
!                           - add Co2,ch4,o2, remove qq and s1.
!                           - calculation of dto3.

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
! o3     o3 mass mixing ratio above the modele top level
! co2    co2 mass mixing ratio at model top layer
! o2     o2 mass mixing ratio at model top layer
! ch4    ch4 mass mixing ratio
! dp     here dp is only the pressure difference, different from
!        that defined in raddriv. so there is a factor 1.02
! dip    interpolation factor for pressure
! dt     layer temperature - 250 k
! dt0    temperature in moon layer - 250 k
!
!Implicites
!
#include "ccc_tracegases.cdk"
#include "ccc2_bands.cdk"
#include "ccc2_bandsh.cdk"
!
       integer   i, im, isl
!*
      if (ib .eq. 1)                                                then
        if (gh)                                                     then
          if (ig .eq. 3)                                            then
            do 100 i = il1, il2
              dto3     =  dt(i) + 23.13
              tau      =  1.02 * ((cs1o3gh(1,ig) + dto3 * &
                         (cs1o3gh(2,ig) + dto3 * cs1o3gh(3,ig))) * &
                          o3(i) +  cs1o2gh3 * o2(i)) * dp(i)
              atten(i) = exp(- tau / rmu(i))
  100       continue
          else
            do 110 i = il1, il2
              dto3     =  dt(i) + 23.13
              tau      =  1.02 * (cs1o3gh(1,ig) + dto3 * &
                         (cs1o3gh(2,ig) + dto3 * cs1o3gh(3,ig))) * &
                          o3(i) * dp(i)
              atten(i) = exp(- tau / rmu(i))
  110       continue
          endif
        else
!
          if (ig .eq. 1)                                            then
            do 120 i = il1, il2
              dto3     =  dt(i) + 23.13
              tau      =  1.02 * ((cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                          dto3 * cs1o3(3,ig))) * o3(i) + &
                          cs1o21 * o2(i)) * dp(i)
              atten(i) =  exp(- tau / rmu(i))
  120       continue
          else
            do 130 i = il1, il2
              dto3     =  dt(i) + 23.13
              tau      =  1.02 * (cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                          dto3 * cs1o3(3,ig))) * o3(i) * dp(i)
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
          im = ig - 1
          isl = 1
          call ccc2_attenue4 (atten, cs2o2gh(1,1,im), o2, dp, dip, dt, &
                         dt0, rmu, inpt, ntl, isl, il1, il2, ilg)
        endif
!
      else if (ib .eq. 3)                                           then
        isl = 1
        call ccc2_attenue4 (atten, cs3co2gh(1,1,ig), co2, dp, dip, dt, dt0, &
                       rmu, inpt, ntl, isl, il1, il2, ilg)
!
      else if (ib .eq. 4)                                           then
        isl = 1
        if (ig .ne. 4 .and. ig .ne. 6 .and. ig .ne. 8)              then
          if (ig .le. 3)  im =  ig
          if (ig .eq. 5)  im =  ig - 1
          if (ig .eq. 7)  im =  ig - 2
          if (ig .eq. 9)  im =  ig - 3
          call ccc2_attenue4 (atten, cs4co2gh(1,1,im), co2, dp, dip, dt, &
                         dt0, rmu, inpt, ntl, isl, il1, il2, ilg)
        else if (ig .eq. 4)                                         then
          call ccc2_attenue4 (atten, cs4ch4gh, ch4, dp, dip, dt, dt0, rmu, &
                         inpt, ntl, isl, il1, il2, ilg)
        else
          do 400 i = il1, il2
            atten(i)   =  1.0
  400     continue
        endif
      endif
!
      return
      end
