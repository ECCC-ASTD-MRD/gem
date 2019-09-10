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
!**S/P LATTENU - ATTENUATION FOR THE DOWNWARD FLUX
!
      subroutine ccc1_lattenu (atten, ib, ig, o3, qq, &
                          dp, dip, dt, dt0, inpt, &
                          il1, il2, ilg, s1, rmu)

      implicit none
!!!#include <arch_specific.hf>
!
      integer il1, il2, ilg, ib, ig, ng, i
      real atten(ilg), o3(ilg), qq(ilg), dp(ilg), dip(ilg), dt(ilg), &
           dt0(ilg), s1(ilg), rmu(ilg)
      integer inpt(ilg),isl
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
!        Calculation of the attenuation for the downward flux above the
!        model top level. since the temperature at 0.005 mb is unknown we
!        assume it is the same as that of model top level
!
!Arguments
!
! atten  for solar: the attenuation factor for downward flux from
!        toa to the model top level; for longwave: the optical
!        / diffuse factor
! dp     here dp is only the pressure difference, different from
!        that defined in raddriv. so there is a factor 1.02
! o3     o3 mass mixing ratio
! qq     water vapor mass mixing ratio
! dip    interpolation factor for pressure between two
!        neighboring standard input data pressure levels
! dt     layer temperature - 250 k
! dt0    temperature in moon layer - 250 k
!----------------------------------------------------------------------
!
!*
#include "ccc1_bandlh.cdk"
!
      if (ib .eq. 1)                                                then
        ng = 3
        isl=2
        call ccc1_attenue (atten, cl1co2gh(1,1,ig), o3, qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
!
      else if (ib .eq. 2)                                           then
        do 200 i = il1, il2
          atten(i) =  0.0
  200   continue
!
      else if (ib .eq. 3)                                           then
        ng = 1
        isl=2
        call ccc1_attenue (atten, cl3h2ogh(1,1,ig), o3, qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
!
      else if (ib .eq. 4)                                           then
        do 400 i = il1, il2
          atten(i) =  0.0
  400   continue
!
      else if (ib .eq. 5)                                           then
        ng = 2
        isl=2
        call ccc1_attenue (atten, cl5o3gh(1,1,ig), o3, qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
!
      else if (ib .eq. 7)                                           then
        ng = 3
        isl=2
        call ccc1_attenue (atten, cl7co2gh(1,1,ig), o3, qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
!
      else if (ib .eq. 8)                                           then
        ng = 1
        isl=2
        call ccc1_attenue (atten, cl8h2ogh(1,1,ig), o3, qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
!
      else if (ib .eq. 9)                                           then
        ng = 1
        isl=2
        call ccc1_attenue (atten, cl9h2ogh(1,1,ig), o3, qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, ng, isl, il1, il2, ilg, s1)
!
      endif
!
      return
      end
