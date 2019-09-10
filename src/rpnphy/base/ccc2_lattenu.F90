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
      subroutine ccc2_lattenu4 (atten, ib, ig, o3, qq, co2, &
                          dp, dip, dt, dt0, inpt, &
                          il1, il2, ilg)

      implicit none
!!!#include <arch_specific.hf>
!
      integer il1, il2, ilg, ib, ig, i
      real atten(ilg), o3(ilg), qq(ilg), co2(ilg), dp(ilg), dip(ilg), &
           dt(ilg), dt0(ilg)
      integer inpt(ilg)
!
!Authors
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
! 001    M. Lazare.(may 05,2006)          - previous version lattenu2 for gcm15e/f:
!                                           - pass integer variables "init" and
!                                            "nit" instead of actual integer
!                                            values, to "attenue" routines.
! 002    M. Lazare,J. Li. (apr 18,2008)          - previous version lattenu3 for gcm15g:
!                                            - calls attenue3 instead of attenue2.*
! 003    M.Lazare,K.Winger,P.Vaillancourt   (Apr 2008) - use integer variables instead of actual integers
!
! 004    J. Li  ( feb 09,2009) -  new version for gcm15h:
!                               - 3d ghg implemented, thus no need
!                                 for "trace" common block or
!                                 temporary work arrays to hold
!                                 mixing ratios of ghg depending on
!                                 a passed, specified option.
!                               - calls attenue4 instead of attenue3.
! 005   P. vaillancourt, A-M Leduc  (Apr 2010) - Integrating the Li/Lazare 2009 version. (lattenu4).
!                                           - Add co2. Remove s1,ng.  rmu becomes local . call attenu4.
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
! co2    co2 mass mixing ratio
! qq     water vapor mass mixing ratio
! dip    interpolation factor for pressure between two
!        neighboring standard input data pressure levels
! dt     layer temperature - 250 k
! dt0    temperature in moon layer - 250 k
!----------------------------------------------------------------------
!
!*
#include "ccc2_bandlh.cdk"

      integer isl
      real rmu(ilg)
!
      if (ib .eq. 1)                                                then
        isl=2

        call ccc2_attenue4 (atten, cl1co2gh(1,1,ig), co2, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, isl, il1, il2, ilg)
!
      else if (ib .eq. 2)                                           then
        do 200 i = il1, il2
          atten(i) =  0.0
  200   continue
!
      else if (ib .eq. 3)                                           then
        isl=2
        call ccc2_attenue4 (atten, cl3h2ogh(1,1,ig), qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, isl, il1, il2, ilg)
!
      else if (ib .eq. 4)                                           then
        do 400 i = il1, il2
          atten(i) =  0.0
  400   continue
!
      else if (ib .eq. 5)                                           then
        isl=2
        call ccc2_attenue4 (atten, cl5o3gh(1,1,ig), o3, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, isl, il1, il2, ilg)
!
      else if (ib .eq. 7)                                           then
        isl=2
        call ccc2_attenue4 (atten, cl7co2gh(1,1,ig), co2, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, isl, il1, il2, ilg)
!
      else if (ib .eq. 8)                                           then
        isl=2
        call ccc2_attenue4 (atten, cl8h2ogh(1,1,ig), qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, isl, il1, il2, ilg)
!
      else if (ib .eq. 9)                                           then
        isl=2
        call ccc2_attenue4 (atten, cl9h2ogh(1,1,ig), qq, dp, dip, dt, dt0, &
                      rmu, inpt, ntl, isl, il1, il2, ilg)
!
      endif
!
      return
      end
