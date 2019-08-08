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
!-------------------------------------- LICENCE END --------------------------

subroutine ccc2_gasoptl7(taug, gw, dp, ib, ig, o3, qq, co2, ch4,an2o, &
     f11,f12,f113,f114, inptr, inpt, mcont, &
     dir, dip, dt, lev1, gh, &
     il1, il2, ilg, lay)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, ib, ig, mcont(ilg), lev1, il1, il2
   real taug(ilg,lay), gw

   real dp(ilg,lay), o3(ilg,lay), qq(ilg,lay), co2(ilg,lay), &
        ch4(ilg,lay),an2o(ilg,lay),f11(ilg,lay), f12(ilg,lay), &
        f113(ilg,lay), f114(ilg,lay), dir(ilg,lay), &
        dip(ilg,lay), dt(ilg,lay)
   real CAN2O,CRMCO2
   integer inptr(ilg,lay), inpt(ilg,lay)
   logical gh
   integer init2

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !         JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Revisions
   ! 001    M.Lazare(May 05,2006)  - previous version gasoptl3 for gcm15e:
   !                               - pass integer variables "init" and
   !                                 "mit" instead of actual integer
   !                                 values, to "tline_" routines.
   ! 002    J. Li(Jun 21,2006)     - previous version gasoptl4 for gcm15f:
   !                               - update the minor gas absorption.
   !                               - calls new tconthl2.
   ! 003    M.Lazare,L.Solheim,J.Li (Apr 18,2008) - previous version gasoptl5 for gcm15g:
   !                               - cosmetic change to add threadprivate
   !                                 for common block "trace", in support
   !                               - calls tline{1,2,3}y instead of
   !                                 tline{1,2,3}x. also calls tlinehc3,
   !                                 tconthl3 instead of tlinehc2,tconthl2.
   ! 004    K.Winger,P.Vaillancourt  (Apr 08) - integer mtl from bandl.cdk
   !
   ! 005    J.Li (Feb 09,2009)     - new version for gcm15h:
   !                               - 3d ghg implemented, thus no need
   !                                 for "trace" common block or
   !                                 temporary work arrays to hold
   !                                 mixing ratios of ghg depending on
   !                                 a passed, specified option.
   !                               - calls tline{1,2,3}z instead of
   !                                 tline{1,2,3}y. also calls tlinehc4,
   !                                 tconthl4 instead of tlinehc3,tconthl3.
   ! 006   P. vaillancourt, A-M. Leduc (Apr 2010 - Integrating the Li/Lazare 2009 version.(gasoptl6)
   !                                - add co2,ch4,ano2,f11,f12,f113,f114
   !
   ! 007    P.Vaillancourt   (Sep 14) - update to gcm17/gasoptl7.f:
   !     * APR 21/2010 - J.LI.     NEW VERSION FOR GCM15I:
   !     *                         - REVISION OF MANY BANDS FOR GREATER
   !     *                           ACCURACY.

   !@bject
   !        Calculation of the optical depths due to nongray gaseous
   !        absorption for the infrared, in each layer for a given band ib
   !        and cumulative probability gw.
   !        From band1 to band4, the solar and infrared interaction is
   !        considered. the total solar energy considered in the infrared
   !        region is 11.9006 w / m^2
   !        for gases with constant mixing ratio:
   !        1 = h2o
   !        2 = o3
   !        3 = co2
   !        4 = ch4
   !        5 = n2o
   !        6 = o2
   !        7 = cfc11
   !        8 = cfc12
   !        tline, etc., deal with line absorption and tcontl and tconthl
   !        deal with water vapor continuum

   !@Arguments
   !
   ! taug   gaseous optical depth
   ! dp     air mass path for a model layer (explained in raddriv).
   ! o3     o3 mass mixing ratio
   ! qq     water vapor mass mixing ratio
   ! an2o   n2o, also co2,ch4,f11,f12,f113,f114 are mass mixing ratios
   ! dir    interpolation factor for mass ratio of h2o / co2
   !        between two neighboring standard input ratios
   ! dip    interpolation factor for pressure between two
   !        neighboring standard input data pressure levels
   ! dt     layer temperature - 250 k
   ! inpr   number of the ratio level for the standard 5 ratios
   ! inpt   number of the level for the standard input data pressures
   ! mcont  the highest level for water vapor continuum calculation
   !----------------------------------------------------------------------

#include "ccc_tracegases.cdk"
#include "ccc2_bandl.cdk"

   !     number of vertical levels in absorber pressure-based coefficient
   !     array ("m" references non-saturated bands active below 1 mb only).

   integer k, i, lc
   real fact

   data init2 /2/

   if (ib .eq. 1) then

      !----------------------------------------------------------------------
      !     band (2500 - 2200 cm^-1), nongray gaseous absorption of h2o and
      !     co2.
      !----------------------------------------------------------------------

      call ccc2_tline2z(taug, cl1h2o, cl1co2, qq, co2, &
           dp, dip, dt, inpt, &
           lev1, gh, mtl, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the n2o effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            CAN2O     =  2.0437 * 1.E-07 * sqrt(AN2O(I,K) * 1.E+07)
            FACT      =  QQ(I,K) / (QQ(I,K) + 8.E+04 * CAN2O)
            TAUG(I,K) =  TAUG(I,K) + (754.9786 + 10141.5049 * FACT * FACT) * &
                 CAN2O * DP(I,K)
         enddo
      enddo

      gw =  gwl1(ig)

   else if (ib .eq. 2) then

      !----------------------------------------------------------------------
      !     band (2200 - 1900 cm^-1), nongray gaseous absorption of h2o + n2o
      !----------------------------------------------------------------------

      call ccc2_tline1z(taug, cl2h2o, qq, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay)

      lc =  3
      call ccc2_tcontl1(taug, cl2cs, cl2cf, qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the n2o effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            CAN2O     =  2.0437 * 1.E-07 * sqrt(AN2O(I,K) * 1.E+07)
            FACT      =  QQ(I,K) / (QQ(I,K) + 72000. * CAN2O)
            TAUG(I,K) =  TAUG(I,K) + (93. + 3500. * FACT * FACT) * &
                 CAN2O * DP(I,K)
         enddo
      enddo

      gw =  gwl2(ig)

   else if (ib .eq. 3) then

      !----------------------------------------------------------------------
      !     band (1900 - 1400 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      call ccc2_tline1z(taug, cl3h2o(1,1,ig), qq, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay)

      lc =  4
      call ccc2_tcontl1(taug, cl3cs(1,1,ig), cl3cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      do K = LEV1, LAY
         do I = IL1, IL2
            TAUG(I,K) =  TAUG(I,K) + CL3CH4(IG) * CH4(I,K) * DP(I,K)
         enddo
      enddo

      gw =  gwl3(ig)

   else if (ib .eq. 4) then

      !----------------------------------------------------------------------
      !     band3 (1100 - 1400 cm^-1), overlapping absorption of h2o, n2o,
      !     ch4 and cfc12. direct mapping method for h2o and ch4 and n2o
      !     cfc are considered as minor gases
      !----------------------------------------------------------------------

      call ccc2_tline3z(taug, cl4h2o(1,1,ig), cl4ch4(1,1,ig), &
           cl4n2o(1,1,ig), qq, ch4, an2o, dp, dip, &
           dt, inpt,lev1, gh, mtl, il1, il2, &
           ilg, lay)

      lc =  4
      call ccc2_tcontl1(taug, cl4cs(1,1,ig), cl4cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the cfc effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            TAUG(I,K) =  TAUG(I,K) + (CL4F12(IG) * F12(I,K) +       &
                 1037.3 * F113(I,K) + 1426.9 * F114(I,K)) * &
                 DP(I,K)
         enddo
      enddo

      gw =  gwl4(ig)

   else if (ib .eq. 5) then

      !----------------------------------------------------------------------
      !     band5 (980 - 1100 cm^-1), overlapping absorption of h2o and o3
      !     direct mapping method. co2 and cfc are simply added
      !----------------------------------------------------------------------
      !
      call ccc2_tline2z(taug, cl5h2o(1,1,ig), cl5o3(1,1,ig), qq, o3, &
           dp, dip, dt, inpt, lev1, gh, mtl, &
           il1, il2, ilg, lay)

      lc =  4
      call ccc2_tcontl1(taug, cl5cs(1,1,ig), cl5cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the co2 + cfc effect
      !     since the interaction of co2 and h2o, qq(i,k) appears in co2
      !     effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            CRMCO2    =  2.3056E-04 * sqrt(CO2(I,K) * 1.E+04)
            TAUG(I,K) =  TAUG(I,K) + ( (0.009 +  0.093 * QQ(I,K) / (QQ(I,K) +   &
                 2.1 * CRMCO2)) * CO2(I,K) + CL5F11(IG) * F11(I,K) +    &
                 CL5F12(IG) * F12(I,K) + 1687.4 * F113(I,K) +           &
                 2924.1 * F114(I,K)) * DP(I,K)
         enddo
      enddo

      gw =  gwl5(ig)

   else if (ib .eq. 6) then

      !----------------------------------------------------------------------
      !     band (800 - 980 cm^-1), nongray gaseous absorption of h2o.
      !     + cfc11 and cfc12
      !----------------------------------------------------------------------

      call ccc2_TLINE3Z(TAUG, CL6H2O(1,1,IG), CL6F11(1,1,IG),           &
           CL6F12(1,1,IG), QQ, F11, F12, DP, DIP, DT, INPT, &
           LEV1, GH, MTL, IL1, IL2, ILG, LAY)

      !----------------------------------------------------------------------
      !     simply add the co2 + cfc effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            TAUG(I,K) =  TAUG(I,K) + ( (0.0074 + 0.0396 * QQ(I,K) /   &
                 (QQ(I,K) + 2.8 * CO2(I,K))) * CO2(I,K) +      &
                 1191. * F113(I,K) + 1098. * F114(I,K) ) * DP(I,K)
         enddo
      enddo

      if (ig .eq. 1) then
         lc =  4
         call ccc2_tcontl1 (taug, cl6cs, cl6cf, qq, dp, dip, dt, &
              lc, inpt, mcont, gh, il1, il2, ilg, lay)
      endif

      gw =  gwl6(ig)

   else if (ib .eq. 7) then

      !----------------------------------------------------------------------
      !     band6 (540 - 800 cm^-1), overlapping absorption of h2o and co2
      !     exact mapping method for h2o and co2, direct mapping for n2o
      !     o3 effect is simply added
      !----------------------------------------------------------------------

      call ccc2_tlinehc4(taug, cl7h2ou(1,1,ig), cl7h2od(1,1,1,ig), &
           cl7co2u(1,1,ig), cl7co2d(1,1,1,ig), qq, co2, dp, dip, &
           dir, dt, inptr, inpt, lev1, il1, il2, ilg, lay)
      !
      call ccc2_tconthl4(taug, cl7cs(1,1,1,ig), cl7cf(1,1,1,ig), qq, dp, dip, &
           dir, dt, inptr, inpt, mcont, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the o3 and n2o effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            CAN2O     =  2.0437 * 1.E-07 * sqrt(AN2O(I,K) * 1.E+07)
            TAUG(I,K) =  TAUG(I,K) + (CL7O3(IG) * O3(I,K) +    &
                 CL7N2O(IG) * CAN2O) * DP(I,K)
         enddo
      enddo

      gw =  gwl7(ig)

   else if (ib .eq. 8) then

      !----------------------------------------------------------------------
      !     band (340 - 540 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------
      !
      call ccc2_tline1z(taug, cl8h2o(1,1,ig), qq, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay)

      if (ig .le. 4) then
         lc =  6
         call ccc2_tcontl1(taug, cl8cs(1,1,ig), cl8cf(1,1,ig), qq, dp, dip,dt, &
              lc, inpt, mcont, gh, il1, il2, ilg, lay)
      endif

      gw =  gwl8(ig)

   else if (ib .eq. 9) then

      !----------------------------------------------------------------------
      !     band (0 - 340 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      call ccc2_tline1z(taug, cl9h2o(1,1,ig), qq, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay)

      lc =  6
      call ccc2_tcontl1(taug, cl9cs(1,1,ig), cl9cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      gw =  gwl9(ig)

   endif

   return
end subroutine ccc2_gasoptl7
