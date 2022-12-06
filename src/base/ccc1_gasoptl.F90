!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

subroutine ccc1_gasoptl2(taug, gw, dp, ib, ig, &
     o3, qq, inptr, inpt, mcont, &
     dir, dip, dt, lev1, gh, &
     il1, il2, ilg, lay, tg)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, ib, ig, lev1, il1, il2, ng2
   integer k, i, ng, lc, ng3
   integer mcont(ilg)
   real taug(ilg,lay), gw, fact

   real dp(ilg,lay), o3(ilg,lay), qq(ilg,lay), dir(ilg,lay), &
        dip(ilg,lay), dt(ilg,lay), tg(ilg,lay)
   integer inptr(ilg,lay), inpt(ilg,lay)
   logical gh
   integer init1,init2

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !         JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Revisions
   ! 001    M. Lazare (May 2005) new version for gcm15e: pass integer
   !        variables "init" and "mtl" instead of actual integer values,
   !        to "tline_" routines.
   ! 002    J.Li (June 2006) update the minor gas absorption
   ! 003    K.Winger,P.Vaillancourt   (Apr 08) - integer mtl from bandl.cdk

   !@Object
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
#include "ccc1_bandl.cdk"

   !     number of vertical levels in absorber pressure-based coefficient
   !     array ("m" references non-saturated bands active below 1 mb only).

   data init1,init2 /1,2/

   if (ib .eq. 1)                                                then

      !----------------------------------------------------------------------
      !     band (2500 - 2200 cm^-1), nongray gaseous absorption of h2o and
      !     co2.
      !----------------------------------------------------------------------

      ng2 =  3
      call ccc1_tline2(taug, cl1h2o, cl1co2, qq, o3, &
           ng2, dp, dip, dt, inpt, &
           lev1, gh, mtl, il1, il2, ilg, lay, tg)
      
      !----------------------------------------------------------------------
      !     simply add the n2o effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            fact      =  qq(i,k) / (qq(i,k) + 8.e+04 * rmn2o)
            taug(i,k) =  taug(i,k) + (754.9786 + 10141.5049 * fact * fact) * &
                 rmn2o * dp(i,k)
         enddo
      enddo

      gw =  gwl1(ig)

   else if (ib .eq. 2)                                           then

      !----------------------------------------------------------------------
      !     band (2200 - 1900 cm^-1), nongray gaseous absorption of h2o + n2o
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1(taug, cl2h2o, qq, ng, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay, tg)

      lc =  3
      call ccc1_tcontl2(taug, cl2cs, cl2cf, qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the n2o effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            fact      =  qq(i,k) / (qq(i,k) + 72000. * rmn2o)
            taug(i,k) =  taug(i,k) + (93. + 3500. * fact * fact) * rmn2o * &
                 dp(i,k)
         enddo
      enddo

      gw =  gwl2(ig)

   else if (ib .eq. 3)                                           then

      !----------------------------------------------------------------------
      !     band (1900 - 1400 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1(taug, cl3h2o(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay, tg)

      lc =  4
      call ccc1_tcontl2(taug, cl3cs(1,1,ig), cl3cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      gw =  gwl3(ig)

   else if (ib .eq. 4)                                           then

      !----------------------------------------------------------------------
      !     band3 (1100 - 1400 cm^-1), overlapping absorption of h2o, n2o,
      !     ch4 and cfc12. direct mapping method for h2o and ch4 and n2o
      !     cfc are considered as minor gases
      !----------------------------------------------------------------------

      ng2 =  4
      ng3 =  5
      call ccc1_tline3(taug, cl4h2o(1,1,ig), cl4ch4(1,1,ig), &
           cl4n2o(1,1,ig), qq, ng2, ng3, dp, dip, dt, inpt, &
           lev1, gh, mtl, il1, il2, ilg, lay)

      lc =  4
      call ccc1_tcontl2(taug, cl4cs(1,1,ig), cl4cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the cfc effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            taug(i,k) =  taug(i,k) + (cl4f12(ig) * rmf12 + 629.0 * rmf113 + &
                 751.5 * rmf114) * dp(i,k)
         enddo
      enddo

      gw =  gwl4(ig)

   else if (ib .eq. 5)                                           then

      !----------------------------------------------------------------------
      !     band5 (980 - 1100 cm^-1), overlapping absorption of h2o and o3
      !     direct mapping method. co2 and cfc are simply added
      !----------------------------------------------------------------------

      ng2 =  2
      call ccc1_tline2(taug, cl5h2o(1,1,ig), cl5o3(1,1,ig), qq, o3, &
           ng2, dp, dip, dt, inpt, &
           lev1, gh, mtl, il1, il2, ilg, lay, tg)

      lc =  4
      call ccc1_tcontl2(taug, cl5cs(1,1,ig), cl5cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      !----------------------------------------------------------------------
      !     simply add the co2 + cfc effect
      !     since the interaction of co2 and h2o, qq(i,k) appears in co2
      !     effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            taug(i,k) =  taug(i,k) + ((0.009 +  0.093 * qq(i,k) / (qq(i,k) + &
                 2.1 * rmco2)) * rmco2 + &
                 cl5f11(ig) * rmf11 + cl5f12(ig) * rmf12 + &
                 1023.0 * rmf113 + 1539.0 * rmf114) * dp(i,k)
         enddo
      enddo

      gw =  gwl5(ig)

   else if (ib .eq. 6)                                           then

      !----------------------------------------------------------------------
      !     band (800 - 980 cm^-1), nongray gaseous absorption of h2o.
      !     + cfc11 and cfc12
      !----------------------------------------------------------------------

      if (ig .eq. 1)                                                then
         ng2 =  7
         ng3 =  8
         call ccc1_tline3(taug, cl6h2o(1,1,ig), cl6f11, &
              cl6f12, qq, ng2, ng3, dp, dip, dt, inpt, &
              lev1, gh, mtl, il1, il2, ilg, lay)

         lc =  4
         call ccc1_tcontl2(taug, cl6cs, cl6cf, qq, dp, dip, dt, &
              lc, inpt, mcont, gh, il1, il2, ilg, lay)

         !----------------------------------------------------------------------
         !     simply add the co2 + cfc effect
         !----------------------------------------------------------------------

         do k = lev1, lay
            do i = il1, il2
               taug(i,k) =  taug(i,k) + ( (0.0074 + 0.0396 * qq(i,k) / &
                    (qq(i,k) + 2.8 * rmco2)) * rmco2 + &
                    722.0 * rmf113 + 578.0 * rmf114 ) * dp(i,k)
            enddo
         enddo
      else
         ng =  1
         call ccc1_tline1(taug, cl6h2o(1,1,ig), qq, ng, dp, dip, &
              dt, inpt, lev1, gh, mtl, init2, &
              il1, il2, ilg, lay, tg)
      endif

      gw =  gwl6(ig)

   else if (ib .eq. 7)                                           then

      !----------------------------------------------------------------------
      !     band6 (540 - 800 cm^-1), overlapping absorption of h2o and co2
      !     exact mapping method for h2o and co2, direct mapping for n2o
      !     o3 effect is simply added
      !----------------------------------------------------------------------

      call ccc1_tlinehc(taug, cl7h2ou(1,1,ig), cl7h2od(1,1,1,ig), &
           cl7co2u(1,1,ig), cl7co2d(1,1,1,ig), qq, dp, dip,dir, &
           dt, inptr, inpt, lev1, il1, il2, ilg, lay)

      call ccc1_tconthl3(taug, cl7cs(1,1,1,ig), cl7cf(1,1,1,ig), qq, dp, dip, &
           dir, dt, inptr, inpt, mcont, il1, il2, ilg, lay)

      ng =  5
      call ccc1_tline1(taug, cl7n2o(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, mtl, init1, &
           il1, il2, ilg, lay, tg)

      !----------------------------------------------------------------------
      !     simply add the o3 effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            taug(i,k) =  taug(i,k) + cl7o3(ig) * o3(i,k) * dp(i,k)
         enddo
      enddo

      gw =  gwl7(ig)

   else if (ib .eq. 8)                                           then

      !----------------------------------------------------------------------
      !     band (340 - 540 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1(taug, cl8h2o(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay, tg)

      if (ig .le. 4)                                                then
         lc =  6
         call ccc1_tcontl2(taug, cl8cs(1,1,ig), cl8cf(1,1,ig), qq, dp, dip,dt, &
              lc, inpt, mcont, gh, il1, il2, ilg, lay)
      endif

      gw =  gwl8(ig)

   else if (ib .eq. 9)                                           then

      !----------------------------------------------------------------------
      !     band (0 - 340 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1(taug, cl9h2o(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, mtl, init2, &
           il1, il2, ilg, lay, tg)

      lc =  6
      call ccc1_tcontl2(taug, cl9cs(1,1,ig), cl9cf(1,1,ig), qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      gw =  gwl9(ig)

   endif

   return
end subroutine ccc1_gasoptl2
