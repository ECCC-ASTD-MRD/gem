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

subroutine ccc1_gasoptlgh3(taug, gwgh, dp, ib, ig, &
     o3, qq, inpt, mcont, &
     dip, dt, lev1, gh, &
     il1, il2, ilg, lay, tg)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, ib, ig, lev1, il1, il2, ng, lc
   integer ng2, ng3, k, i
   integer mcont(ilg)
   real taug(ilg,lay), gwgh

   real dp(ilg,lay), o3(ilg,lay), qq(ilg,lay), &
        dip(ilg,lay), dt(ilg,lay), tg(ilg,lay)
   integer inpt(ilg,lay)
   logical gh

   integer initaug

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)
   !@Revisions
   ! 001    M.Lazarre,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers
   !@Object OPTICAL DEPTHS CALCULATION
   !        The same as gasoptl but for intervals close to 1 in the
   !        accumulated probability space
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
   ! taug  gaseous optical depth
   ! dp    air mass path for a model layer (explained in raddriv).
   ! o3    o3 mass mixing ratio
   ! qq    water vapor mass mixing ratio
   ! dip   interpolation factor for pressure between two
   !       neighboring standard input data pressure levels
   ! dt    layer temperature - 250 k
   ! inpt  number of the level for the standard input data pressures
   !----------------------------------------------------------------------

#include "ccc1_bandlh.cdk"

   !     * initaug is a switch used in tline1y (as "iplus") which
   !     * initializes taug to zero if its value is two. this is what
   !     * we require throughout this routine.

   data initaug /2/


   if (ib .eq. 1)                                                then

      !----------------------------------------------------------------------
      !     band (2500 - 2200 cm^-1), nongray gaseous absorption of co2.
      !----------------------------------------------------------------------

      ng =  3
      call ccc1_tline1 (taug, cl1co2gh(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay, tg)

      gwgh =  gwl1gh(ig)

   else if (ib .eq. 2)                                           then

      !----------------------------------------------------------------------
      !     band (2200 - 1900 cm^-1), nongray gaseous absorption of h2o
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1 (taug, cl2h2ogh, qq, ng, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay, tg)

      lc =  3
      call ccc1_tcontl2 (taug, cl2csgh, cl2cfgh, qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      gwgh =  gwl2gh(ig)

   else if (ib .eq. 3)                                           then

      !----------------------------------------------------------------------
      !     band (1900 - 1400 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1 (taug, cl3h2ogh(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay, tg)

      if (ig .eq. 1)                                                then
         lc =  4
         call ccc1_tcontl2 (taug, cl3csgh, cl3cfgh, qq, dp, dip, dt, &
              lc, inpt, mcont, gh, il1, il2, ilg, lay)

      endif

      gwgh =  gwl3gh(ig)

   else if (ib .eq. 4)                                           then

      !----------------------------------------------------------------------
      !     band3 (1100 - 1400 cm^-1), overlapping absorption of h2o, n2o,
      !     and ch4. direct mapping method for h2o and ch4 and n2o
      !----------------------------------------------------------------------

      ng2 =  4
      ng3 =  5
      call ccc1_tline3 (taug,cl4h2ogh(1,1,ig),cl4ch4gh(1,1,ig), &
           cl4n2ogh(1,1,ig), qq, ng2, ng3, dp, dip, dt, inpt, &
           lev1, gh, ntl, il1, il2, ilg, lay)

      gwgh =  gwl4gh(ig)

   else if (ib .eq. 5)                                           then

      !----------------------------------------------------------------------
      !     band5 (980 - 1100 cm^-1), overlapping absorption of h2o and o3
      !     direct mapping method
      !----------------------------------------------------------------------

      ng2 =  2
      call ccc1_tline2 (taug, cl5h2ogh(1,1,ig), cl5o3gh(1,1,ig), qq, o3, &
           ng2, dp, dip, dt, inpt, &
           lev1, gh, ntl, il1, il2, ilg, lay, tg)

      if (ig .le. 2)                                                then
         lc =  4
         call ccc1_tcontl2 (taug, cl5csgh(1,1,ig), cl5cfgh(1,1,ig), qq, dp, dip, &
              dt, lc, inpt, mcont, gh, il1, il2, ilg, lay)
      endif

      gwgh =  gwl5gh(ig)

      !----------------------------------------------------------------------
      !     band (800 - 980 cm^-1), no gh
      !----------------------------------------------------------------------

   else if (ib .eq. 7)                                           then

      !----------------------------------------------------------------------
      !     band6 (540 - 800 cm^-1), overlapping absorption of h2o and co2
      !     direct mapping method. for ig > 4, the contribution by h2o is
      !     very small.
      !----------------------------------------------------------------------

      if (ig .le. 4)                                                then
         ng2 =  3
         call ccc1_tline2 (taug, cl7h2ogh(1,1,ig), cl7co2gh(1,1,ig), qq, o3, &
              ng2, dp, dip, dt, inpt, &
              lev1, gh, ntl, il1, il2, ilg, lay, tg)

         !----------------------------------------------------------------------
         !     simply add the o3 effect
         !----------------------------------------------------------------------

         if (ig .le. 2)                                              then
            do  k = 1, lay
               do  i = il1, il2
                  taug(i,k) =  taug(i,k) + cl7o3gh(ig) * o3(i,k) * dp(i,k)
               enddo
            enddo
         endif
      else

         ng =  3
         call ccc1_tline1 (taug, cl7co2gh(1,1,ig), qq, ng, dp, dip, &
              dt, inpt, lev1, gh, ntl, initaug, &
              il1, il2, ilg, lay, tg)
      endif

      gwgh =  gwl7gh(ig)

   else if (ib .eq. 8)                                           then

      !----------------------------------------------------------------------
      !     band (340 - 540 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1 (taug, cl8h2ogh(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay, tg)

      gwgh =  gwl8gh(ig)

   else if (ib .eq. 9)                                           then

      !----------------------------------------------------------------------
      !     band (0 - 340 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      ng =  1
      call ccc1_tline1 (taug, cl9h2ogh(1,1,ig), qq, ng, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay, tg)

      gwgh =  gwl9gh(ig)

   endif

   return
end subroutine ccc1_gasoptlgh3
