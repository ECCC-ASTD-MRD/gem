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

subroutine ccc2_gasoptlgh7(taug, gwgh, dp, ib, ig, o3, qq, &
     co2, ch4, an2o, inpt, mcont, &
     dip, dt, lev1, gh, &
     il1, il2, ilg, lay)

   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, ib, ig, mcont(ilg), lev1, il1, il2
   real taug(ilg,lay), gwgh

   real dp(ilg,lay), o3(ilg,lay), qq(ilg,lay), co2(ilg,lay), &
        ch4(ilg,lay), an2o(ilg,lay), dip(ilg,lay), &
        dt(ilg,lay)
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
   ! 001    M. Lazare (May 05/2006) - previous version gasoptlgh3 for gcm15e/f:
   !                                - calls new versions of:
   !                                  tline1x,tline2x,tline3y,tcontl1,
   !                                  tconthl1.
   ! 002    M. Lazare, J. Li (Apr 18,2008) - previous version gasoptlgh4 for gcm15g:
   !                                - cosmetic change to use scalar variable
   !                                 "initaug" (=2) in calls to tline1y
   !                                 instead of the actual number itself.
   !                                 similar cosmetic change to use ntl(=28)
   !                                 scalar variable instead of the actual
   !                                 number itself in calls to all "tline_"
   !                                 routines.
   !                               - calls tline{1,2,3}y instead of tline{1,2,3}x.
   ! 003    M.Lazarre,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers
   ! 004    J.Li (Feb 09,2009)     - new version for gcm15h:
   !                               - 3d ghg implemented, thus no need
   !                                 for "trace" common block or
   !                                 temporary work arrays to hold
   !                                 mixing ratios of ghg depending on
   !                                 a passed, specified option.
   !                               - calls tline{1,2,3}z instead of
   !                                 tline{1,2,3}y.
   ! 005   P. vaillancourt, A-M. Leduc (Apr 2010 - Integrating the Li/Lazare 2009 version.(gasoptlgh5)
   !                                 - add co2,ch4,an2o
   ! 006   P. vaillancourt (Sep 2014) - Update to gcm17/gasoptlgh6: - only change is in content of bandlh.cdk
   !     * APR 22/2010 - J.LI.     NEW VERSION FOR GCM15I:
   !     *                         - ADD ONE EXTRA TERM TO BANDL4GH FOR
   !     *                           GREATER ACCURACY.
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
   ! an2o  n2o, also co2,ch4 and o3 are mass mixing ratios
   ! qq    water vapor mass mixing ratio
   ! dip   interpolation factor for pressure between two
   !       neighboring standard input data pressure levels
   ! dt    layer temperature - 250 k
   ! inpt  number of the level for the standard input data pressures
   !----------------------------------------------------------------------

#include "ccc2_bandlh.cdk"

   !     * initaug is a switch used in tline1y (as "iplus") which
   !     * initializes taug to zero if its value is two. this is what
   !     * we require throughout this routine.

   integer  i, k , lc

   data initaug /2/


   if (ib .eq. 1)                                                then

      !----------------------------------------------------------------------
      !     band (2500 - 2200 cm^-1), nongray gaseous absorption of co2.
      !----------------------------------------------------------------------

      call ccc2_tline1z (taug, cl1co2gh(1,1,ig), co2, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay)

      gwgh =  gwl1gh(ig)

   else if (ib .eq. 2)                                           then

      !----------------------------------------------------------------------
      !     band (2200 - 1900 cm^-1), nongray gaseous absorption of h2o
      !----------------------------------------------------------------------

      call ccc2_tline1z (taug, cl2h2ogh, qq, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay)

      lc =  3
      call ccc2_tcontl1 (taug, cl2csgh, cl2cfgh, qq, dp, dip, dt, &
           lc, inpt, mcont, gh, il1, il2, ilg, lay)

      gwgh =  gwl2gh(ig)

   else if (ib .eq. 3)                                           then

      !----------------------------------------------------------------------
      !     band (1900 - 1400 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      call ccc2_tline1z (taug, cl3h2ogh(1,1,ig), qq, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay)

      if (ig .eq. 1)                                                then
         lc =  4
         call ccc2_tcontl1 (taug, cl3csgh, cl3cfgh, qq, dp, dip, dt, &
              lc, inpt, mcont, gh, il1, il2, ilg, lay)

      endif

      gwgh =  gwl3gh(ig)

   else if (ib .eq. 4)                                           then

      !----------------------------------------------------------------------
      !     band3 (1100 - 1400 cm^-1), overlapping absorption of h2o, n2o,
      !     and ch4. direct mapping method for h2o and ch4 and n2o
      !----------------------------------------------------------------------

      call ccc2_tline3z (taug,cl4h2ogh(1,1,ig),cl4ch4gh(1,1,ig), &
           cl4n2ogh(1,1,ig), qq, ch4, an2o, dp, dip, dt, inpt, &
           lev1, gh, ntl, il1, il2, ilg, lay)

      gwgh =  gwl4gh(ig)

   else if (ib .eq. 5)                                           then

      !----------------------------------------------------------------------
      !     band5 (980 - 1100 cm^-1), overlapping absorption of h2o and o3
      !     direct mapping method
      !----------------------------------------------------------------------

      call ccc2_tline2z (taug, cl5h2ogh(1,1,ig), cl5o3gh(1,1,ig), qq, o3, &
           dp, dip, dt, inpt, &
           lev1, gh, ntl, il1, il2, ilg, lay)

      if (ig .le. 2)                                                then
         lc =  4
         call ccc2_tcontl1 (taug, cl5csgh(1,1,ig), cl5cfgh(1,1,ig), qq, dp, dip, &
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
         call ccc2_tline2z (taug, cl7h2ogh(1,1,ig), cl7co2gh(1,1,ig), qq, co2, &
              dp, dip, dt, inpt, &
              lev1, gh, ntl, il1, il2, ilg, lay)

         !----------------------------------------------------------------------
         !     simply add the o3 effect
         !----------------------------------------------------------------------

         if (ig .le. 2)                                              then
            do k = 1, lay
               do i = il1, il2
                  taug(i,k) =  taug(i,k) + cl7o3gh(ig) * o3(i,k) * dp(i,k)
               enddo
            enddo
         endif
      else

         call ccc2_tline1z (taug, cl7co2gh(1,1,ig), co2, dp, dip, &
              dt, inpt, lev1, gh, ntl, initaug, &
              il1, il2, ilg, lay)
      endif

      gwgh =  gwl7gh(ig)

   else if (ib .eq. 8)                                           then

      !----------------------------------------------------------------------
      !     band (340 - 540 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      call ccc2_tline1z (taug, cl8h2ogh(1,1,ig), qq, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay)

      gwgh =  gwl8gh(ig)

   else if (ib .eq. 9)                                           then

      !----------------------------------------------------------------------
      !     band (0 - 340 cm^-1), nongray gaseous absorption of h2o.
      !----------------------------------------------------------------------

      call ccc2_tline1z (taug, cl9h2ogh(1,1,ig), qq, dp, dip, &
           dt, inpt, lev1, gh, ntl, initaug, &
           il1, il2, ilg, lay)

      gwgh =  gwl9gh(ig)

   endif

   return
end subroutine ccc2_gasoptlgh7
