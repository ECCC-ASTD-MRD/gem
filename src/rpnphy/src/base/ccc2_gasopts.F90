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

subroutine ccc2_gasopts5(taug, gw, dp, ib, ig, o3, qq, co2, ch4, o2, &
     inpt, mcont, dip, dt,rmu3, lev1, gh, &
     il1, il2, ilg, lay)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, ib, ig, lev1, il1, il2, mcont(ilg)
   real taug(ilg,lay)
   real gw, xx

   real dp(ilg,lay), o3(ilg,lay), qq(ilg,lay), co2(ilg,lay), &
        ch4(ilg,lay), o2(ilg,lay), dip(ilg,lay), &
        dt(ilg,lay), rmu3(ilg)
   integer inpt(ilg,lay),lc
   logical gh

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Revisions
   ! 001    M. Lazare (May 05,2006)- previous version gasopts2 for gcm15e/f:
   !                               - pass integer variables "init" and
   !                                 "mit" instead of actual integer
   !                                 values, to "tline_" routines.
   ! 002    M.Lazare,L.Solheim,J.Li (Apr 18,2008) - previous version gasopts3 for gcm15g:
   !                              cosmetic change to add threadprivate
   !                              for common block "trace", in support
   !                               - calls tline{1,2,3}y instead of
   !                                 tline{1,2,3}x.
   !                              - using temperture dependent o3 absorption
   !                                 data, adding ch4 in solar range, using
   !                                 kuruz solar function.
   ! 003    M.Lazare,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers

   ! 004    J. Li (Feb 09,2009)   -  new version for gcm15h:
   !                               - 3d ghg implemented, thus no need
   !                                 for "trace" common block or
   !                                 temporary work arrays to hold
   !                                 mixing ratios of ghg depending on
   !                                 a passed, specified option.
   !                               - calls tline{1,2,3}z instead of
   !                                 tline{1,2,3}y.
   ! 005   P. vaillancourt, A-M. Leduc (Apr 2010 - Integrating the Li/Lazare 2009 version.(gasopts4)
   !                                - add co2,ch4,o2,calculation of dto3
   !                                - remove ng,ng2,urbf
   ! 006   P. vaillancourt, (Sep 2014) - Update to gcm17/gasopts5:
   !     * MAY 01,2012 - J.LI.     NEW VERSION FOR GCM16:
   !     *                         - INCLUDE WATER VAPOUR CONTINUUM
   !     *                           (BANDS 2-4).

   !@Object
   !        Calculation of the optical depths due to nongray gaseous
   !        absorption for the solar, in each layer for a given band ib and
   !        cumulative probability gw.
   !        relative solar energy in each solar band are
   !        band 1:   622.8483
   !        band 1gh:   7.5917
   !        band 2:   430.0919
   !        band 2gh:   8.9036
   !        band 3:   238.6979
   !        band 3gh:   7.4453
   !        band 4:    33.4129
   !        band 4gh:   7.0384
   !
   !        total relative solar energy in from 0.2 - 4 um is
   !        1356.0300 w / m^2, plus 11.9096 w / m^2 in 4 - 10 um.
   !        total 1367.9396 w / m^2
   !
   !        this subroutine only calculates taug below 1 mb

   !@Arguments
   ! taug   gaseous optical depth
   ! dp     air mass path for a model layer (explained in raddriv).
   ! o3     o3 mass mixing ratio
   ! qq     water vapor mass mixing ratio
   ! co2,   ch4, o2 are all mass mixing ratios
   ! dip    interpolation factor for pressure between two
   !        neighboring standard input data pressure levels
   ! dt     layer temperature - 250 k
   ! inpt   number of the level for the standard input data pressures
   ! rmu3   a factor of solar zenith angle, given in raddriv
   !----------------------------------------------------------------------


#include "ccc2_bands.cdk"
#include "ccc_tracegases.cdk"

   integer i, k, m , init
   real dto3

   if (ib .eq. 1) then

      !----------------------------------------------------------------------
      !     band (14500 - 50000 cm^-1), nongray gaseous absorption of o3,
      !     h2o and o2.
      !     relative solar energy 630.4401 wm^-2.
      !     ig9 (50000-43000)  uvc                           1.21100 (wm^-2)
      !     ig8 (43000-37500)  uvc                           3.17570
      !     ig7 (37500-35700)  uvc                           3.20501
      !     uvc all included in gh part

      !     ig6 (35700-34200)  uvb                           4.73084
      !     ig5 (34200-32185)  uvb                           10.14919
      !     ig4 (32185-31250)  uvb, j value: 32185 cm^-1     6.70594
      !     ig3 (31250-25000)  uva                           83.43346
      !     ig2 (25000-19000)  par                         236.97212
      !     ig1 (19000-14500)  par                         280.85678
      !     par: photosynthetic active radiation
      !     note the spectral structure is slightly diff from Li & Barker
      !     (2005 JAS)
      !
      !     the effect of h2o and o2 is added with simple method
      !----------------------------------------------------------------------

      if (ig .eq. 1) then
         do k = lev1, lay
            do i = il1, il2
               m =  inpt(i,k)
               if (inpt(i,k).ge.950) m = inpt(i,k) - 1000

               if (m .lt. 7) then
                  xx      = (cs1o21 - 0.881e-05 * rmu3(i)) * o2(i,k)
               else
                  xx      = (0.108e-04 - 0.881e-05 * rmu3(i)) * o2(i,k)
               endif

               if (m .lt. 15) then
                  xx      =  xx+ (0.199e-02 - 0.952e-03 * rmu3(i)) * qq(i,k)
               else
                  xx      =  xx+ (0.208e-02 - 0.952e-03 * rmu3(i)) * qq(i,k)
               endif

               dto3      = dt(i,k) + 23.13
               taug(i,k) = ((cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                    dto3 * cs1o3(3,ig))) * o3(i,k) + xx) * dp(i,k)
            enddo
         enddo
      else
         do k = lev1, lay
            do i = il1, il2
               dto3      =  dt(i,k) + 23.13
               taug(i,k) =  (cs1o3(1,ig)+ dto3 *( cs1o3(2,ig) + &
                    dto3 * cs1o3(3,ig))) * o3(i,k) * dp(i,k)
            enddo
         enddo
      endif

      gw =  gws1(ig)

   else if (ib .eq. 2) then

      !----------------------------------------------------------------------
      !     band (8400 - 14500 cm^-1), nongray gaseous absorption of h2o,
      !     o2 and o3
      !     relative solar energy 430.0919 w m^-2
      !----------------------------------------------------------------------

      if (ig .le. 2) then
         call ccc2_tline2z(taug, cs2h2o(1,1,ig), cs2o2(1,1,ig), qq, o2, &
              dp, dip, dt, inpt, lev1, gh, mtl, &
              il1, il2, ilg, lay)
      else
         init=2
         call ccc2_tline1z(taug, cs2h2o(1,1,ig), qq,  dp, dip, &
              dt, inpt, lev1, gh, mtl, init, &
              il1, il2, ilg, lay)
      endif

      !----------------------------------------------------------------------
      !     simply add o3 effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            taug(i,k)   =  taug(i,k) + cs2o3(ig) * o3(i,k) * dp(i,k)
         enddo
      enddo


      !----------------------------------------------------------------------C
      !     WATER VAPOUR CONTINUUM                                           C
      !----------------------------------------------------------------------C

      LC =  5
      call ccc2_TCONTL1(TAUG, CS2CS(1,1,IG), CS2CF(1,1,IG), QQ, DP, DIP, &
           DT, LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)

      gw =  gws2(ig)

   else if (ib .eq. 3) then

      !----------------------------------------------------------------------
      !     band (4200 - 8400 cm^-1), nongray gaseous absorption of h2o
      !     and ch4
      !     relative solar energy  238.6979 wm^-2
      !----------------------------------------------------------------------

      if (ig .le. 2) then
         call ccc2_tline3z(taug, cs3h2o(1,1,ig), cs3co2(1,1,ig), &
              cs3ch4(1,1,ig), qq, co2,ch4, dp, dip, dt, inpt, &
              lev1, gh, mtl, il1, il2, ilg, lay)
      else
         call ccc2_tline2z(taug, cs3h2o(1,1,ig), cs3co2(1,1,ig), qq, co2, &
              dp, dip, dt, inpt, lev1, gh, mtl, &
              il1, il2, ilg, lay)
      endif

      !----------------------------------------------------------------------C
      !     WATER VAPOUR CONTINUUM                                           C
      !----------------------------------------------------------------------C

      LC =  5
      call ccc2_TCONTL1(TAUG, CS3CS(1,1,IG), CS3CF(1,1,IG), QQ, DP, DIP, &
           DT, LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)

      gw =  gws3(ig)

   else if (ib .eq. 4) then

      !----------------------------------------------------------------------
      !     band (2500 - 4200 cm^-1), nongray gaseous absorption of h2o
      !     and co2
      !     relative solar energy 33.4129 w m^-2
      !----------------------------------------------------------------------

      if (ig .le.2) then
         call ccc2_tline3z(taug, cs4h2o(1,1,ig), cs4co2(1,1,ig), &
              cs4ch4(1,1,ig), qq, co2,ch4,dp, dip, dt, inpt, &
              lev1, gh, mtl, il1, il2, ilg,lay)
      else
         call ccc2_tline2z(taug, cs4h2o(1,1,ig), cs4co2(1,1,ig), qq,co2, &
              dp, dip, dt, inpt, lev1, gh, mtl, &
              il1, il2, ilg,lay)
      endif

      !----------------------------------------------------------------------C
      !     WATER VAPOUR CONTINUUM                                           C
      !----------------------------------------------------------------------C

      LC =  5
      call ccc2_TCONTL1(TAUG, CS4CS(1,1,IG), CS4CF(1,1,IG), QQ, DP, DIP, &
           DT, LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)

      gw =  gws4(ig)

   endif

   return
end subroutine ccc2_gasopts5
