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

subroutine ccc1_gasopts(taug, gw, dp, ib, ig, o3, qq, inpt, dip, dt, &
     rmu3, lev1, gh, il1, il2, ilg, lay, urbf)
   implicit none
!!!#include <arch_specific.hf>
 
   integer ilg, lay, ib, ig, lev1, il1, il2, k, i, m,ng2, ng
   real taug(ilg,lay)
   real gw, xx
 
   real dp(ilg,lay), o3(ilg,lay), qq(ilg,lay), dip(ilg,lay), &
        dt(ilg,lay), rmu3(ilg), urbf(ilg,lay)
   integer inpt(ilg,lay)
   logical gh
   integer init
 
   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)
 
   !@Revisions
   ! 001    M.Lazarre,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers

   !@Object
   !        Calculation of the optical depths due to nongray gaseous
   !        absorption for the solar, in each layer for a given band ib and
   !        cumulative probability gw.
   !        relative solar energy in each solar band are
   !        band 1:   621.6958
   !        band 1gh:   6.8988
   !        band 2:   429.5020
   !        band 2gh:   8.7245
   !        band 3:   242.8750
   !        band 3gh:   4.0330
   !        band 4:    30.9570
   !        band 4gh:   9.6167
   !
   !        total relative solar energy in from 0.2 - 4 um is
   !        1354.3029 w / m^2, plus 11.9006 w / m^2 in 4 - 10 um.
   !        total 1366.2035 w / m^2
   !
   !        minor gas:
   !        3 = co2
   !        6 = o2
   !        this subroutine only calculates taug below 1 mb
 
   !@Arguments
   ! taug   gaseous optical depth
   ! dp     air mass path for a model layer (explained in raddriv).
   ! o3     o3 mass mixing ratio
   ! qq     water vapor mass mixing ratio
   ! dip    interpolation factor for pressure between two
   !        neighboring standard input data pressure levels
   ! dt     layer temperature - 250 k
   ! inpt   number of the level for the standard input data pressures
   ! rmu3   a factor of solar zenith angle, given in raddriv
   !----------------------------------------------------------------------

#include "ccc1_bands.cdk"
#include "ccc_tracegases.cdk"

   if (ib .eq. 1) then

      !----------------------------------------------------------------------
      !     band (14500 - 50000 cm^-1), nongray gaseous absorption of o3,
      !     h2o and o2.
      !     relative solar energy 629.0629 wm^-2.
      !     ig9 (50000-42000)  uvc                           1.57122 (wm^-2)
      !     ig8 (42000-37400)  uvc                           2.78822
      !     ig7 (37400-35700)  uvc                           2.78486
      !     uvc all included in gh part

      !     ig6 (35700-34000)  uvb                           5.54284
      !     ig4 (34000-32185)  uvb, uva                      9.71153
      !     ig3 (32185-30300)  uva  j value: 32185 cm^-1    16.03334
      !     ig2 (30300-25000)  uva                          73.07372
      !     ig1 (25000-20000)  par                         187.18623
      !     ig1 (20000-14500)  par                         330.38562
      !     par: photosynthetic active radiation

      !     the effect of h2o and o2 is added with simple method
      !----------------------------------------------------------------------

      if (ig .eq. 1) then

         do k = lev1, lay
            do i = il1, il2
               m =  inpt(i,k)
               if (inpt(i,k).ge.950) m = inpt(i,k) - 1000

               if (m .lt. 7) then
                  xx      = (cs1o21 - 0.881e-05 * rmu3(i)) * rmo2
               else
                  xx      = (0.108e-04 - 0.881e-05 * rmu3(i)) * rmo2
               endif

               if (m .lt. 15) then
                  xx      =  xx+ (0.199e-02 - 0.952e-03 * rmu3(i)) * qq(i,k)
               else
                  xx      =  xx+ (0.208e-02 - 0.952e-03 * rmu3(i)) * qq(i,k)
               endif

               taug(i,k) = (cs1o3(ig) * o3(i,k) + xx) * dp(i,k)
            enddo
         enddo
      else
         do k = lev1, lay
            do i = il1, il2
               taug(i,k) =  cs1o3(ig) * o3(i,k) * dp(i,k)
            enddo
         enddo
      endif

      gw =  gws1(ig)

   else if (ib .eq. 2) then

      !----------------------------------------------------------------------
      !     band (8400 - 14500 cm^-1), nongray gaseous absorption of h2o,
      !     o2 and o3
      !     relative solar energy 429.5020 w m^-2
      !----------------------------------------------------------------------

      if (ig .le. 2) then
         ng2 =  6
         call ccc1_tline2 (taug, cs2h2o(1,1,ig), cs2o2(1,1,ig), qq, o3, &
              ng2, dp, dip, dt, inpt, &
              lev1, gh, mtl, il1, il2, ilg, lay, urbf)
      else
         ng =  1
         init=2
         call ccc1_tline1 (taug, cs2h2o(1,1,ig), qq, ng, dp, dip, &
              dt, inpt, lev1, gh, mtl, init, &
              il1, il2, ilg, lay, urbf)
      endif

      !----------------------------------------------------------------------
      !     simply add o3 effect
      !----------------------------------------------------------------------

      do k = lev1, lay
         do i = il1, il2
            taug(i,k)   =  taug(i,k) + cs2o3(ig) * o3(i,k) * dp(i,k)
         enddo
      enddo

      gw =  gws2(ig)

   else if (ib .eq. 3) then

      !----------------------------------------------------------------------
      !     band (4200 - 8400 cm^-1), nongray gaseous absorption of h2o
      !     and co2
      !     relative solar energy  242.8750 wm^-2
      !----------------------------------------------------------------------

      ng2 =  3
      call ccc1_tline2 (taug, cs3h2o(1,1,ig), cs3co2(1,1,ig), qq, o3, &
           ng2, dp, dip, dt, inpt, &
           lev1, gh, mtl, il1, il2, ilg, lay, urbf)

      gw =  gws3(ig)

   else if (ib .eq. 4) then

      !----------------------------------------------------------------------
      !     band (2200 - 4200 cm^-1), nongray gaseous absorption of h2o
      !     and co2
      !     relative solar energy 30.9570 w m^-2
      !----------------------------------------------------------------------

      ng2 =  3
      call ccc1_tline2 (taug, cs4h2o(1,1,ig), cs4co2(1,1,ig), qq, o3, &
           ng2, dp, dip, dt, inpt, &
           lev1, gh, mtl, il1, il2, ilg, lay, urbf)

      gw =  gws4(ig)

   endif

   return
end subroutine ccc1_gasopts
