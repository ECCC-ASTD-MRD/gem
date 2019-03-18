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

subroutine initown2(ni, nk, trnch)
   use sfc_options
   use sfcbus_mod
   implicit none
#include <arch_specific.hf>

   integer ni, nk, trnch

   !@Author Aude Lemonsu (March 2004)
   !@Object Transfer and initialize initial fields for TEB
   !@Arguments
   !       - Input -
   ! ni       horizontal dimension
   ! nk       vertical dimension
   !

   real, parameter :: undef = 999
   integer i, k

   include "initown_ptr.cdk"
#include "initown_ptr_as.cdk"

   ! Initialize  urban parameters and prognostic var.

   DO_I: do i=1,ni
      if(zblden(i).le.0.0.or.zpaven(i).le.0.0)then
         zurban(i) = 0.0
      else
         zurban (i) = min(zblden(i)+zpaven(i), 1.0)
      endif
      if(zurban(i).gt.0.0) then
         zbld  (i) = zblden (i) / zurban (i)
         zpav  (i) = zpaven (i) / zurban (i)
      else
         zbld  (i) = 0.0
         zpav  (i) = 0.0
      endif

! min and max
      zwall_o_hor  (i) = max(zwall_o_horen(i), 0.009)
      zbld_height  (i) = min(max(zbld_heighten(i), 3.0), 30.0)
      zz0_town     (i) = min(max(zz0_townen   (i), 0.1), 15.0)

      zz0_road     (i) = zz0_roaden   (i)
      zz0_roof     (i) = zz0_roofen   (i)

      zemis_road   (i) = zemis_roaden (i)
      zemis_roof   (i) = zemis_roofen (i)
      zemis_wall   (i) = zemis_wallen (i)
      zalb_road    (i) = zalb_roaden  (i)
      zalb_roof    (i) = zalb_roofen  (i)
      zalb_wall    (i) = zalb_wallen  (i)

      if (zurban(i) > 0.0) then
         ! calculate canyon aspect ration and sky view factors.
         zcan_hw_ratio (i) = 0.5 * zwall_o_hor(i) / (1 - zbld(i))
         zsvf_road     (i) = sqrt(zcan_hw_ratio(i)**2+1.) -zcan_hw_ratio(i)
         zsvf_wall     (i) =  0.5 *  (zcan_hw_ratio(i) + 1.&
              - sqrt(zcan_hw_ratio(i)**2+1.))/zcan_hw_ratio(i)
      endif

      do k=1,3
         zhc_road (i,k) = zhc_roaden (i,k)
         ztc_road (i,k) = ztc_roaden (i,k)
         zd_road  (i,k) = zd_roaden  (i,k)

         ztc_roof (i,k) = ztc_roofen (i,k)
         zhc_roof (i,k) = zhc_roofen (i,k)
         zd_roof  (i,k) = zd_roofen  (i,k)


         zhc_wall (i,k) = zhc_wallen (i,k)
         ztc_wall (i,k) = ztc_wallen (i,k)
         zd_wall  (i,k) = zd_wallen  (i,k)
      end do
      zhc_roof (i,1)    = max(zhc_roofen(i,1), 1.e6)

      zh_industry  (i) = zbld (i)*  0.0
      zh_traffic   (i) = zpav (i)*  0.0
      zle_industry (i) = zbld (i)*  0.0
      zle_traffic  (i) = zpav (i)*  0.0

      zws_roof (i) = zws_roofen (i)
      zws_road (i) = zws_roaden (i)
      zti_road (i) = zti_roaden (i)
      zti_bld  (i) = zti_blden  (i)


      do k=1,3
         zt_roof (i,k) = zt_roofen (i,k)
         zt_road (i,k) = zt_roaden (i,k)
         zt_wall (i,k) = zt_wallen (i,k)
      end do

      zt_canyon (i) =  zt_canyonen (i)
      zq_canyon (i) =  zq_canyonen (i)

      zsroof_wsnow  (i) =  zsroof_wsnowen (i)
      if(zsroof_wsnow(i) .lt. critsnow ) then
         zsroof_wsnow  (i) = 0.0
         zsroof_t      (i) = undef
         zsroof_rho    (i) = undef
         zsroof_alb    (i) = undef
         zsroof_emis   (i) = undef
         zsroof_ts     (i) = undef
      else
         zsroof_t      (i) =  zsroof_ten     (i)
         zsroof_rho    (i) =  zsroof_rhoen   (i)
         zsroof_alb    (i) =  zsroof_alben   (i)
         zsroof_emis   (i) =  zsroof_emisen  (i)
         zsroof_ts     (i) =  zsroof_tsen    (i)
      endif

      zsroad_wsnow  (i) =  zsroad_wsnowen (i)
      if( zsroad_wsnow (i) .lt. critsnow ) then
         zsroad_wsnow  (i) =  0.0
         zsroad_t      (i) = undef
         zsroad_rho    (i) = undef
         zsroad_alb    (i) = undef
         zsroad_emis   (i) = undef
         zsroad_ts     (i) = undef
      else
         zsroad_t      (i) =  zsroad_ten     (i)
         zsroad_rho    (i) =  zsroad_rhoen   (i)
         zsroad_alb    (i) =  zsroad_alben   (i)
         zsroad_emis   (i) =  zsroad_emisen  (i)
         zsroad_ts     (i) =  zsroad_tsen    (i)
         zsnodp        (i) =  0.0
      endif

   enddo DO_I


   ! Mask (some) urban variables outside of urban areas with 999.
   !TODO: use of where could be more efficient then if inside a loop
   DO_I2: do i=1,ni

      if (zurban(i) <= 0.) then
         do k=1,3
            zt_roof(i,k) = undef
            zt_road(i,k) = undef
            zt_wall(i,k) = undef
         end do
         zt_canyon     (i) = undef
         zq_canyon     (i) = undef
         zsvf_road     (i) = undef
         zsvf_wall     (i) = undef
         zws_roof      (i) = undef
         zws_road      (i) = undef
         zti_bld       (i) = undef
         zti_road      (i) = undef
         zt_canyon     (i) = undef
         zq_canyon     (i) = undef
         zsroof_wsnow  (i) = undef
         zsroof_t      (i) = undef
         zsroof_rho    (i) = undef
         zsroof_alb    (i) = undef
         zsroof_emis   (i) = undef
         zsroof_ts     (i) = undef
         zsroad_wsnow  (i) = undef
         zsroad_t      (i) = undef
         zsroad_rho    (i) = undef
         zsroad_alb    (i) = undef
         zsroad_emis   (i) = undef
         zsroad_ts     (i) = undef
      endif

   end do DO_I2

   return
end subroutine initown2
