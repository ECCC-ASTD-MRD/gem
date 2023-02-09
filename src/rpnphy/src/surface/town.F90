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

!/@*
subroutine town2(bus, bussiz, ptsurf, ptsurfsiz, dt, kount, n, m, nk)
   use iso_c_binding
   use, intrinsic :: iso_fortran_env, only: INT64
   use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
   use sfclayer,   only : sl_prelim,sl_sfclayer,SL_OK
   use modd_town,      only : xtown,                            &
        xq_town,                               &
        xu_canyon,                             &
        xrn_roof,xh_roof,xle_roof,xles_roof,   &
        xgflux_roof,xrunoff_roof,              &
        xrn_road,xh_road,xle_road,xles_road,   &
        xgflux_road,xrunoff_road,              &
        xrn_wall,xh_wall,xle_wall,xgflux_wall, &
        xrnsnow_roof,xhsnow_roof,xlesnow_roof, &
        xgsnow_roof,xmelt_roof,                &
        xrnsnow_road,xhsnow_road,xlesnow_road, &
        xgsnow_road,xmelt_road,                &
        xrn,xh,xle,xgflux,xevap,xrunoff,       &
        xch,xri,xustar,                        &
        XTRAD_IN,XTRAD_SUN,XTRAD_SHADE,        &
        XTRAD_RFSUN,XTRAD_RFSHADE,             &
        XTGLOBE_SUN,XTGLOBE_SHADE,             &
        XTGLOBE_RFSUN,XTGLOBE_RFSHADE,         &
        XTWETB,XTWETB_ROOF,                    &
        XUTCI_IN,XUTCI_OUTSUN,XUTCI_OUTSHADE,  &
        XUTCI_RFSUN,XUTCI_RFSHADE,             &
        XWBGT_OUTSUN,XWBGT_OUTSHADE,           &
        XWBGT_RFSUN,XWBGT_RFSHADE,             &
        XUTCIC_IN,XUTCIC_OUTSUN,               &
        XUTCIC_OUTSHADE,XUTCIC_RFSUN,          &
        XUTCIC_RFSHADE,                        &
        XTRFZT,XTRDZT,XURDZU,                  &
        XQ1,XQ2,XQ3,XQ4,XQ5,XQ6,XQ7,           &
        XQ8,XQ9,XQ10,XQ11,XQ12,XQ13

   use modd_teb,       only : xzs, xbld, xbld_height, xz0_town,      &
        xz0_roof,xz0_road,                     &
        xwall_o_hor, xcan_hw_ratio,            &
        xsvf_road,xsvf_wall,                   &
        xalb_roof, xalb_road, xalb_wall,       &
        xemis_roof, xemis_road, xemis_wall,    &
        xhc_roof, xtc_roof, xd_roof,           &
        xhc_road, xtc_road, xd_road,           &
        xhc_wall, xtc_wall, xd_wall,           &
        nroof_layer, nroad_layer, nwall_layer, &
        xh_traffic, xle_traffic,               &
        xh_industry, xle_industry,             &
        xt_roof, xt_road, xt_wall,             &
        xws_roof, xws_road,                    &
        xt_canyon, xq_canyon,                  &
        xti_road, xti_bld,                     &
        xwsnow_roof,  xtsnow_roof, xrsnow_roof,&
        xasnow_roof, xesnow_roof, xtssnow_roof,&
        xwsnow_road,  xtsnow_road, xrsnow_road,&
        xasnow_road, xesnow_road, xtssnow_road
   use modd_csts
   use modi_coupling_teb2, only: coupling_teb2
   use modi_sunpos
   use sfc_options, only: atm_tplus, atm_external, jdateo, zu, zt, impflx &
        ,urb_diagwind, urb_diagtemp, sl_Lmin_town, vamin
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   !@Object Choose the surface scheme for towns
   !@Arguments
   !               - Input/Output -
   ! bus           bus of surface variables
   !               - input -
   ! bussiz        size of the surface bus
   ! ptsurf        surface pointers
   ! ptsurfsiz     dimension of ptsurf
   ! kount         number of timestep
   ! dt            timestep
   ! n             running length
   ! m             horizontal dimension
   ! nk            vertical dimension

   integer bussiz, kount
   real dt
   real, target :: bus(bussiz)
   integer ptsurfsiz
   integer ptsurf(ptsurfsiz)

   !@Author Aude Lemonsu (April 2004)
   !*@/
   ! Last version Sylvie Leroyer (2018)

   include "tebcst.cdk"
! surface pointer and bus : definition
   include "town_ptr.cdk"

   integer surflen
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1

   real, external :: juliand
   integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200 !#TODO: use value from mu_jdate_mod

   integer, parameter :: indx_sfc = indx_urb
   real,    parameter :: xundef   = 999.

   integer :: n, m, nk, i, hh, mn, ss

   real    :: julien

   integer,dimension (63) :: alloc_status
   integer                :: kyear       ! current year
   integer                :: kmonth      ! current month
   integer                :: kday        ! current day
   real                   :: ptime       ! current time
   real                   :: ptstep      ! timestep

   real,          dimension(1)     :: psw_bands   ! middle wavelength of each band
   real,          dimension(n)     :: prhoa       ! air density at forcing level
   real,          dimension(n)     :: psnow       ! snow rate
   real,          dimension(n)     :: prain       ! rain rate
   real,          dimension(n)     :: psftq       ! flux of water vapor
   real,          dimension(n)     :: psfu        ! zonal momentum flux
   real,          dimension(n)     :: psfv        ! meridional momentum flux
   real,          dimension(n)     :: pemis       ! emissivity
   real,          dimension(n)     :: ptrad       ! radiative temperature
   real,          dimension(n,1)   :: pdir_alb    ! direct albedo for each band
   real,          dimension(n,1)   :: psca_alb    ! diffuse albedo for each band
   real,          dimension(n,1)   :: zdir_sw    ! direct sw for each band
   real,          dimension(n,1)   :: zsca_sw    ! diffuse sw for each band
   real,          dimension(n)     :: zvdir       ! direction of the wind
   real,          dimension(n)     :: zvmod       ! module of the wind
   real,          dimension(n)     :: lat,lon     ! latitude and longitude
   real,          dimension(n)     :: zuzu
   real, target, dimension(n) ::  zday, zheure, zmin
!!$LOGICAL diagwind_interp
!!$diagwind_interp=.FALSE.
!---------------------------------------------------------------------------

   !# in  offline mode the t-step 0 is (correctly) not performed
   if (atm_external .and. kount == 0) return

   surflen = m

! surface pointer and bus : link to bus
#include "town_ptr_as.cdk"

   nroof_layer = roof_layer
   nroad_layer = road_layer
   nwall_layer = wall_layer


! TEB specific pointer and bus : allocate
!     6. Diagnostic variables :
   allocate( xq_town      (n)             , stat=alloc_status(1) )
   allocate( xles_roof    (n)             , stat=alloc_status(2) )
   allocate( xrunoff_roof (n)             , stat=alloc_status(3) )
   allocate( xles_road    (n)             , stat=alloc_status(4) )
   allocate( xrunoff_road (n)             , stat=alloc_status(5) )
   allocate( xrnsnow_roof (n)             , stat=alloc_status(6) )
   allocate( xhsnow_roof  (n)             , stat=alloc_status(7) )
   allocate( xlesnow_roof (n)             , stat=alloc_status(8) )
   allocate( xgsnow_roof  (n)             , stat=alloc_status(9) )
   allocate( xmelt_roof   (n)             , stat=alloc_status(10) )
   allocate( xrnsnow_road (n)             , stat=alloc_status(11) )
   allocate( xhsnow_road  (n)             , stat=alloc_status(12) )
   allocate( xlesnow_road (n)             , stat=alloc_status(13) )
   allocate( xgsnow_road  (n)             , stat=alloc_status(14) )
   allocate( xmelt_road   (n)             , stat=alloc_status(15) )
   allocate( xevap        (n)             , stat=alloc_status(16) )
   allocate( xrunoff      (n)             , stat=alloc_status(17) )
   allocate( xch          (n)             , stat=alloc_status(18) )
   allocate( xri          (n)             , stat=alloc_status(19) )
   allocate( xustar       (n)             , stat=alloc_status(20) )
   allocate( xz0_roof     (n)             , stat=alloc_status(21) )
   allocate( xz0_road     (n)             , stat=alloc_status(22) )
!     7.heat-stress variables
   allocate(XTRAD_IN          (N)  ,stat=alloc_status(23)    )
   allocate(XTRAD_SUN         (N)  ,stat=alloc_status(24)    )
   allocate(XTRAD_SHADE       (N) ,stat=alloc_status(25)     )
   allocate(XTRAD_RFSUN       (N)  ,stat=alloc_status(26)    )
   allocate(XTRAD_RFSHADE     (N)   ,stat=alloc_status(27)   )
   allocate(XUTCI_IN          (N)   ,stat=alloc_status(28)   )
   allocate(XUTCI_OUTSUN      (N)  ,stat=alloc_status(29)    )
   allocate(XUTCI_OUTSHADE    (N)   ,stat=alloc_status(30)   )
   allocate(XUTCI_RFSUN       (N)   ,stat=alloc_status(31)   )
   allocate(XUTCI_RFSHADE     (N)  ,stat=alloc_status(32)   )
   allocate(XWBGT_OUTSUN      (N)  ,stat=alloc_status(33)    )
   allocate(XWBGT_OUTSHADE    (N)   ,stat=alloc_status(34)   )
   allocate(XWBGT_RFSUN       (N)   ,stat=alloc_status(35)   )
   allocate(XWBGT_RFSHADE     (N)  ,stat=alloc_status(36)   )
   allocate(XUTCIC_IN         (N)   ,stat=alloc_status(37)  )
   allocate(XUTCIC_OUTSUN     (N)   ,stat=alloc_status(38)  )
   allocate(XUTCIC_OUTSHADE   (N)   ,stat=alloc_status(39)  )
   allocate(XUTCIC_RFSUN      (N)  ,stat=alloc_status(40)   )
   allocate(XUTCIC_RFSHADE    (N)  ,stat=alloc_status(41)   )
   allocate(XTGLOBE_SUN        (N)  ,stat=alloc_status(42) )
   allocate(XTGLOBE_SHADE      (N)  ,stat=alloc_status(43) )
   allocate(XTGLOBE_RFSUN      (N)  ,stat=alloc_status(44) )
   allocate(XTGLOBE_RFSHADE    (N)  ,stat=alloc_status(45) )
   allocate(XTWETB             (N)  ,stat=alloc_status(46) )
   allocate(XTWETB_ROOF        (N)  ,stat=alloc_status(47) )
   allocate(XTRFZT     (N)  ,stat=alloc_status(48) )
   allocate(XTRDZT     (N)  ,stat=alloc_status(49) )
   allocate(XURDZU     (N)  ,stat=alloc_status(50) )
   allocate(XQ1     (N)  ,stat=alloc_status(51) )
   allocate(XQ2     (N)  ,stat=alloc_status(52) )
   allocate(XQ3     (N)  ,stat=alloc_status(53) )
   allocate(XQ4     (N)  ,stat=alloc_status(54) )
   allocate(XQ5     (N)  ,stat=alloc_status(55) )
   allocate(XQ6     (N)  ,stat=alloc_status(56) )
   allocate(XQ7     (N)  ,stat=alloc_status(57) )
   allocate(XQ8     (N)  ,stat=alloc_status(58) )
   allocate(XQ9     (N)  ,stat=alloc_status(59) )
   allocate(XQ10    (N)  ,stat=alloc_status(60) )
   allocate(XQ11    (N)  ,stat=alloc_status(61) )
   allocate(XQ12    (N)  ,stat=alloc_status(62) )
   allocate(XQ13    (N)  ,stat=alloc_status(63) )

! TEB specific pointer and bus : Initialisation
!     6. Diagnostic variables :
   xq_town         = xundef
   xles_roof       = xundef
   xrunoff_roof    = xundef
   xles_road       = xundef
   xrunoff_road    = xundef
   xrnsnow_roof    = xundef
   xhsnow_roof     = xundef
   xlesnow_roof    = xundef
   xgsnow_roof     = xundef
   xmelt_roof      = xundef
   xrnsnow_road    = xundef
   xhsnow_road     = xundef
   xlesnow_road    = xundef
   xgsnow_road     = xundef
   xmelt_road      = xundef
   xevap           = xundef
   xrunoff         = xundef
   xch             = xundef
   xri             = xundef
   xustar          = xundef
!     7. heat-stress variables
    XTRAD_IN   = XUNDEF
    XTRAD_SUN   = XUNDEF
    XTRAD_SHADE   = XUNDEF
    XTRAD_RFSUN   = XUNDEF
    XTRAD_RFSHADE   = XUNDEF
    XUTCI_IN   = XUNDEF
    XUTCI_OUTSUN   = XUNDEF
    XUTCI_OUTSHADE   = XUNDEF
    XUTCI_RFSUN   = XUNDEF
    XUTCI_RFSHADE   = XUNDEF
    XWBGT_OUTSUN   = XUNDEF
    XWBGT_OUTSHADE   = XUNDEF
    XWBGT_RFSUN   = XUNDEF
    XWBGT_RFSHADE   = XUNDEF
    XUTCIC_IN   = XUNDEF
    XUTCIC_OUTSUN   = XUNDEF
    XUTCIC_OUTSHADE   = XUNDEF
      XUTCIC_RFSUN  = XUNDEF
      XUTCIC_RFSHADE= XUNDEF
      XTGLOBE_SUN   = XUNDEF
      XTGLOBE_SHADE   = XUNDEF
      XTGLOBE_RFSUN   = XUNDEF
      XTGLOBE_RFSHADE   = XUNDEF
      XTWETB         = XUNDEF
      XTWETB_ROOF    = XUNDEF
      XQ1 = XUNDEF
      XQ2 = XUNDEF
      XQ3 = XUNDEF
      XQ4 = XUNDEF
      XQ5 = XUNDEF
      XQ6 = XUNDEF
      XQ7 = XUNDEF
      XQ8 = XUNDEF
      XQ9 = XUNDEF
      XQ10 = XUNDEF
      XQ11 = XUNDEF
      XQ12 = XUNDEF
      XQ13 = XUNDEF

   call ini_csts

!  Time
   julien       = real(jdate_day_of_year(jdateo + kount*int(dt) + MU_JDATE_HALFDAY))
   call mu_js2ymdhms(jdateo, kyear, kmonth, kday, hh, mn, ss)
   ptime        = hh*3600. + mn*60 + ss + dt*(kount)
   kday         = kday + int(ptime/86400.)

   ptime        = amod(ptime,3600*24.)
   ptstep       = dt
   psw_bands    = 0.

!   Zenithal angle computation
!   zdlat in rad
     do i=1,n
        lat(i) = zdlat(i)*180./XPI
        lon(i) = zdlon(i)*180./XPI
      end do
      zday(:)   = julien
      zheure(:) = int(ptime/3600.)*1.
      zmin(:)   = (ptime/3600.-int(ptime/3600.))*60.

      call sunpos(kyear,kmonth,kday,ptime,lon,lat,ztsun,zzenith,zazim)

! convert snow and rain rates
      do i=1,n
        psnow         (i) = zsnowrate(i) *1000.
        prain         (i) = zrainrate(i) *1000.
      enddo

!       General variables
!       -----------------
! TEB specific pointer and bus : link to bus
!     1. Geometric parameters :
      xbld          (1:n) => bus(x(bld         ,1,1) : )
      xbld_height   (1:n) => bus(x(bld_height  ,1,1) : )
      xz0_town      (1:n) => bus(x(z0_town     ,1,1) : )
      xz0_roof      (1:n) => bus(x(z0_roof     ,1,1) : )
      xz0_road      (1:n) => bus(x(z0_road     ,1,1) : )
      xwall_o_hor   (1:n) => bus(x(wall_o_hor  ,1,1) : )
      xcan_hw_ratio (1:n) => bus(x(can_hw_ratio,1,1) : )
      xsvf_road     (1:n) => bus(x(svf_road    ,1,1) : )
      xsvf_wall     (1:n) => bus(x(svf_wall    ,1,1) : )
!     2. Radiative properties :
      xalb_roof     (1:n) => bus(x(alb_roof    ,1,1) : )
      xalb_road     (1:n) => bus(x(alb_road    ,1,1) : )
      xalb_wall     (1:n) => bus(x(alb_wall    ,1,1) : )
      xemis_roof    (1:n) => bus(x(emis_roof   ,1,1) : )
      xemis_road    (1:n) => bus(x(emis_road   ,1,1) : )
      xemis_wall    (1:n) => bus(x(emis_wall   ,1,1) : )
!     4. Anthropogenic fluxes :
      xh_traffic    (1:n) => bus(x(h_traffic   ,1,1) : )
      xle_traffic   (1:n) => bus(x(le_traffic  ,1,1) : )
      xh_industry  (1:n) => bus(x(h_industry  ,1,1) : )
      xle_industry  (1:n) => bus(x(le_industry ,1,1) : )
!     4. Pronostic variables :
      xws_roof      (1:n) => bus(x(ws_roof     ,1,1) : )
      xws_road      (1:n) => bus(x(ws_road     ,1,1) : )
      xt_canyon     (1:n) => bus(x(t_canyon    ,1,1) : )
      xq_canyon     (1:n) => bus(x(q_canyon    ,1,1) : )
      xti_bld       (1:n) => bus(x(ti_bld      ,1,1) : )
      xti_road      (1:n) => bus(x(ti_road     ,1,1) : )
      xhc_roof (1:n,1:nroof_layer) => bus(x(hc_roof ,1,1) : )
      xtc_roof (1:n,1:nroof_layer) => bus(x(tc_roof ,1,1) : )
      xd_roof  (1:n,1:nroof_layer) => bus(x(d_roof  ,1,1) : )
      xt_roof  (1:n,1:nroof_layer) => bus(x(t_roof  ,1,1) : )
      xhc_road (1:n,1:nroad_layer) => bus(x(hc_road ,1,1) : )
      xtc_road (1:n,1:nroad_layer) => bus(x(tc_road ,1,1) : )
      xd_road  (1:n,1:nroad_layer) => bus(x(d_road  ,1,1) : )
      xt_road  (1:n,1:nroad_layer) => bus(x(t_road  ,1,1) : )
      xhc_wall (1:n,1:nwall_layer) => bus(x(hc_wall ,1,1) : )
      xtc_wall (1:n,1:nwall_layer) => bus(x(tc_wall ,1,1) : )
      xd_wall  (1:n,1:nwall_layer) => bus(x(d_wall  ,1,1) : )
      xt_wall  (1:n,1:nwall_layer) => bus(x(t_wall  ,1,1) : )
!     5. Snow variables :
      xwsnow_roof (1:n) => bus(x(sroof_wsnow,1,1) : )
      xtsnow_roof (1:n) => bus(x(sroof_t    ,1,1) : )
      xrsnow_roof (1:n) => bus(x(sroof_rho  ,1,1) : )
      xasnow_roof (1:n) => bus(x(sroof_alb  ,1,1) : )
      xesnow_roof (1:n) => bus(x(sroof_emis ,1,1) : )
      xtssnow_roof(1:n) => bus(x(sroof_ts   ,1,1) : )
      xwsnow_road (1:n) => bus(x(sroad_wsnow,1,1) : )
      xtsnow_road (1:n) => bus(x(sroad_t    ,1,1) : )
      xrsnow_road (1:n) => bus(x(sroad_rho  ,1,1) : )
      xasnow_road (1:n) => bus(x(sroad_alb  ,1,1) : )
      xesnow_road (1:n) => bus(x(sroad_emis ,1,1) : )
      xtssnow_road(1:n) => bus(x(sroad_ts   ,1,1) : )
!     6. Diagnostic variables :
      xu_canyon   (1:n) => bus(x( u_canyon,1,1)  : )
      xrn         (1:n) => bus(x( rn_town ,1,1)  : )
      xh          (1:n) => bus(x( h_town  ,1,1)  : )
      xle         (1:n) => bus(x( le_town ,1,1)  : )
      xgflux      (1:n) => bus(x( g_town  ,1,1)  : )
      xrn_roof    (1:n) => bus(x( rn_roof ,1,1)  : )
      xh_roof     (1:n) => bus(x( h_roof  ,1,1)  : )
      xle_roof    (1:n) => bus(x( le_roof ,1,1)  : )
      xgflux_roof (1:n) => bus(x( g_roof  ,1,1)  : )
      xrn_road    (1:n) => bus(x( rn_road ,1,1)  : )
      xh_road     (1:n) => bus(x( h_road  ,1,1)  : )
      xle_road    (1:n) => bus(x( le_road ,1,1)  : )
      xgflux_road (1:n) => bus(x( g_road  ,1,1)  : )
      xrn_wall    (1:n) => bus(x( rn_wall ,1,1)  : )
      xh_wall     (1:n) => bus(x( h_wall  ,1,1)  : )
      xle_wall    (1:n) => bus(x( le_wall ,1,1)  : )
      xgflux_wall (1:n) => bus(x( g_wall  ,1,1)  : )
!     7. Heat-stress variables :
      xtrad_in       (1:n) => bus( x( yradin,1,1)         : )
      xtrad_sun      (1:n) => bus( x( yradsun,1,indx_urb)        : )
      xtrad_shade    (1:n) => bus( x( yradshade,1,indx_urb)      : )
      xtrad_rfsun    (1:n) => bus( x(yradrfsun ,1,1)      : )
      xtrad_rfshade  (1:n) => bus( x( yradrfshade,1,1)    : )
      xutci_in      (1:n) => bus( x( yutciin,1,1)        : )
      xutci_outsun     (1:n) => bus( x( yutcisun,1,indx_urb)       : )
      xutci_outshade   (1:n) => bus( x( yutcishade,1,indx_urb)     : )
      xutci_rfsun   (1:n) => bus( x( yutcirfsun,1,1)     : )
      xutci_rfshade (1:n) => bus( x( yutcirfshade,1,1)   : )
      xwbgt_outsun     (1:n) => bus( x( ywbgtsun,1,indx_urb)       : )
      xwbgt_outshade   (1:n) => bus( x( ywbgtshade,1,indx_urb)     : )
      xwbgt_rfsun   (1:n) => bus( x( ywbgtrfsun,1,1)     : )
      xwbgt_rfshade (1:n) => bus( x( ywbgtrfshade,1,1)   : )
      xutcic_in     (1:n) => bus( x( yutcicin,1,1)       : )
      xutcic_outsun    (1:n) => bus( x(yutcicsun ,1,1)      : )
      xutcic_outshade  (1:n) => bus( x(yutcicshade ,1,1)    : )
      xutcic_rfsun  (1:n) => bus( x(yutcicrfsun ,1,1)    : )
      xutcic_rfshade(1:n) => bus( x(yutcicrfshade ,1,1)   : )
      xtglobe_sun     (1:n) => bus( x( ytglbsun,1,indx_urb)       : )
      xtglobe_shade   (1:n) => bus( x( ytglbshade,1,indx_urb)     : )
      xtwetb    (1:n) => bus( x(ytwetb ,1,indx_urb)         : )
      xtglobe_rfsun    (1:n) => bus( x( ytglbrfsun,1,1)         : )
      xtglobe_rfshade   (1:n) => bus( x( ytglbrfshade,1,1)         : )
      xtwetb_roof   (1:n) => bus( x( ytwetbrf,1,1)         : )
      xtrfzt   (1:n) => bus( x( ytrfzt,1,1)         : )
      xtrdzt   (1:n) => bus( x( ytrdzt,1,1)         : )
      xurdzu   (1:n) => bus( x( yurdzu,1,1)         : )
      xQ1   (1:n) => bus( x( yQ1,1,indx_urb)         : )
      xQ2   (1:n) => bus( x( yQ2,1,indx_urb)         : )
      xQ3   (1:n) => bus( x( yQ3,1,indx_urb)         : )
      xQ4   (1:n) => bus( x( yQ4,1,indx_urb)         : )
      xQ5   (1:n) => bus( x( yQ5,1,indx_urb)         : )
      xQ6   (1:n) => bus( x( yQ6,1,indx_urb)         : )
      xQ7   (1:n) => bus( x( yQ7,1,indx_urb)         : )
      xQ8   (1:n) => bus( x( yQ8,1,1)         : )
      xQ9   (1:n) => bus( x( yQ9,1,1)         : )
      xQ10  (1:n) => bus( x( yQ10,1,1)        : )
      xQ11  (1:n) => bus( x( yQ11,1,1)        : )
      xQ12  (1:n) => bus( x( yQ12,1,1)        : )
      xQ13  (1:n) => bus( x( yQ13,1,1)        : )


      do i=1,n
    !# in  offline mode no forcing of diffuse (scattered) radiation
    if (atm_external) then
    !# about 13-15 % (Leroyer et al. 2018, urban climate), we should consider a diurnal cycle sunset/sunrise
        zdir_sw  (i,1) = 0.85*psol_sw(i,1)
        zsca_sw  (i,1) = 0.15*psol_sw(i,1)
    else
        zdir_sw  (i,1) = pdir_sw(i,1)
        zsca_sw  (i,1) = psca_sw(i,1)
    endif
    enddo

       ! coherence between solar zenithal angle and radiation
       !
       where (sum(zdir_sw+zsca_sw,2)>0.)
         zzenith = min (zzenith,xpi/2.-0.01)
       elsewhere
         zzenith = max (zzenith,xpi/2.)
       end where

   !     PRELIM
   !     ----
      I = SL_PRELIM(PTA,PQA,PU,PV,PPS,PUREF,MIN_WIND_SPEED=VAMIN,SPD_AIR=ZVMOD,DIR_AIR=ZVDIR, &
           RHO_AIR=PRHOA)

   IF (I /= SL_OK) THEN
      call physeterror('town', 'error returned by sl_prelim()')
      return
   endif

!-----------------------------------------------------------------------------

      call coupling_teb2(ptstep, &
              ztsun, zzenith, pzref, puref, pu, pv, pqa, pta, &
              prhoa, prain, psnow, plw, zdir_sw, zsca_sw, psw_bands, pps,  &
              ppa, psftq,  zfc, psfu, psfv,                                &
              ptrad, pdir_alb, psca_alb, pemis, zdlat)

!-----------------------------------------------------------------------------

      do i=1,n
!       Variables a agreger
!       -------------------
        zz0(i) = xz0_town   (i)
!   evaluation of the thermal roughness lenght for the diagnostic (first guess)
!   will come from the bus later when z0 is read in geophy (not this version)

        zz0t(i) =  MAX( xz0_town(i)                      *          &
                7.4 * exp( - 1.29 *(  xz0_town(i)        *          &
                0.4 * (pu(i)**2 + pv(i) **2 )**0.5                  &
                /log((puref(i)+xbld_height(i)/3.)/xz0_town(i))      &
                /1.461e-5)**0.25), 1.e-20 )
        ztsurf (i) = ptrad      (i)
        ztsrad (i) = ptrad      (i)
        zqsurf (i) = xq_town    (i)
        zalvis (i) = pdir_alb   (i,1)
        zsnodp (i) = xbld(i)      * (xwsnow_roof(i)/xrsnow_roof(i))  +   &
                    (1.-xbld(i))  * (xwsnow_road(i)/xrsnow_road(i))
        zfv    (i) = psftq      (i) * xlvtt
!  runoff -- aggregated with all other surface tiles. Convert to mm
        zrunofftot(i) = xrunoff(i) * dt
        zalscatw (i)  = psca_alb   (i,1)
        zemtw (i)     = pemis   (i)
        ztsradtw (i ) = ptrad   (i)
      end do

!-----------------------------------------------------------------------------
!   Compute town surface layer var. (for alfat alfaq ctu...)
!-----------------------------------------------------------------------------
   !# Compute town surface layer var. (need for  udiag, zbm & the implicit condition)
   i = sl_sfclayer(pthetaa,pqa,zvmod,zvdir,puref,pzref,ztsurf,zqsurf,  &
        zz0,zz0t,zdlat,zfcor,optz0=0,L_min=sl_Lmin_town,hghtm_diag=zu, &
        hghtt_diag=zt,ilmo=zilmo,h=zhst,ue=zfrv,flux_t=zftemp,flux_q=zfvap,  &
        coefm=zbm,coeft=zbt,u_diag=zudiag,v_diag=zvdiag)

   if (i /= SL_OK) then
      call physeterror('town', 'error returned by sl_sfclayer(number 1)')
      return
   endif

      do i=1,n
         zalfat(i) = -zfc(i)/(xcpd*prhoa(i))
         zalfaq(i) = -psftq(i)
        if (.not.impflx) zbt(i) = 0.
        if (impflx) then
!         zalfat(i) = zalfat(i) - zbt(i) * pthetaa(i)
!         zalfaq(i) = zalfaq(i) - zbt(i) * pqa(i)
         zalfat(i) = zalfat(i) - zbt(i) * pthetaa(i)
         zalfaq(i) = zalfaq(i) - zbt(i) * pqa(i)
        endif
      end do

!-----------------------------------------------------------------------------
!   Sreen-level Diagnostics
!-----------------------------------------------------------------------------
      ! => temperature (default diag is t_diag from sl_sfclayer number 2 and T_can)
      if( urb_diagtemp ) then 
      !  ztdiag = xtrdzt
        ztdiag = xt_canyon
      else
      ! sl_sfclayer between  road and mid-canyon for zt => compute ztdiag
     i = sl_sfclayer(xt_canyon,xq_canyon,xu_canyon,zvdir,xbld_height/2.0,xbld_height/2.0,xt_road(:,1),xq_canyon, &
        xz0_road,xz0_road/10.0,zdlat,zfcor,optz0=8,L_min=sl_Lmin_town,hghtm_diag=zu,hghtt_diag=zt,      &
        t_diag=ztdiag)

      if (i /= SL_OK) then
         call physeterror('town', 'error 2 returned by sl_sfclayer(number 2)')
         return
      endif
      endif
      ! => humidity
        zqdiag = xq_canyon

      ! special conditions
      do i=1,n
      if (xbld_height(i) .le. (zt*2.0) ) then 
        ztdiag(i) = xt_canyon(i)
      endif
      ! => wind (default diag is u_diag from sl_sfclayer number 1) (option urb_diagwind for adaptation with street morphology)
      ! if bldh>15m => u_canyon
      if (xbld_height(i) .ge. (zu*3.0/2.0) ) then
        zudiag(i) = xu_canyon(i) *cos(zvdir(i))
        zvdiag(i) = xu_canyon(i) *sin(zvdir(i))
      endif

      if( urb_diagwind ) then
      !  linear interpolation between two cases du/dz=Utop-Ucan / H-2H/3
      !        => U(z=10)= U10 = (3 zu/H-2) Utop + (-3 zu/H+3) Ucan  
      if( xbld_height(i) .gt. zu .and. xbld_height(i) .lt. (zu*3.0/2.0) )  then 
        zuzu(i) =  (3.0 * zu /xbld_height(i) -2.) * zvmod(i)                  &
           * log( (     xbld_height(i)/3.) / xz0_town (i))   &
           / log( (puref(i) + xbld_height(i)/3.) / xz0_town   (i))   &
           + (- 3.*zu /xbld_height(i) +3.) * xu_canyon(i)
        zudiag(i) = zuzu(i) *COS(zvdir(i))
        zvdiag(i) = zuzu(i) *SIN(zvdir(i))
      elseif (xbld_height(i) .le. zu ) then
      !        => log law. above roof level 
      zuzu(i) =  zvmod(i)                                               &
           * log( (     zu    - 2 * xbld_height(i)/3.) / xz0_town(i))   &
           / log( ( puref(i)  + xbld_height(i)/3.) / xz0_town(i))
        zudiag(i) = zuzu(i) *cos(zvdir(i))
        zvdiag(i) = zuzu(i) *sin(zvdir(i))
      endif
      endif
      end do

      ! Fill surface type-specific diagnostic values
      zqdiagtyp = zqdiag
      ztdiagtyp = ztdiag
      zudiagtyp = zudiag
      zvdiagtyp = zvdiag
      zqdiagtypv = zqdiag
      ztdiagtypv = ztdiag
      zudiagtypv = zudiag
      zvdiagtypv = zvdiag

      call fillagg ( bus, bussiz, ptsurf, ptsurfsiz, indx_urb, surflen )

!     6. diagnostic variables
      deallocate( xq_town           , stat=alloc_status(1) )
      deallocate( xles_roof         , stat=alloc_status(2) )
      deallocate( xrunoff_roof      , stat=alloc_status(3) )
      deallocate( xles_road         , stat=alloc_status(4) )
      deallocate( xrunoff_road      , stat=alloc_status(5) )
      deallocate( xrnsnow_roof      , stat=alloc_status(6) )
      deallocate( xhsnow_roof       , stat=alloc_status(7) )
      deallocate( xlesnow_roof      , stat=alloc_status(8) )
      deallocate( xgsnow_roof       , stat=alloc_status(9) )
      deallocate( xmelt_roof        , stat=alloc_status(10) )
      deallocate( xrnsnow_road      , stat=alloc_status(11) )
      deallocate( xhsnow_road       , stat=alloc_status(12) )
      deallocate( xlesnow_road      , stat=alloc_status(13) )
      deallocate( xgsnow_road       , stat=alloc_status(14) )
      deallocate( xmelt_road        , stat=alloc_status(15) )
      deallocate( xevap             , stat=alloc_status(16) )
      deallocate( xrunoff           , stat=alloc_status(17) )
      deallocate( xch               , stat=alloc_status(18) )
      deallocate( xri               , stat=alloc_status(19) )
      deallocate( xustar            , stat=alloc_status(20) )
      deallocate( xz0_roof          , stat=alloc_status(21) )
      deallocate( xz0_road          , stat=alloc_status(22) )
!     7. heat-stress indices
      deallocate( XTRAD_IN        ,stat=alloc_status(23)  )
      deallocate( XTRAD_SUN      ,stat=alloc_status(24)   )
      deallocate( XTRAD_SHADE     ,stat=alloc_status(25)  )
      deallocate( XTRAD_RFSUN      ,stat=alloc_status(26) )
      deallocate( XTRAD_RFSHADE    ,stat=alloc_status(27) )
      deallocate( XUTCI_IN        ,stat=alloc_status(28)  )
      deallocate( XUTCI_OUTSUN     ,stat=alloc_status(29) )
      deallocate( XUTCI_OUTSHADE   ,stat=alloc_status(30) )
      deallocate( XUTCI_RFSUN     ,stat=alloc_status(31)  )
      deallocate( XUTCI_RFSHADE    ,stat=alloc_status(32) )
      deallocate( XWBGT_OUTSUN     ,stat=alloc_status(33) )
      deallocate( XWBGT_OUTSHADE   ,stat=alloc_status(34) )
      deallocate( XWBGT_RFSUN     ,stat=alloc_status(35)  )
      deallocate( XWBGT_RFSHADE    ,stat=alloc_status(36) )
      deallocate( XUTCIC_IN        ,stat=alloc_status(37))
      deallocate( XUTCIC_OUTSUN    ,stat=alloc_status(38) )
      deallocate( XUTCIC_OUTSHADE  ,stat=alloc_status(39) )
      deallocate( XUTCIC_RFSUN     ,stat=alloc_status(40) )
      deallocate( XUTCIC_RFSHADE   ,stat=alloc_status(41) )
      deallocate( XTGLOBE_SUN      ,stat=alloc_status(42) )
      deallocate( XTGLOBE_SHADE    ,stat=alloc_status(43) )
      deallocate( XTGLOBE_RFSUN    ,stat=alloc_status(44) )
      deallocate( XTGLOBE_RFSHADE  ,stat=alloc_status(45) )
      deallocate( XTWETB           ,stat=alloc_status(46) )
      deallocate( XTWETB_ROOF      ,stat=alloc_status(47) )
      deallocate( XTRFZT          ,stat=alloc_status(48) )
      deallocate( XTRDZT          ,stat=alloc_status(49) )
      deallocate( XURDZU          ,stat=alloc_status(50) )
      deallocate( XQ1          ,stat=alloc_status(51) )
      deallocate( XQ2          ,stat=alloc_status(52) )
      deallocate( XQ3          ,stat=alloc_status(53) )
      deallocate( XQ4          ,stat=alloc_status(54) )
      deallocate( XQ5          ,stat=alloc_status(55) )
      deallocate( XQ6          ,stat=alloc_status(56) )
      deallocate( XQ7          ,stat=alloc_status(57) )
      deallocate( XQ8          ,stat=alloc_status(58) )
      deallocate( XQ9          ,stat=alloc_status(59) )
      deallocate( XQ10         ,stat=alloc_status(60) )
      deallocate( XQ11         ,stat=alloc_status(61) )
      deallocate( XQ12         ,stat=alloc_status(62) )
      deallocate( XQ13         ,stat=alloc_status(63) )

      !--------------------------------------------------------------------
   return
end subroutine town2
