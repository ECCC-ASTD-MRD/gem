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
subroutine water2(bus, bussiz, ptsurf, ptsurfsiz, lcl_indx, kount, &
     n, m, nk)
   use tdpack
   use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
   use cpl_itf     , only: cpl_update
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   !@Object Calculate:   - surface roughness length (Z0) over open water
   !                       (not covered by ice) using Charnock's relation with a
   !                       BETA1 parameter (optional sea state dependency);
   !                     - surface fluxes of momentum, heat and moisture.
   !@Arguments
   !             - Input/Output -
   ! BUS         Bus for the WATER surface scheme
   !             - Input -
   ! BUSSIZ      dimension of bus
   ! PTSURF      surface pointers
   ! PTSURFSIZ   dimension of ptsurf
   ! KOUNT       timestep number
   ! N           horizontal dimension (row length)
   ! M           horizontal dimensions of fields
   !             (not used for the moment)
   ! NK          vertical dimension

   integer :: bussiz, kount, N, M, NK
   real, target :: bus(bussiz)
   integer :: ptsurfsiz
   integer :: ptsurf(ptsurfsiz), lcl_indx(2,n)

   !@Author J. Mailhot, S. Belair and B. Bilodeau (Dec 1998)
   !Revisions
   ! 001      M. Carrera and V. Fortin (Nov 2007) - Compute total runoff
   !                as precipitation - evaporation (no storage)
   !@Notes
   !          Z0 = BETA1*USTAR**2/GRAV (in metres) with minimum value
   !          Z0MIN and a maximum value Z0MAX
   !*@/
   include "sfcinput.cdk"

   integer :: surflen
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1

   real, parameter :: z0max = 5.e-3
   real, parameter :: LAKE_MASK_THRESHOLD = 0.001
   real, parameter :: am = 0.11        !Roughness adjustment coefficient (Beljaars qjrms 1994)
   real, parameter :: beta1 = 0.018    !Charnock's coefficient
   real, parameter :: k_visc = 1.5e-5  !Kinematic viscosity estimate at ~20oc
   real, parameter :: zfrvmin = 0.01   !Minimum friction velocity for inverse relations
   real, parameter :: rho_w = 1025.    !Density of seawater (Rudnick and Martin; Dyn. of Atm. and Oceans, 2002)
   real, parameter :: diff_w = 1.45e-7 !Thermal moelcular diffusivity of seawater (National Physical Laboratory tables)
   real, parameter :: cond_w = 0.58    !Thermal conducitivity of seawater (National Physical Laboratory tables)
   real, parameter :: cv_w = 3950.     !Volumetric heat capacity of seawater (National Physical Laboratory tables)
   real, parameter :: nu = 0.3         !Warm layer profile shape parameter (1 = linear)
   real, parameter :: base_depth = 3.  !Depth (m) of the "base" temperature without diurnal variability
   real, parameter :: skin_depth_min = 0.0002 !Minimum depth of the cold skin (Tu and Tsuang; GRL 2005)
   real, parameter :: ref_depth = 10.  !Reference depth (m) for the Artale et al. (JGR 2002) forumulation of Saunders' constant
   real, parameter :: emis = 0.984     !Sea surface emissivity (Konda et al.; J. of Oceanography 1994)
   character(len = *), parameter :: saunders = 'artale' !Calculation strategy for Saunders' constant ('fairall','artale')
   real, parameter :: dsst_epsilon = 1.e-5 !Convergence threshold for warm layer increment
   real, parameter :: delh_tresh = 1.e-10  !Treshold to compute ctu in ocean coupled mode
   real, parameter :: zt_rho = 1.5     !Height used to compute air density in ocean coupled mode
   integer, parameter :: max_iter = 10 !Maximum number of iterations for convergence of warm layer increment
   integer, parameter ::INDX_SFC = INDX_WATER
   logical, parameter :: WATER_TDIAGLIM = .false.

   real,pointer,dimension(:) ::  alvis_wat, cmu, ctu, fc_wat, mlac
   real,pointer,dimension(:) ::  fv_wat, zemisr, zalwater
   real,pointer,dimension(:) ::  hst_wat, hu, ilmo_wat
   real,pointer,dimension(:) ::  ps, qs, th, ts, tt, uu, vv
   real,pointer,dimension(:) ::  z0h, z0m, zdsst
   real,pointer,dimension(:) ::  zalfaq, zalfat, zdlat, zfcor
   real,pointer,dimension(:) ::  zftemp, zfvap, zqdiag, zrainrate
   real,pointer,dimension(:) ::  zrunofftot, zsnowrate, ztdiag
   real,pointer,dimension(:) ::  ztsurf, ztsrad, zudiag, zvdiag
   real,pointer,dimension(:) ::  zfrv, zzusl, zztsl
   real,pointer,dimension(:) ::  zflusolis,zfdsi,zskin_depth,zskin_inc
   real,pointer,dimension(:) :: zqdiagtyp, ztdiagtyp, zudiagtyp, zvdiagtyp
   real,pointer,dimension(:) :: zqdiagtypv, ztdiagtypv, zudiagtypv, zvdiagtypv
   real,pointer,dimension(:) ::  zfsd, zfsf, zcoszeni
   real,pointer,dimension(:) ::  zutcisun, zutcishade, zwbgtsun, zwbgtshade
   real,pointer,dimension(:) ::  zradsun, zradshade, ztglbsun, ztglbshade
   real,pointer,dimension(:) ::  ztwetb, zq1, zq2, zq3, zq4, zq5, zq6, zq7


   real, dimension(n) :: z0m_adjust,rhoa,vmod,vdir,vmodd,vmod0
   real, dimension(n) :: alpha_w,rho_a,frv_w,frv_a,visc_w,skin_solar,q_bal,skin_q,sst
   real, dimension(n) :: warm_increment,lambda,ud,vdi,gamma
   real, dimension(n) :: this_inc,prev_inc,denom
   real, dimension(n) :: my_ta,my_qa
   real, dimension(n) :: zu10,zusr          ! wind at 10m and sensor level
   real, dimension(n) :: zref_sw_surf, zemit_lw_surf, zzenith
   real, dimension(n) :: zusurfzt, zvsurfzt, zqd

   integer I
   real qsat_o_salty, delh, delq, zw, z1, z2, re
   logical            :: cplupd

   !** ------------------------------------------------------------------

   SURFLEN = M

   alvis_wat(1:n) => bus( x(alvis,1,indx_sfc) : )
   cmu      (1:n) => bus( x(bm,1,1)           : )
   ctu      (1:n) => bus( x(bt,1,indx_sfc)    : )
   mlac     (1:n) => bus( x(ml,1,1)           : )
   fc_wat   (1:n) => bus( x(fc,1,indx_sfc)    : )
   fv_wat   (1:n) => bus( x(fv,1,indx_sfc)    : )
   hst_wat  (1:n) => bus( x(hst,1,indx_sfc)   : )
   ilmo_wat (1:n) => bus( x(ilmo,1,indx_sfc)  : )
   qs       (1:n) => bus( x(qsurf,1,indx_sfc) : )
   ts       (1:n) => bus( x(twater,1,1)       : )
   z0h      (1:n) => bus( x(z0t,1,indx_sfc)   : )
   z0m      (1:n) => bus( x(z0,1,indx_sfc)    : )
   zalfaq   (1:n) => bus( x(alfaq,1,1)        : )
   zalfat   (1:n) => bus( x(alfat,1,1)        : )
   zalwater (1:n) => bus( x(alwater,1,1)      : )
   zdlat    (1:n) => bus( x(dlat,1,1)         : )
   zdsst    (1:n) => bus( x(dsst,1,1)         : )
   zemisr   (1:n) => bus( x(emisr,1,1)        : )
   zfcor    (1:n) => bus( x(fcor,1,1)         : )
   zfdsi    (1:n) => bus( x(fdsi,1,1)         : )
   zflusolis(1:n) => bus( x(flusolis,1,1)     : )
   zftemp   (1:n) => bus( x(ftemp,1,indx_sfc) : )
   zfvap    (1:n) => bus( x(fvap,1,indx_sfc)  : )
   zrainrate(1:n) => bus( x(rainrate,1,1)     : )
   zrunofftot(1:n) => bus( x(runofftot,1,indx_sfc) : )
   zsnowrate(1:n) => bus( x(snowrate,1,1)     : )
   zskin_depth(1:n) => bus( x(skin_depth,1,1) : )
   zskin_inc(1:n) => bus( x(skin_inc,1,1)     : )
   ztsurf   (1:n) => bus( x(tsurf,1,indx_sfc) : )
   ztsrad   (1:n) => bus( x(tsrad,1,1)        : )
   zudiag   (1:n) => bus( x(udiag,1,1)        : )
   zudiagtyp(1:n) => bus( x(udiagtyp,1,indx_sfc) : )
   zudiagtypv(1:n) => bus( x(udiagtypv,1,indx_sfc) : )
   zvdiag   (1:n) => bus( x(vdiag,1,1)        : )
   zvdiagtyp(1:n) => bus( x(vdiagtyp,1,indx_sfc) : )
   zvdiagtypv(1:n) => bus( x(vdiagtypv,1,indx_sfc) : )
   ztdiag   (1:n) => bus( x(tdiag,1,1)        : )
   ztdiagtyp(1:n) => bus( x(tdiagtyp,1,indx_sfc) : )
   ztdiagtypv(1:n) => bus( x(tdiagtypv,1,indx_sfc) : )
   zqdiag   (1:n) => bus( x(qdiag,1,1)        : )
   zqdiagtyp(1:n) => bus( x(qdiagtyp,1,indx_sfc) : )
   zqdiagtypv(1:n) => bus( x(qdiagtypv,1,indx_sfc) : )
   zfrv     (1:n) => bus( x(frv,1,indx_sfc)   : )
   zzusl    (1:n) => bus( x(zusl,1,1)         : )
   zztsl    (1:n) => bus( x(ztsl,1,1)         : )
   zcoszeni (1:n)   => bus( x(cang,1,1)        : )
   zfsd (1:n)       => bus( x(fsd,1,1)         : )
   zfsf (1:n)       => bus( x(fsf,1,1)         : )
   zutcisun (1:n)   => bus( x(yutcisun,1,indx_sfc)    : )
   zutcishade (1:n) => bus( x(yutcishade,1,indx_sfc)  : )
   zwbgtsun (1:n)   => bus( x(ywbgtsun,1,indx_sfc)    : )
   zwbgtshade (1:n) => bus( x(ywbgtshade,1,indx_sfc)  : )
   zradsun (1:n)    => bus( x(yradsun,1,indx_sfc)     : )
   zradshade (1:n)  => bus( x(yradshade,1,indx_sfc)   : )
   ztglbsun (1:n)   => bus( x(ytglbsun,1,indx_sfc)    : )
   ztglbshade (1:n) => bus( x(ytglbshade,1,indx_sfc)  : )
   ztwetb (1:n)     => bus( x(ytwetb,1,indx_sfc)      : )
   zq1 (1:n)        => bus( x(yq1,1,indx_sfc)         : )
   zq2 (1:n)        => bus( x(yq2,1,indx_sfc)         : )
   zq3 (1:n)        => bus( x(yq3,1,indx_sfc)         : )
   zq4 (1:n)        => bus( x(yq4,1,indx_sfc)         : )
   zq5 (1:n)        => bus( x(yq5,1,indx_sfc)         : )
   zq6 (1:n)        => bus( x(yq6,1,indx_sfc)         : )
   zq7 (1:n)        => bus( x(yq7,1,indx_sfc)         : )

   if (atm_tplus) then
      hu       (1:n) => bus( x(huplus,1,nk)      : )
      ps       (1:n) => bus( x(pplus,1,1)        : )
      th       (1:n) => bus( x(thetaap,1,1)      : )
      tt       (1:n) => bus( x(tplus,1,nk)       : )
      uu       (1:n) => bus( x(uplus,1,nk)       : )
      vv       (1:n) => bus( x(vplus,1,nk)       : )
   else
      hu       (1:n) => bus( x(humoins,1,nk)     : )
      ps       (1:n) => bus( x(pmoins,1,1)       : )
      th       (1:n) => bus( x(thetaa,1,1)       : )
      tt       (1:n) => bus( x(tmoins,1,nk)      : )
      uu       (1:n) => bus( x(umoins,1,nk)      : )
      vv       (1:n) => bus( x(vmoins,1,nk)      : )
   endif

   !------------------------------------------------------------------------


   ! 0.     Precompute local SST from "base" SST (ts) and diurnal components
   ! -----------------------------------------------------------------------
   sst = ts + zdsst + zskin_inc !add dirunal warm-layer and cool skin increments


   ! 1.     Saturated specific humidity at the water surface
   ! -------------------------------------------------------

   ! Uses FOQSA instead of FOQST to take into account saturation
   ! with respect to sea water (liquid between 0 and -1.8 C)
   qsat_o_salty = 1.0
   if (salty_qsat) qsat_o_salty = 0.98
   do I=1,N
      !        QS(I) = FOQSA(TS(I),PS(I)) * (1.-MLAC(I))*qsat_o_salty
      QS(I) = FOQSA(SST(I),PS(I))
      if(MLAC(I) .le. LAKE_MASK_THRESHOLD) QS(I) = qsat_o_salty*QS(I)
   end do


   ! 2.     Calculate roughness lengths based on generalized Charnock's relation
   ! ---------------------------------------------------------------------------

   if (kount.gt.0 .or. any('frv'==phyinread_list_s(1:phyinread_n))) then

      ! Precompute the Beljaars (QJRMS 1995, pp255-270) adjustment to the Charnock 
      ! relationship to increase the roughness length at low wind speeds (low ZFRV).
      select case (z0mtype)
      case ('BELJAARS') 
         Z0M_ADJUST = AM*K_VISC/max(ZFRV,ZFRVMIN)
      case ('CHARNOCK','WRF1','WRF2')
         Z0M_ADJUST = 0.
      case DEFAULT
         call physeterror('water', 'Unsupported z0mtype: '//trim(z0mtype))
         return
      end select

      ! Use Charnock's relation to estimate roughness length
      do I=1,N
         Z0M(I) = max( min( BETA1*ZFRV(I)**2/GRAV + Z0M_ADJUST(I),Z0MAX ) , Z0MIN )
      enddo

      ! WRF isftcflx=1/2
      if (z0mtype == 'WRF1' .or. z0mtype == 'WRF2') then
         do i=1,n
            zw = min(1., (zfrv(i)/1.06)**0.3)
            z1 = 0.011*zfrv(i)**2/GRAV + 1.59e-5
            z2 = 10./exp(9.5*max(zfrv(i),zfrvmin)**(-1./3.)) + 1.65e-6/max(zfrv(i),0.01)
            z0m(i) = max(1.27e-7, min(zw*z2 + (1-zw)*z1, 2.85e-3))
         enddo
      endif
      
   endif

   select case (z0ttype)
   case ('DEACU12')
      z0h = min(2.e-5/max(zfrv,ZFRVMIN),1.e-4)
   case ('ECMWF')
      z0h = min(6.e-6/max(zfrv,ZFRVMIN),5.e-4)
   case ('MOMENTUM')
      z0h = z0m
   case ('WRF1')
      z0h = 1e-4
   case ('WRF2')
      do i=1,n
         re = zfrv(i) * z0m(i)/k_visc
         z0h(i) = max(zfrv(i), Z0MIN) * exp(-karman*(7.3 * re**0.25 * 7.3**0.5 - 5.))
      enddo
   case DEFAULT
      call physeterror('water', 'Unsupported z0ttype: '//trim(z0ttype))
      return
   end select

   ! Note:  For |lat| >= Z0TLAT(2)  Charnock's (or Deacu) relation is used
   !        For |lat| <= Z0TLAT(1)  Z0HCON is used.
   !        For Z0TLAT(1) < |lat| < Z0TLAT(2)
   !        we do a linear interpolation between Charnock(or Deacu) and Z0HCON.

   do I=1,N

      if (abs(ZDLAT(I)) .ge. Z0TLAT(2)) then
         Z0H(I) = Z0H(I)
      else if (abs(ZDLAT(I)) .le. Z0TLAT(1)) then
         Z0H(I) = Z0HCON
      else
         Z0H(I)=( ((abs(ZDLAT(I))-Z0TLAT(1))/(Z0TLAT(2)-Z0TLAT(1))) &
              *(Z0H(I)-Z0HCON) ) + Z0HCON
      endif

   end do

   ! 3.     Calculate the surface transfer coefficient and fluxes
   ! ------------------------------------------------------------

   i = sl_prelim(tt,hu,uu,vv,ps,zzusl,rho_air=rho_a,spd_air=vmod0,dir_air=vdir,min_wind_speed=VAMIN)
   if (i /= SL_OK) then
      call physeterror('water', 'error returned by surface layer precalculations')
      return
   endif
   
   i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,sst,qs,z0m,z0h,zdlat,zfcor, &
        hghtm_diag=zu,hghtt_diag=zt,coefm=cmu,coeft=ctu,flux_t=zftemp, &
        flux_q=zfvap,ilmo=ilmo_wat,ue=zfrv,h=hst_wat,t_diag=ztdiag,    &
        q_diag=zqdiag,u_diag=zudiag,v_diag=zvdiag,tdiaglim=WATER_TDIAGLIM, &
        L_min=sl_Lmin_water,spdlim=vmod)
   if (i /= SL_OK) then
      call physeterror('water', 'error returned by surface layer calculations')
      return
   endif
   if (sl_Lmin_water > 0.) then
      zudiag = zudiag * vmod0 / vmod
      zvdiag = zvdiag * vmod0 / vmod
   endif
      
   ! Fill surface type-specific diagnostic values
   zqdiagtyp = zqdiag
   ztdiagtyp = ztdiag
   zudiagtyp = zudiag
   zvdiagtyp = zvdiag
   zqdiagtypv = zqdiag
   ztdiagtypv = ztdiag
   zudiagtypv = zudiag
   zvdiagtypv = zvdiag

   ! Also compute the air density at a chosen height (for now, the chosen height is
   ! at the screen level, i.e. at zt = 1.5m)

   i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,sst,qs,z0m,z0h,zdlat,zfcor, &
        hghtt_diag=1.5,t_diag=my_ta,q_diag=my_qa,tdiaglim=WATER_TDIAGLIM, &
        L_min=sl_Lmin_water,spdlim=vmod)
   if (i /= SL_OK) then
      call physeterror('water', 'error returned by surface layer density')
      return
  endif



   ! 4.     Finalize the fluxes
   ! --------------------------

!VDIR NODEP
   do I=1,N
      
      ZTSURF   (I) = SST (I)
      ZTSRAD   (I) = SST (I)
      
      ZALFAT   (I) = - CTU(I) * ( SST(I)-TH(I) )
      ZALFAQ   (I) = - CTU(I) * ( QS(I)-HU(I) )
      if (.not.IMPFLX) CTU (I) = 0.
      RHOA(i) = PS(I)/(RGASD * my_ta(I)*(1.+DELTA*my_qa(I)))
      !         RHOA(i) = PS(I)/(RGASD * ZTDIAG(I)*(1.+DELTA*ZQDIAG(I)))
      ! Compute runoff as total precipitation in kg/m2/s (or mm/s) minus evaporation.
      ! ZALFAQ being negative for an upward flux, we need to add rho*zalfaq    
      ZRUNOFFTOT (I) = (1000.*(ZRAINRATE(I) + ZSNOWRATE(I)) &
           + RHOA(I)*ZALFAQ(I)) * DELT
      FC_WAT(I)    = -CPD *RHOA(i)*ZALFAT(I)
      FV_WAT(I)    = -CHLC*RHOA(i)*ZALFAQ(I)

      if (IMPFLX) then
         ZALFAT   (I) = - CTU(I) *  SST(I)
         ZALFAQ   (I) = - CTU(I) *  QS(I)
      endif

      ! Set emissivity of water for the radiation scheme
      select case (water_emiss)
      case DEFAULT
         zemisr(i) = water_emiss_const
      end select
 
      ! Set albedo of water to predicted value if available
      if (update_alwater) alvis_wat(i) = zalwater(i)

   end do

   if (cplocn) then
      ! Update with fluxes and diagnostic variables from ocean model
      cplupd=.false.
      call cpl_update (vmod(1:n), 'UVO', lcl_indx, n, u=UU, v=VV, cplu=cplupd)
      if (cplupd) then
         call cpl_update (FC_WAT(1:n), 'SHO' , lcl_indx, n)
         call cpl_update (FV_WAT(1:n), 'LHO' , lcl_indx, n)
         call cpl_update (ZUDIAG(1:n), 'ZUO' , lcl_indx, n)
         call cpl_update (ZVDIAG(1:n), 'ZVO' , lcl_indx, n)
         call cpl_update (ZTDIAG(1:n), 'ZTO' , lcl_indx, n)
         call cpl_update (ZQDIAG(1:n), 'ZQO' , lcl_indx, n)
         call cpl_update (ZTSURF(1:n),   'TMW' , lcl_indx, n)
         call cpl_update (ZTSRAD(1:n),   'T4O' , lcl_indx, n)
         call cpl_update (QS(1:n),       'QSO' , lcl_indx, n)
         call cpl_update (ILMO_WAT(1:n), 'ILO' , lcl_indx, n)
         call cpl_update (Z0M   (1:n),   'ZMO' , lcl_indx, n)
         call cpl_update (Z0H   (1:n),   'ZHO' , lcl_indx, n)
         
         ! Derives other variables from cpl_update output
         ! assuming FC_WAT and FV_WAT were computed with
         ! RHOA at the same level
         if (zt_rho == zt) then
            my_ta(1:n)=ztdiag(1:n)
            my_qa(1:n)=zqdiag(1:n)
         else
            call physeterror('water', 'attempt to use an inconsistent density level')
            return
         endif
         do I=1,N
            RHOA(I) = PS(I)/(RGASD * my_ta(I)*(1.+DELTA*my_qa(I)))
            ZALFAT(I)= FC_WAT(I)/(-CPD *RHOA(I))
            ZALFAQ(I)= FV_WAT(I)/(-CHLC*RHOA(I))
            ZFTEMP(I) = -ZALFAT(I)
            ZFVAP(I)  = -ZALFAQ(I)
            SST(I)    = ZTSURF(I)
            if (IMPFLX) then
               ! CTU consistent with ocean model fluxes (as possible)
               ! Uncertainties in CTU from interpolation/agregation
               ! will be compensated by ZALFAT and ZALFAQ
               delh=SST(I)-TH(I) ; delq=QS(I)-HU(I)
               CTU(I) = 0.5*( -ZALFAT(I)/sign(max(abs(delh),delh_tresh),delh) &
                    -ZALFAQ(I)/sign(max(abs(delq),delh_tresh),delq) )
               CTU(I) = max(0.,CTU(I)) ! CTU<0 should never occur except under vicinity
               ! of null fluxes agregation/interpolation anyway
               ! e.g. sign(ZALFAT) /= sign(delh) or sign(ZALFAQ) /= sign(delq)
               ! ZALFAT and ZALFAQ consistent with FC_WAT and FV_WAT
               ZALFAT(I) = ZALFAT(I) - CTU(I)*TH(I)
               ZALFAQ(I) = ZALFAQ(I) - CTU(I)*HU(I)
            else
               CTU(I) = 0. !Redundant but less dangerous
            endif
         end do
         ! Diagnostic ustar based on ocean model stress (tau=rho*ustar**2)
         call cpl_update (ZFRV  (1:n), 'FRO' , lcl_indx, n, rho=RHOA, vmod=VMOD, cmu=CMU)
         ! Updated consistent CMU needed by implicit scheme using the following relation
         ! ==> CM=ustar/vmod (ustar**2=tau/rho=CM*CM*vmod**2)
      endif
   endif


   ! 5.     Prognostic evolution of SST (based on Zeng and Beljaars; GRL 2005)
   ! -------------------------------------------------------------------------

   DIURNAL_SST: if ( diusst == 'FAIRALL' ) then

      if (cplocn .and. diusst_ocean) then
         call physeterror('water', 'cplocn=true, diusst NOT taken care of yet')
         return
      endif

      ! Preliminary calculations
      i = sl_prelim(tt,hu,uu,vv,ps,zzusl,rho_air=rho_a,spd_air=vmod0,dir_air=vdir,min_wind_speed=VAMIN)
      
      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,sst,qs,z0m,z0h,zdlat,zfcor,hghtm_diag=10., &
           ue=frv_a,u_diag=ud,v_diag=vdi,tdiaglim=WATER_TDIAGLIM,L_min=sl_Lmin_water,spdlim=vmod)
      vmodd = sqrt(ud**2 + vdi**2)
      if (sl_Lmin_water > 0.) vmodd = vmodd * vmod0 / vmod
      alpha_w = max(8.75e-6*(sst-tcdk+9.),0.) !expansion coefficient estimate from Hayes et al. (JGR, 1991)
      visc_w = 1.8e-6 - 4.e-8*(sst-tcdk) !kinematic viscosity (from IITC 2011 recommendations at 10oC)
      frv_w = frv_a * sqrt(rho_a/rho_w) !friction velocity of water
      q_bal = -(fv_wat + fc_wat + (emis*stefan*sst**4 - zfdsi)) !surface energy balance (without shortwave)

      ! Cool skin parameterization 

      ! Estimate the depth and cooling effect of the cold skin layer (typically 1mm-1cm and 0.3K)
      ! Solar radiation absorbed in the skin layer
      where (zskin_depth < skin_depth_min)
         skin_solar = 0.1 !estimate proposed by Fairall et al. (JGR 1996) section 2.4
      elsewhere
         skin_solar = 0.067 + 11.*zskin_depth - (6.6e-5/zskin_depth) * (1.-exp(-zskin_depth/8.e-4))
      endwhere
      ! Updated skin layer depth (cold skin only exists when net heat is being extracted from the surface)
      skin_q = q_bal + (1.-alvis_wat)*zflusolis*skin_solar
      ! Compute Saunders' constant of proportionality (lambda)
      if (saunders == 'fairall') then !follow Fariall et al. (GRL; 1996)
         where (skin_q <= 0.)
            lambda = 6.*(1. + ((-16.*grav*alpha_w*visc_w**3)/(frv_w**4*diff_w**2*rho_w*cv_w) * (skin_q))**0.75)**(-1./3.)
         elsewhere
            lambda = 6.
         endwhere
      elseif (saunders == 'artale') then !follow Artale et al. (JGR; 2002)
         where (vmodd > 10.)
            gamma = 6.
         elsewhere (vmodd > 7.5)
            gamma = 1.6*vmodd - 10.
         elsewhere
            gamma = 0.2*vmodd + 0.5
         endwhere
         lambda = 86400.*frv_w*cond_w/(gamma*rho_w*cv_w*ref_depth*visc_w)
      else
         call physeterror('water', 'No Saunders constants defined for '//trim(saunders))
         return
      endif
      zskin_depth = lambda * visc_w / (sqrt(rho_a/rho_w)*frv_a)
      zskin_depth = max(min(zskin_depth,0.01),skin_depth_min) !limit skin depth to physical values
      ! Compute cool skin temperature increment (difference SST - T(zskin_depth))
      zskin_inc = zskin_depth / cond_w * (skin_q)
      if (saunders == 'fairall') zskin_inc = min(zskin_inc,0.)
      if (.not.diusst_lakes) where (mlac(:) > LAKE_MASK_THRESHOLD) zskin_inc(:) = 0.
      if (.not.diusst_ocean) where (mlac(:) < LAKE_MASK_THRESHOLD) zskin_inc(:) = 0.

      ! Warm layer parameterization

      ! Estimate the evolving temperature departure in the warm layer implicitly (secant iterations)
      this_inc = zdsst ; prev_inc = zdsst+10.*dsst_epsilon ; i = 0
      do while (maxval(abs(this_inc-prev_inc)) > dsst_epsilon .and. i < max_iter)
         i = i+1
         where (abs(this_inc-prev_inc) > dsst_epsilon)
            denom = warm_func(this_inc)-warm_func(prev_inc)
            where (abs(denom) > 1.e-8)
               ! Use local slope to estimate new increment value
               warm_increment = this_inc - warm_func(this_inc)*((this_inc-prev_inc)/denom)
               prev_inc = this_inc
               this_inc = warm_increment
            elsewhere
               ! Attempt a reset of the increment difference if the function is locally flat
               warm_increment = this_inc
               prev_inc = this_inc+10.*dsst_epsilon
            endwhere
         endwhere
      enddo
      if (diusst_lakes .and. diusst_warmlayer_lakes) &
           where (mlac(:) > LAKE_MASK_THRESHOLD) zdsst(:) = warm_increment(:)
      if (diusst_ocean) where (mlac(:) < LAKE_MASK_THRESHOLD) zdsst(:) = warm_increment(:)
      
   endif DIURNAL_SST

   !--------------------------------------
   !   6.     Heat Stress Indices
   !------------------------------------
   !#TODO: at least 4 times identical code in surface... separeted s/r to call
   IF_THERMAL_STRESS: if (thermal_stress) then

   ! Compute wind at z=zt
   i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,sst,qs,z0m,z0h,zdlat,zfcor, &
        hghtm_diag=zt,hghtt_diag=zt,u_diag=zusurfzt,v_diag=zvsurfzt, &
        tdiaglim=WATER_TDIAGLIM,L_min=sl_Lmin_water,spdlim=vmod)

   if (sl_Lmin_water > 0.) then
      zusurfzt = zusurfzt * vmod0 / vmod
      zvsurfzt = zvsurfzt * vmod0 / vmod
   endif

   if (i /= SL_OK) then
      call physeterror('water', 'error 3 returned by surface layer calculations')
      return
   endif

      do I=1,N

         if (abs(zzusl(i)-zu) <= 2.0) then
            zu10(i) = sqrt(uu(i)**2+vv(i)**2)
         else
            zu10(i) = sqrt(zudiag(i)**2+zvdiag(i)**2)
         endif

         ! wind  at SensoR level zubos at z=zt
         if( (abs(zusurfzt(i)) >= 0.1) .and. (abs(zvsurfzt(i)) >= 0.1)) then
         zusr(i) = sqrt( zusurfzt(i)**2 + zvsurfzt(i)**2)
         else
         zusr(i) = zu10(i)
         endif

            zqd(i) =max( ZQDIAG(i) , 1.e-6)  

         zref_sw_surf(i) = alvis_wat(i) * zflusolis(i)
         zemit_lw_surf(i)  = (1. - zemisr(i)) * zfdsi(i) + zemisr(i)*STEFAN   &
              *ztsurf(i)**4
         zzenith(i) = acos(zcoszeni(i))      
         if (zflusolis(i) > 0.0) then
            zzenith(i) = min(zzenith(i), pi/2.)
         else
            zzenith(i) = max(zzenith(i), pi/2.)
         endif

      end do

      call SURF_THERMAL_STRESS(ZTDIAG, zqd,            &
           ZU10,ZUSR,  ps,                             &
           ZFSD, ZFSF, ZFDSI, ZZENITH,                 &
           ZREF_SW_SURF,ZEMIT_LW_SURF,                 &
           Zutcisun ,Zutcishade,                       &
           zwbgtsun, zwbgtshade,                       &
           zradsun, zradshade,                         &
           ztglbsun, ztglbshade, ztwetb,               &
           ZQ1, ZQ2, ZQ3, ZQ4, ZQ5,                    &
           ZQ6,ZQ7, N)
   endif IF_THERMAL_STRESS
   !--------------------------------------

   ! FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
   call FILLAGG ( BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_WATER, SURFLEN )

   return

contains

   ! Warm layer function for secant iterations
   function warm_func(inc) result(fval)
      ! Argument declarations
      real, dimension(:), intent(in) :: inc             !warm layer increment
      real, dimension(size(inc)) :: fval                !warm layer function result
      ! Local variables
      real, dimension(size(inc)) :: base_bflux,mo_len,stab_func
      ! Compute buoyancy flux at base depth
      where (inc > 0.)
         base_bflux = sqrt(nu*grav*alpha_w/(5.*base_depth)) * rho_w*cv_w*frv_w**2*sqrt(inc)
      elsewhere
         base_bflux = grav*alpha_w*(q_bal+0.64*zflusolis)
      endwhere
      base_bflux = sign(max(abs(base_bflux),1.e-4),base_bflux)
      mo_len = rho_w*cv_w*frv_w**3/(karman*base_bflux) !Monin-Obukhov length in the warm layer
      mo_len = sign(max(abs(mo_len),1.e-4),mo_len)
      ! Compute stability function at base depth
      where (mo_len > 0.)
         stab_func = 1. + 5.*(base_depth/mo_len)
      elsewhere
         stab_func = (1. - 16.*(base_depth/mo_len))**(-0.5)
      endwhere
      ! Compute warm layer function value
      fval = inc - delt*((q_bal + 0.64*zflusolis) / (base_depth*rho_w*cv_w*nu/(nu+1)) - &
           (nu+1)*karman*frv_w/(base_depth*stab_func) * inc) - zdsst
   end function warm_func

end subroutine water2
