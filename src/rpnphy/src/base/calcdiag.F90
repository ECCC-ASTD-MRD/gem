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

module calcdiag
   implicit none
   private
   public :: calcdiag1

contains

   !/@*
   subroutine calcdiag1(pvars, dt, kount, ni, nk)
      !@Object Calculates averages and accumulators of tendencies and diagnostics
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use tdpack_const, only: CHLC, CPD, GRAV, TCDK, RAUW, CAPPA
      use integrals, only: int_profile, INT_OK
      use pbl_utils, only: blheight
      use sfclayer, only: sl_prelim, sl_sfclayer, SL_OK
      use phy_options
      use phybusidx
      use phymem, only: phyvar
      use phybudget, only: pb_compute, pb_residual
      use phy_status, only: PHY_OK
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - input -
      ! kount    timestep number
      ! dt       length of timestep
      ! n        horizontal running length
      ! nk       vertical dimension

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: kount, ni, nk
      real, intent(in) :: dt

      !@Author B. Bilodeau Feb 2003
      !*@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

      include "surface.cdk"
      include "physteps.cdk"

      real, parameter :: EC_Z0M_GRASS=0.03             !Threshold "flat grass" roughness for ECMWF diagnostics
      real, parameter :: EC_Z0T_GRASS=0.003            !Grass thermodynamic roughness for ECMWF diagnostics
      real, parameter :: EC_Z_ROUGH=40.                !Fixed height for ECMWF diagnostics in rough terrain
      real, parameter :: EC_MIN_LAND=0.1               !Minimum land fraction for soil-only ECMWF diagnostics
      character(len=*), parameter :: EC_INTERP='cubic' !Type of vertical interpolation for ECMWF diagnostics

      logical :: lmoyhr, laccum, lreset, lavg, lkount0, lacchr, is_mp, is_consun, &
           is_pcptype_nil, is_pcptype_sps, is_pcptype_b3d, is_radia, is_fluvert_nil, is_fluvert_sfc
      integer :: i, k, moyhr_steps, istat, istat1, nkm1
      real :: moyhri, tempo, tempo2, sol_stra, sol_conv, liq_stra, liq_conv, sol_mid, liq_mid
      real :: w1, w2, w_p3v5, w_etccdiag, w_my2
      real, dimension(ni) :: uvs, vmod, vdir, th_air, hblendm, ublend, &
           vblend, z0m_ec, z0t_ec, esdiagec, tsurfec, qsurfec, zrtrauw
      real(REAL64), dimension(ni) :: en0, pw0, en1, pw1
      real, dimension(ni,nk) :: presinv, t2inv, tiinv, lwcp, iwcp
      real, dimension(ni,nk-1) :: q_grpl, iiwc, prest

      real, target :: dummy1Dni(ni), dummy2Dnk(ni,nk)

#define PHYPTRDCL
#include "calcdiag_ptr.hf"

!# if (V1 >= V2) w=1 ; else w=0
#define MASK_GE(V1, V2) max(0., sign(1., V1 - V2))
!# if (V1 < V2) w=1 ; else w=0
#define MASK_LT(V1, V2) (1. - MASK_GE(V1, V2))
!# if (V1 <= V2) w=1 ; else w=0
#define MASK_LE(V1, V2) max(0., sign(1., V2 - V1))
!# if (V1 > V2) w=1 ; else w=0
#define MASK_GT(V1, V2) (1. - MASK_LE(V1, V2))
      
#define ISOUT1(ON)       (debug_alldiag_L .or. any(phystepoutlist_S(1:nphystepoutlist) == ON))
#define ISOUT2(ON,ONMOY) (debug_alldiag_L .or. any(phystepoutlist_S(1:nphystepoutlist) == ON) .or. any(phyoutlist_S(1:nphyoutlist) == ONMOY))

      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'calcdiag [BEGIN]')
      if (timings_L) call timing_start_omp(470, 'calcdiag', 46)

      nkm1 = nk-1

#undef PHYPTRDCL
#define PHYPTRASSIGNLOCAL
#include "calcdiag_ptr.hf"
      
      dummy1Dni = 0.
      dummy2Dnk = 0.
      
      call init2nan(uvs, vmod, vdir, th_air, hblendm, ublend)
      call init2nan(vblend, z0m_ec, z0t_ec, esdiagec, tsurfec, qsurfec, zrtrauw)
      call init2nan(en0, pw0, en1, pw1)
      call init2nan(presinv, t2inv, tiinv, q_grpl, iiwc, lwcp, iwcp)

      w_p3v5 = 0.
      if (stcond == 'MP_P3') w_p3v5 = 1.
      w_etccdiag = 0.
      if (etccdiag) w_etccdiag = 1.
      w_my2 = 0.
      if (stcond(1:6) == 'MP_MY2') w_my2 = 1.
      is_mp = (stcond(1:3) == 'MP_')
      is_consun = (stcond == 'CONSUN' .or. stcond == 'S2')
      is_pcptype_nil = (pcptype(1:3)=='NIL')
      is_pcptype_sps = any(pcptype == (/'SPS_W19', 'SPS_FRC'/))
      is_pcptype_b3d = (pcptype == 'BOURGE3D')
      is_radia = (radia(1:8) == 'CCCMARAD')
      is_fluvert_nil = (fluvert == 'NIL')
      is_fluvert_sfc = (fluvert == 'SURFACE')
      
      lkount0 = (kount == 0)
      lmoyhr = (moyhr > 0)
      lacchr = .false.
      if (acchr > 0) then
         lacchr = (mod(step_driver-1, acchr) == 0)
      elseif (acchr == 0) then
         lacchr = (step_driver-1 == 0)
      endif
      laccum = (lmoyhr .or. dynout)
      lreset = .false.
      lavg   = .false.
      moyhri = 1.
      if (lmoyhr) then
         lreset = (mod((step_driver-1),moyhr) == 0)
         lavg   = (mod(step_driver,moyhr) == 0)
         !# Compute the averaging interval (inverse), with all available data
         !  averaged when moving through driver step = 0
         moyhr_steps = moyhr
         if (step_driver == 0 .and. .not.lkount0) moyhr_steps = min(moyhr,kount)
         if (lavg) moyhri = 1./float(moyhr_steps)
      endif

      !# Final PBL height
      istat = blheight(zh,ztplus,zhuplus,zuplus,zvplus,zgzmom,zgztherm,zsigt, &
           ztsurf,zqsurf,zpplus,zz0_ag,zz0t_ag,zdlat,zfcor,ni,nk-1)

      ! Compute winds at the lowest thermodynamic level
      COMPUTE_SLT_WINDS: if (slt_winds) then         
         istat = sl_prelim(ztplus(:,nkm1),zhuplus(:,nkm1),zuplus(:,nkm1),zvplus(:,nkm1), &
              zpplus,zgzmom(:,nkm1),spd_air=vmod,dir_air=vdir,min_wind_speed=VAMIN)
         if (istat /= SL_OK) then
            call physeterror('calcdiag', 'Problem preparing surface layer calculations')
            return
         endif
         th_air = ztplus(:,nkm1)*zsigt(:,nkm1)**(-CAPPA)
         istat = sl_sfclayer(th_air,zhuplus(:,nkm1),vmod,vdir,zgzmom(:,nkm1), &
              zgztherm(:,nkm1),ztsurf,zqsurf,zz0_ag,zz0t_ag,zdlat,zfcor, &
              hghtm_diag_row=zgztherm(:,nkm1),u_diag=zuslt,v_diag=zvslt)
         if (istat /= SL_OK) then
            call physeterror('calcdiag', 'Problem with surface layer wind diagnostic')
            return
         endif
      endif COMPUTE_SLT_WINDS

      ! Sum dissipative wind tendencies from applicable sources (for SKEB)
      zudis(:,:) = zugwd(:,:) + zugno(:,:) + zufcp(:,:)
      zvdis(:,:) = zvgwd(:,:) + zvgno(:,:) + zvfcp(:,:)

      !****************************************************************
      !     Screen-level fields
      !     ---------------------------

      ! Derived screen-level fields
      if (.not.(lkount0 .and. .not.is_fluvert_nil) .or. is_fluvert_sfc) then

         ! Clip the screen level relative humidity to a range from 0-1
         call mfohr4(zrhdiag, zqdiag, ztdiag, zpplus, ni, 1, ni ,satuco)
         zrhdiag(1:ni) = max(min(zrhdiag(1:ni), 1.0), 0.)

         ! Screen level dewpoint depression
         call mhuaes3(zesdiag, zqdiag, ztdiag, zpplus, .false., ni, 1, ni)
         zesdiag(1:ni) = max(zesdiag(1:ni),0.)
         ztdew(1:ni) = ztdiag(1:ni) - zesdiag(1:ni)
         call mhuaes3(zesdiag, zqdiagstn, ztdiagstn, zpplus, .false., ni, 1, ni)
         zesdiag(1:ni) = max(zesdiag(1:ni),0.)
         ztddiagstn(1:ni) = ztdiagstn(1:ni) - zesdiag(1:ni)
         call mhuaes3(zesdiag, zqdiagstnv, ztdiagstnv, zpplus, .false., ni, 1, ni)
         zesdiag(1:ni) = max(zesdiag(1:ni),0.)
         ztddiagstnv(1:ni) = ztdiagstnv(1:ni) - zesdiag(1:ni)         

      endif

      ! ECMWF screen-level calculations performed on request
      ECMWF_SCREEN: if (ecdiag) then

         ! Prepare surface and height information for ECMWF diagnostic formulation
         where (zz0_ag(:) > EC_Z0M_GRASS)
            hblendm(:) = EC_Z_ROUGH
            z0m_ec(:) = EC_Z0M_GRASS
            z0t_ec(:) = EC_Z0T_GRASS
         elsewhere
            hblendm(:) = zgzmom(:,nkm1)
            z0m_ec(:) = zz0_ag(:)
            z0t_ec(:) = min(zz0t_ag(:),EC_Z0T_GRASS)
         endwhere
         where (zsfcwgt_soil(:) > .1)
            tsurfec(:) = ztsurf_soil(:)
            qsurfec(:) = zqsurf_soil(:)
         elsewhere
            tsurfec(:) = ztsurf(:)
            qsurfec(:) = zqsurf(:)
         endwhere

         ! Anemometer-level winds computed using the ECMWF formulation
         call vte_intvertx3(ublend, zuplus, zgzmom, hblendm, ni, nk, 1, 'UU', EC_INTERP)
         call vte_intvertx3(vblend, zvplus, zgzmom, hblendm, ni, nk, 1, 'VV', EC_INTERP)
         if ( sl_prelim(ztplus(:,nkm1), zhuplus(:,nkm1), ublend, vblend, zpplus, hblendm,  &
              spd_air=vmod, dir_air=vdir, min_wind_speed=VAMIN) /= SL_OK ) then
            call physeterror('calcdiag', 'Problem preparing EC screen-level calculations')
            return
         endif
         th_air(:) = ztplus(:,nkm1)*zsigt(:,nkm1)**(-CAPPA)
         if ( sl_sfclayer(th_air, zhuplus(:,nkm1), vmod, vdir, hblendm, zgztherm(:,nkm1),  &
              tsurfec, qsurfec, z0m_ec, z0t_ec, zdlat, zfcor, hghtt_diag=zt,               &
              hghtm_diag=zu, t_diag=ztdiagec, q_diag=zqdiagec, u_diag=zudiagec,            &
              v_diag=zvdiagec) /= SL_OK )  then
            call physeterror('calcdiag', 'Problem with EC screen-level diagnostic')
            return
         endif
         call mhuaes3(esdiagec, zqdiagec, ztdiagec, zpplus, .false., ni, 1, ni)
         do i=1,ni
            esdiagec(i) = max(esdiagec(i), 0.)
            ztddiagec(i) = ztdiagec(i) - esdiagec(i)
         enddo

      endif ECMWF_SCREEN
      
      !****************************************************************
      !     PRECIPITATION RATES AND ACCUMULATIONS
      !     -------------------------------------

      !# Set precipitation accumulators to zero at the beginning and after every
      !  acchr hours, and by default (acchr=0) as the model step goes through 0.
      IF_RESET_PRECIP: if (lkount0 .or. lacchr) then

         do i = 1, ni
            zasc(i)  = 0.
            zascs(i) = 0.
            zalc(i)  = 0.
            zalcm(i) = 0.
            zalcs(i) = 0.
            zass(i)  = 0.
            zals(i)  = 0.
            zpc(i)   = 0.
            zpy(i)   = 0.
            zpz(i)   = 0.
            zae(i)   = 0.
            zpr(i)   = 0.
            zazr(i)  = 0.
            zrn(i)   = 0.
            zaip(i)  = 0.
            zsn(i)   = 0.

            zals_rn1(i)  = 0.
            zals_rn2(i)  = 0.
            zals_fr1(i)  = 0.
            zals_fr2(i)  = 0.
            zass_sn1(i)  = 0.
            zass_sn2(i)  = 0.
            zass_sn3(i)  = 0.
            zass_pe1(i)  = 0.
            zass_pe2(i)  = 0.
            zass_pe2l(i)  = 0.
            zass_snd(i)   = 0.
            zass_mx(i)    = 0.
            zass_s2l(i)   = 0.
            zass_ws(i)    = 0.
            zsw_dhmax(i) = 0.
         enddo
      endif IF_RESET_PRECIP

      IF_KOUNT_NOT_0: if (.not.lkount0) then

         ! Diagnostics on precipitation type.
         if (pcptype == 'BOURGE' .or. (is_mp .and. is_pcptype_nil)) then
            !Note: 'bourge2' is always called if microphysics scheme is used and pcptyp=nil,
            !       but it is only applied to implicit (convection) precipitation
            call bourge2(zfneige1d, zfip1d, ztplus, zsigw, zpplus, ni, nk-1)
         else if (is_pcptype_b3d) then
            call bourge1_3d(zfneige2d, zfip2d, ztplus, zsigw, zpplus, ni, nk-1)
            zfneige1d(1:ni) => zfneige2d(:,nk)
            zfip1d(1:ni)    => zfip2d(:,nk)
         endif

         IF_BOURG3D: if (is_pcptype_b3d .and. is_consun) then
            !AZR3D: Accumulation des precipitations verglaclacantes en 3D
            !AIP3D: Accumulation des precipitations re-gelees en 3D
            do k = 1,nk-1
               do i = 1, ni
                  ! Flux de consun1
                  tempo    = (zsnoflx(i,k)+zrnflx(i,k))*.001
                  sol_stra = max(0.,zfneige2d(i,k)*tempo)
                  liq_stra = max(0.,tempo-sol_stra)
                  ! Flux de kfc
                  tempo    = (zkfcrf(i,k)+zkfcsf(i,k))*.001
                  sol_conv = max(0.,zfneige2d(i,k)*tempo)
                  liq_conv = max(0.,tempo-sol_conv)
                  ! Flux from mid-level convection
                  tempo    = (zkfmrf(i,k)+zkfmsf(i,k))*.001
                  sol_mid = max(0.,zfneige2d(i,k)*tempo)
                  liq_mid = max(0.,tempo-sol_mid)
                  tempo    = liq_stra+liq_conv+liq_mid
!!$                  if (ztplus(i,k) < TCDK) then
!!$                     zzr3d(i,k) = (1.-zfip2d(i,k))*tempo*dt
!!$                  else
!!$                     zzr3d(i,k) = 0.
!!$                  endif
                  w1 = MASK_LT(ztplus(i,k), TCDK)  !# ztplus(i,k) < TCDK
                  zzr3d(i,k) = w1 * ((1.-zfip2d(i,k))*tempo*dt)
                  zazr3d(i,k) = zazr3d(i,k) + zzr3d(i,k)
                  zpl3d(i,k) = zfip2d(i,k)*tempo*dt
                  zaip3d(i,k) = zaip3d(i,k) + zpl3d(i,k)
               enddo
            enddo
         endif IF_BOURG3D

!VDIR NODEP
         DO_NI: do i = 1,ni

            ! Running time maximum hail size (ice P3 microphysics)
            if (sw_dhmax > 0 .and. w_p3v5 > 0.) zsw_dhmax(i) = max(zsw_dhmax(i),a_diag_dhmax(i,nkm1))

            !taux des precipitations de la convection profonde
            zry(i) = ztsc(i) + ztlc(i)

            !taux des precipitations de la convection restreinte
            zrz(i) = ztscs(i) + ztlcs(i)

            !taux des precipitations liquides implicites
            zrlc(i) = ztlc(i) + ztlcm(i) + ztlcs(i)

            !taux des precipitations solides implicites
            zrsc(i) = ztsc(i) + ztscs(i)

            !taux total des precipitations implicites
            zrc(i) = zrlc(i) + zrsc(i)

            !Precipitation types from Milbrandt-Yau cloud microphysics scheme:
            IF_MP: if (is_mp)  then

               !diagnose freezing drizzle/rain based on sfc temperature
               !---override rn/fr partition from my2 [to be removed]
!!$               if (stcond(1:6) == 'MP_MY2') then
                  ztls_rn1(i) = ztls_rn1(i) + w_my2 * ztls_fr1(i)
                  ztls_rn2(i) = ztls_rn2(i) + w_my2 * ztls_fr2(i)
!!$               endif

               if  (ztplus(i,nk) < TCDK) then
                  ztls_fr1(i) = ztls_rn1(i)
                  ztls_fr2(i) = ztls_rn2(i)
                  ztls_rn1(i) = 0.
                  ztls_rn2(i) = 0.
               endif
!!$               w1 = MASK_LT(ztplus(i,k), TCDK)  !# ztplus(i,k) < TCDK
!!$               ztls_fr1(i) = (1.-w1) * ztls_fr1(i) + w1 * ztls_rn1(i)
!!$               ztls_fr2(i) = (1.-w1) * ztls_fr2(i) + w1 * ztls_rn2(i)
!!$               ztls_rn1(i) = (1.-w1) * ztls_rn1(i)
!!$               ztls_rn2(i) = (1.-w1) * ztls_rn2(i)

               !tls:  rate of liquid precipitation (sum of rain and drizzle)
               ztls(i) = ztls_rn1(i) + ztls_rn2(i) + ztls_fr1(i) + ztls_fr2(i)

               !tss:  rate of solid precipitation (sum of all fozen precipitation types)
               ztss(i) = ztss_sn1(i) + ztss_sn2(i) +ztss_sn3(i) + ztss_pe1(i) + ztss_pe2(i)
               ztss(i) = ztss(i) + w_p3v5 * ztss_ws(i)

               !tss_mx:  rate of mixed precipitation (liquid and solid simultaneously [>0.01 mm/h each])
!!$               if (ztls(i) > 2.78e-9 .and. ztss(i) > 2.78e-9) then   ![note: 2.78e-9 m/s = 0.01 mm/h]
!!$                  ztss_mx(i) = ztls(i) + ztss(i)
!!$               else
!!$                  ztss_mx(i) = 0.
!!$               endif
               w1 = MASK_GT(ztls(i), 2.78e-9)  !# ztls(i) > 2.78e-9
               w2 = MASK_GT(ztss(i), 2.78e-9)  !# ztss(i) > 2.78e-9
               ztss_mx(i) = w1 * w2 * (ztls(i) + ztss(i))


               !als_rn1:  accumulation of liquid drizzle
               zals_rn1(i) = zals_rn1(i) + ztls_rn1(i) * dt

               !als_rn2:  accumulation of liquid rain
               zals_rn2(i) = zals_rn2(i) + ztls_rn2(i) * dt

               !als_fr1:  accumulation of freezing drizzle
               zals_fr1(i) = zals_fr1(i) + ztls_fr1(i) * dt

               !als_fr2:  accumulation of freezing rain
               zals_fr2(i) = zals_fr2(i) + ztls_fr2(i) * dt

               !ass_sn1:  accumulation of ice crystals
               zass_sn1(i) = zass_sn1(i) + ztss_sn1(i) * dt

               !ass_sn2:  accumulation of snow
               zass_sn2(i) = zass_sn2(i) + ztss_sn2(i) * dt

               !ass_ws:  accumulation of wet snow
               zass_ws(i)  = zass_ws(i) + w_p3v5 * ztss_ws(i) * dt

               !ass_sn3:  accumulation of graupel
               zass_sn3(i) = zass_sn3(i) + ztss_sn3(i) * dt

               !ass_pe1:  accumulation of ice pellets
               zass_pe1(i) = zass_pe1(i) + ztss_pe1(i) * dt

               !ass_pe2:  accumulation of hail
               zass_pe2(i) = zass_pe2(i) + ztss_pe2(i) * dt

               !ass_pe2l:  accumulation of hail (large only)
               zass_pe2l(i) = zass_pe2l(i) + ztss_pe2l(i) * dt

               !ass_snd:  accumulation of total unmelted snow (i+s+g)
               zass_snd(i) = zass_snd(i) + ztss_snd(i) * dt

               !ass_mx:  accumulation of mixed precipitation
               zass_mx(i) = zass_mx(i) + ztss_mx(i) * dt

               !ass_s2l:  solid-to-liquid ratio for accumulated total "snow" (i+s+g)
               tempo =  max(zass_sn1(i) + zass_sn2(i) + zass_sn3(i), 1.e-18)
               zass_s2l(i) = zass_snd(i)/tempo

            endif IF_MP

            !taux des precipitations, grid-scale condensation scheme
            zrr(i) = ztss(i) + ztls(i)

            !taux total
            zrt(i) = zrc(i) + zrr(i)

            !asc : accumulation des precipitations solides de la convection profonde
            zasc(i) = zasc(i) + ztsc(i) * dt

            !ascs : accumulation des precipitations solides de la convection restreinte
            zascs(i) = zascs(i) + ztscs(i) * dt

            !alc : accumulation des precipitations liquides de la convection profonde
            zalc(i) = zalc(i) + ztlc(i) * dt

            !alcm : accumulation of liquid precipitation from mid-level convection
            zalcm(i) = zalcm(i) + ztlcm(i) * dt

            !alcs : accumulation des precipitations liquides de la convection restreinte
            zalcs(i) = zalcs(i) + ztlcs(i) * dt

            !ass : accumulation of total solid precipitation, grid-scale condensation scheme
            zass(i) = zass(i) + ztss(i) * dt

            !als : accumulation des precipitations liquides, grid-scale condensation scheme
            zals(i) = zals(i) + ztls(i) * dt

            !pc : accumulation des precipitations implicites
            zpc(i) = zalc(i) + zasc(i) + zalcs(i) + zascs(i) + zalcm(i)

            !py : accumulation des precipitations de la convection profonde
            zpy(i) = zalc(i) + zasc(i)

            !pz : accumulation des precipitations de la convection restreinte
            zpz(i) = zalcs(i) + zascs(i)

            !acm: accumulation des precipitations de la mid-level convective scheme
            zacm(i) = zalcm(i)

            !ae : accumulation des precipitations, grid-scale condensation scheme
            zae(i) = zals(i) + zass(i)

            !pr : accumulation des precipitations totales
            zpr(i) = zpc(i) + zae(i)

            if (is_consun) then
               tempo    = ztls(i) + ztss(i)
               sol_stra = max(0., zfneige1d(i)*tempo)
               liq_stra = max(0., tempo-sol_stra)
               tempo    = ztlc(i) + ztlcm(i) + ztlcs(i) + ztsc(i) + ztscs(i)
               sol_conv = max(0., zfneige1d(i)*tempo)
               liq_conv = max(0., tempo-sol_conv)
            elseif (is_pcptype_nil .or. is_pcptype_sps) then
               sol_stra = ztss(i)
               liq_stra = ztls(i)
               sol_conv = ztsc(i) + ztscs(i)
               liq_conv = ztlc(i) + ztlcm(i) + ztlcs(i)
            else  !i.e. IF stcond==mp_ .and. pcptype==bourge or bourge3d
               tempo    = ztls(i) + ztss(i)
               sol_stra = max(0., zfneige1d(i)*tempo)
               liq_stra = max(0., tempo-sol_stra)
               tempo    = ztlc(i) + ztlcm(i) + ztlcs(i) + ztsc(i) + ztscs(i)
               sol_conv = max(0., zfneige1d(i)*tempo)
               liq_conv = max(0., tempo-sol_conv)
            endif

            !RN : ACCUMULATION DES PRECIPITATIONS de PLUIE
            !AZR: ACCUMULATION DES PRECIPITATIONS VERGLACLACANTES
            !AIP: ACCUMULATION DES PRECIPITATIONS RE-GELEES
            !SN : ACCUMULATION DES PRECIPITATIONS de neige

            zzr(i) = 0.
            zzrflag(i) = 0.
            IF_MP_EXPL: if (is_mp .and. is_pcptype_nil) then

               tempo    = ztlc(i) + ztlcm(i) + ztlcs(i) + ztsc(i) + ztscs(i)   !total convective (implicit)
               sol_conv = max(0., zfneige1d(i)*tempo)
               liq_conv = max(0., tempo-sol_conv)
               tempo2   = zazr(i)
               !add explicit + diagnostic portion of rain/freezing rain from convective schemes:
!!$               w1 = MASK_LT(ztplus(i,nk), TCDK)  !# ztplus(i,nk) < TCDK
               if  (ztplus(i,nk) < TCDK) then
                  zazr(i) = zazr(i)                                  &
                       + (ztls_fr1(i)+ztls_fr2(i))*dt           & !from microphysics
                       + (1.-zfip1d(i))*liq_conv*dt               !from convective schemes
               else
                  zrn (i) = zrn(i)                                   &
                       + (ztls_rn1(i)+ztls_rn2(i))*dt           & !from microphysics
                       + (1.-zfip1d(i))*liq_conv*dt               !from convective schemes
               endif
               if (zazr(i) > tempo2) zzrflag(i) = 1.
!!$               w1 = MASK_GT(zazr(i), tempo2)
!!$               zzrflag(i) = (1.-w1) * zzrflag(i) + w1
               !add explicit + diagnostic portion of ice pellets from convective schemes:
               zaip(i) = zaip(i)                                     &
                    + ztss_pe1(i)*dt                            &  !from microphysics
                    + zfip1d(i)*tempo*dt                           !from convective schemes
               !note: Hail from M-Y (tss_pe2) is not included in total ice pellets (aip)
               !add explicit + diagnostic portion of snow from convective schemes:
               zsn(i)  = zsn(i)                                      &
                    + (ztss_sn1(i)+ztss_sn2(i)+ztss_sn3(i))*dt  &  !from microphysics
                    + sol_conv*dt                                  !from convective schemes
               zsn(i) = zsn(i) + w_p3v5 * ztss_ws(i)*dt
               !note: contribution to SN from microphysics is from ice, snow,
               !      and graupel, instantaneous solid-to-liquid ratio for
               !      pcp rate total snow (i+s+g):
               tempo       =  max(ztss_sn1(i)+ztss_sn2(i)+ztss_sn3(i)+ w_p3v5*ztss_ws(i), 1.e-18)
               ztss_s2l(i) = ztss_snd(i)/tempo

            else

               !diagnostic partitioning of rain vs. freezing rain:
               tempo =  liq_stra + liq_conv
!!$               w1 = MASK_LT(ztplus(i,nk), TCDK)  !# ztplus(i,nk) < TCDK
               if  (ztplus(i,nk) < TCDK) then
                  zzr(i) = (1.-zfip1d(i))*tempo*dt
                  zazr(i) = zazr(i) + zzr(i)
                  if (tempo > 0.) zzrflag(i) = 1.
               else
                  zrn(i) = zrn(i) + (1.-zfip1d(i))*tempo*dt
                  if (is_pcptype_b3d .and. tempo > 0.) zzrflag(i) = 1.
               endif
               !diagnostic ice pellets:
               zpl(i) = zfip1d(i)*tempo*dt
               zaip(i) = zaip(i) + zpl(i)
               !diagnostic snow:
               zsni(i) = (sol_stra+sol_conv)*dt
               zsn (i) = zsn (i) + zsni(i)

            endif IF_MP_EXPL

         enddo DO_NI

      endif IF_KOUNT_NOT_0

      if (lrefract) then
         call refractivity2(zdct_bh, zdct_count, zdct_lvl, zdct_lvlmax, &
              zdct_lvlmin, zdct_sndmax, zdct_sndmin, zdct_str, zdct_thick, &
              zdct_trdmax, zdct_trdmin, zdct_moref, &
              zpplus, zgztherm, zsigm, ztplus, zhuplus, ni, nk)
      endif

      if (llight) then
         if (stcond(1:5)=='MP_P3') then
            do k = 1, nk-1
               do i = 1, ni
                  q_grpl(i,k) = a_qi_4(i,k) + a_qi_5(i,k)
                  iiwc(i,k)   = a_qi_1(i,k) + a_qi_2(i,k) + a_qi_3(i,k) +     &
                       a_qi_4(i,k) + a_qi_5(i,k) + a_qi_6(i,k)
               enddo
            enddo
            call lightning2(zfoudre,zp0_plus,zsigm,ztplus,zwplus,q_grpl,iiwc,ni,nk)
         elseif (stcond(1:6)=='MP_MY2') then
            do k = 1, nk-1
               do i = 1, ni
                  q_grpl(i,k) = zqgplus(i,k)
                  iiwc(i,k)   = zqiplus(i,k) + zqnplus(i,k) + zqgplus(i,k)
               enddo
            enddo
            call lightning2(zfoudre,zp0_plus,zsigm,ztplus,zwplus,q_grpl,iiwc,ni,nk)
         endif
      endif

      !****************************************************************
      !     Energy budget diagnostics
      !     ---------------------------------

      zrtrauw = zrt*RAUW
      
      ! Compute Q1, Q2 apparent heat source / moisture sink on request
      Q1_BUDGET: if (ebdiag) then
         do k=1,nk
            do i = 1, ni
               presinv(i,nk-(k-1)) = zpplus(i)*zsigt(i,k)
               t2inv(i,nk-(k-1)) = zt2(i,k)
               tiinv(i,nk-(k-1)) = zti(i,k)
            enddo
         enddo
         istat  = int_profile(zt2i,t2inv,presinv,presinv(:,nk),presinv(:,1))
         istat1 = int_profile(ztii,tiinv,presinv,presinv(:,nk),presinv(:,1))
         if (istat /= INT_OK .or. istat1 /= INT_OK) then
            call physeterror('calcdiag', 'Problem in radiative tendency integrals')
            return
         endif
         
         do i = 1, ni
            zq1app(i) = CPD*(zt2i(i)+ztii(i))/GRAV + CHLC*zrtrauw(i) + zfc_ag(i)
            zq2app(i) = CHLC*(zrtrauw(i) - zflw(i))
         enddo
      endif Q1_BUDGET

      IF_CONS: if (lcons) then
         ! Compute physics budget residual
         if (conephy > 0 .or. conqphy > 0) then
            en0(:) = dble(zcone0(:))
            pw0(:) = dble(zconq0(:))
            zconephy = 0.
            if (pb_residual(zconephy, zconqphy, en0, pw0, pvars, &
                 delt, nkm1, F_rain=zrtrauw, F_shf=zfc_ag, F_wvf=zflw, &
                 F_rad=znetrad) /= PHY_OK) then
               call physeterror('calcdiag', &
                    'Problem computing physics budget residual')
               return
            endif
         endif

         ! Compute full-model budget residual
         if (conetot > 0 .or. conqtot > 0) then
            en1(:) = dble(zcone1(:))
            pw1(:) = dble(zconq1(:))
            zconetot = 0.
            if (pb_residual(zconetot, zconqtot, en1, pw1, pvars, &
                 delt, nkm1, F_rain=zrtrauw, F_shf=zfc_ag, F_wvf=zflw, &
                 F_rad=znetrad) /= PHY_OK) then
               call physeterror('calcdiag', &
                    'Problem computing full-model budget residual')
               return
            endif
         endif

         ! Compute post-physics budget state
         if (cone1 > 0 .or. conq1 > 0) then
            if (pb_compute(zcone1, zconq1, en1, pw1, pvars, nkm1) /= PHY_OK) then
               call physeterror('phystepinit', &
                    'Problem computing post-physics budget state')
               return
            endif
            if (cone1 > 0) zcone1(:) = real(en1(:))
            if (conq1 > 0) zconq1(:) = real(pw1(:))
         endif
         
      endif IF_CONS
   
      ! Compute surface energy budget (
      if (ISOUT2('FL', 'AG')) then
         zfl(:) = zfns(:) - zfv_ag(:) - zfc_ag(:)
      else
         zfl = 0.
      endif

      !****************************************************************
      !     Moisture diagnostics
      !     --------------------
      do k=1,nkm1
         prest(:,k) = zpplus(:) * zsigt(:,k) 
      enddo
      ! Calculate HRL (liquid), including diag level nk
      if (hrl > 0 .and. ISOUT1('HRL')) then
         call mfohr4(zhrl, zhuplus, ztplus, prest, ni, nkm1, ni, .false.)
         call mfohr4(zhrl(:,nk), zqdiag, ztdiag, zpplus, ni, 1, ni ,.false.)
      endif
      ! Calculate HRI or HRL (ice or liquid), including diag level nk
      if (hrli > 0 .and. ISOUT1('HRLI')) then
         call mfohr4(zhrli, zhuplus, ztplus, prest, ni, nkm1, ni, .true.)
         call mfohr4(zhrli(:,nk), zqdiag, ztdiag, zpplus, ni, 1, ni ,.true.)
      endif

      !****************************************************************
      !     AVERAGES
      !     --------

      !set averages to zero every moyhr hours
      IF_ACCUM_0: if (laccum) then

         !pre-calculate screen wind modulus
         uvs(1:ni) = sqrt(zudiag(1:ni)*zudiag(1:ni) + zvdiag(1:ni)*zvdiag(1:ni))

         IF_RESET_0: if (lkount0 .or. lreset) then

            do i = 1, ni
               zflwm(i)   = 0.
               zfshm(i)   = 0.
               zfcmy(i)   = 0.0
               zfvmm(i)   = 0.0
               zfvmy(i)   = 0.0
               zkshalm(i) = 0.0
               ziwvm(i)   = 0.0
               zq1appm(i) = 0.0
               zq2appm(i) = 0.0
               ztlwpm(i)  = 0.0
               zt2im(i)   = 0.0
               ztiim(i)   = 0.0
               ztiwpm(i)  = 0.0
               ztccm(i)  = 0.0
               ztslm(i)  = 0.0
               ztsmm(i)  = 0.0
               ztshm(i)  = 0.0
               ztzlm(i)  = 0.0
               ztzmm(i)  = 0.0
               ztzhm(i)  = 0.0
               zhrsmax(i) = zrhdiag(i)
               zhrsmin(i) = zrhdiag(i)
               zhusavg(i) = 0.0
               ztdiagavg(i)= 0.0
               zp0avg(i)  = 0.0
               zuvsavg(i) = 0.0
               zuvsmax(i) = uvs(i)
               zkmidm(i) = 0.0
            enddo

            !minimum and maximum temperature and temp. tendencies
            do k = 1, nk
               do i = 1, ni
                  zttmin(i,k) = ztplus(i,k)
                  zttmax(i,k) = ztplus(i,k)
                  ztadvmin(i,k) = ztadv(i,k)
                  ztadvmax(i,k) = ztadv(i,k)
               enddo
            enddo

            if (llinoz .and. out_linoz) then
               do i = 1, ni
                  ! 2D Ozone
                  zo3tcm(i) = 0.0
                  zo3ctcm(i) = 0.0
               enddo
               do k = 1, nk
                  do i = 1, ni
                     ! 3D Ozone
                     zo3avg (i,k) = 0.0
                     zo3ccolm(i,k) = 0.0
                     zo3colm(i,k) = 0.0
                     zo1chmtdm(i,k) = 0.0
                     zo4chmtdm(i,k) = 0.0
                     zo6chmtdm(i,k) = 0.0
                     zo7chmtdm(i,k) = 0.0
                     zo3chmtdm(i,k) = 0.0
                  enddo
               enddo
            end if

            if (llingh .and. out_linoz) then
               ! 2D GHG
               do i = 1, ni
                  zch4tcm(i) = 0.0
                  zn2otcm(i) = 0.0
                  zf11tcm(i) = 0.0
                  zf12tcm(i) = 0.0
               enddo
               ! 3D GHG
               do k = 1, nk
                  do i = 1, ni
                     zch4avg(i,k) = 0.0
                     zn2oavg(i,k) = 0.0
                     zf11avg(i,k) = 0.0
                     zf12avg(i,k) = 0.0
                     zch4colm(i,k) = 0.0
                     zn2ocolm(i,k) = 0.0
                     zf11colm(i,k) = 0.0
                     zf12colm(i,k) = 0.0
                     zch4chmtdm(i,k) = 0.0
                     zn2ochmtdm(i,k) = 0.0
                     zf11chmtdm(i,k) = 0.0
                     zf12chmtdm(i,k) = 0.0
                  enddo
               enddo
            end if

            !#TODO: split by scheme and avoid computation if scheme not active
            do k = 1, nk-1
               do i = 1, ni
                  zccnm(i,k)  = 0.0
                  ztim(i,k)   = 0.0
                  zt2m(i,k)   = 0.0
                  ztgwdm(i,k) = 0.0
                  zutofdm(i,k) = 0.0
                  zvtofdm(i,k) = 0.0
                  zttofdm(i,k) = 0.0
                  zugwdm(i,k) = 0.0
                  zvgwdm(i,k) = 0.0
                  zugnom(i,k) = 0.0
                  zvgnom(i,k) = 0.0
                  ztgnom(i,k) = 0.0
                  zudifvm(i,k) = 0.0
                  zvdifvm(i,k) = 0.0
                  ztdifvm(i,k) = 0.0
                  zqdifvm(i,k) = 0.0
                  ztadvm(i,k)  = 0.0
                  zuadvm(i,k)  = 0.0
                  zvadvm(i,k)  = 0.0
                  zqadvm(i,k)  = 0.0
                  zqmetoxm(i,k) = 0.0
                  zhushalm(i,k) = 0.0
                  ztshalm(i,k)  = 0.0
                  zlwcm(i,k)    = 0.0
                  ziwcm(i,k)    = 0.0
                  zlwcradm(i,k) = 0.0
                  ziwcradm(i,k) = 0.0
                  zcldradm(i,k) = 0.0
                  ztphytdm(i,k) = 0.0
                  zhuphytdm(i,k)= 0.0
                  zuphytdm(i,k) = 0.0
                  zvphytdm(i,k) = 0.0
                  zzctem(i,k)  = 0.0
                  zzstem(i,k)  = 0.0
                  zzcqem(i,k)  = 0.0
                  zzcqcem(i,k) = 0.0
                  zzsqem(i,k)  = 0.0
                  zzsqcem(i,k) = 0.0
                  zmqem(i,k) = 0.0
                  zmtem(i,k) = 0.0
                  zumim(i,k) = 0.0
                  zvmim(i,k) = 0.0
               enddo
            enddo

            !#Note: all convec str in the if below must have same len
            if (any(convec == (/ &
                 'KFC     ', &
                 'KFC2    ', &
                 'BECHTOLD'  &
                 /))) then
               do i = 1, ni
                  zabekfcm(i)  = 0.0
                  zcapekfcm(i) = 0.0
                  zcinkfcm(i)  = 0.0
                  zwumkfcm(i)  = 0.0
                  zzbaskfcm(i) = 0.0
                  zztopkfcm(i) = 0.0
                  zkkfcm(i)    = 0.0
               enddo

               do k = 1, nk
                  do i = 1, ni
                     ztfcpm(i,k)   = 0.0
                     zhufcpm(i,k)  = 0.0
                     zqckfcm(i,k)  = 0.0
                     zumfkfcm(i,k) = 0.0
                     zdmfkfcm(i,k) = 0.0
                     zudcm(i,k) = 0.0
                     zvdcm(i,k) = 0.0
                     zuscm(i,k) = 0.0
                     zvscm(i,k) = 0.0
                     zufcp1m(i,k) = 0.0
                     zufcp2m(i,k) = 0.0
                     zufcp3m(i,k) = 0.0
                     zsufcpm(i,k) = 0.0
                     zvfcp1m(i,k) = 0.0
                     zvfcp2m(i,k) = 0.0
                     zvfcp3m(i,k) = 0.0
                     zsvfcpm(i,k) = 0.0
                  enddo
               enddo
            
            endif

            if (.not.is_fluvert_nil) then
               do i = 1, ni
                  zsumf(i) = 0.
                  zsvmf(i) = 0.
                  zfqm(i)  = 0.
               enddo
            endif

            ! Reset conservation averages
            if (lcons) then
               do i = 1, ni
                  zconedynm(i) = 0.
                  zconecndm(i) = 0.
                  zconedcm(i) = 0.
                  zconemcm(i) = 0.
                  zconemom(i) = 0.
                  zconepblm(i) = 0.
                  zconephym(i) = 0.
                  zconeradm(i) = 0.
                  zconescm(i)  = 0.
                  zconetotm(i) = 0.
                  zconqdynm(i) = 0.
                  zconqcndm(i) = 0.
                  zconqdcm(i)  = 0.
                  zconqmcm(i)  = 0.
                  zconqmom(i)  = 0.
                  zconqpblm(i) = 0.
                  zconqphym(i) = 0.
                  zconqradm(i) = 0.
                  zconqscm(i) = 0.
                  zconqtotm(i) = 0.
               enddo
            endif

         endif IF_RESET_0

      endif IF_ACCUM_0


      IF_ACCUM_1: if (laccum .and. .not.lkount0) then

         do i = 1, ni
            dummy1Dni(i) = 0.
            zflwm    (i) = (zflwm   (i) + zflw(i) ) * moyhri
            zfshm    (i) = (zfshm   (i) + zfsh(i) ) * moyhri
            zfcmy    (i) = (zfcmy   (i) + zfc_ag(i)) * moyhri
            zfvmm    (i) = (zfvmm   (i) + zfvap_ag(i)) * moyhri
            zfvmy    (i) = (zfvmy   (i) + zfv_ag(i)) * moyhri
            zkshalm  (i) = (zkshalm (i) + zkshal(i)) * moyhri
            ziwvm    (i) = (ziwvm   (i) + ziwv  (i)) * moyhri
            zq1appm  (i) = (zq1appm (i) + zq1app(i)) * moyhri
            zq2appm  (i) = (zq2appm (i) + zq2app(i)) * moyhri
            ztlwpm   (i) = (ztlwpm  (i) + ztlwp (i)) * moyhri
            zt2im    (i) = (zt2im   (i) + zt2i  (i)) * moyhri
            ztiim    (i) = (ztiim   (i) + ztii  (i)) * moyhri
            ztiwpm   (i) = (ztiwpm  (i) + ztiwp (i)) * moyhri
            ztccm(i) = (w_etccdiag * (ztccm(i) + ztcc(i))) * moyhri
            ztslm(i) = (w_etccdiag * (ztslm(i) + ztcsl(i))) * moyhri
            ztsmm(i) = (w_etccdiag * (ztsmm(i) + ztcsm(i))) * moyhri
            ztshm(i) = (w_etccdiag * (ztshm(i) + ztcsh(i))) * moyhri
            ztzlm(i) = (w_etccdiag * (ztzlm(i) + ztczl(i))) * moyhri
            ztzmm(i) = (w_etccdiag * (ztzmm(i) + ztczm(i))) * moyhri
            ztzhm(i) = (w_etccdiag * (ztzhm(i) + ztczh(i))) * moyhri
            zhrsmax  (i) = max(zhrsmax  (i) , zrhdiag (i))
            zhrsmin  (i) = min(zhrsmin  (i) , zrhdiag (i))
            zhusavg  (i) = (    zhusavg  (i) + zqdiag  (i)) * moyhri
            ztdiagavg(i) = (    ztdiagavg(i) + ztdiag  (i)) * moyhri
            zp0avg   (i) = (    zp0avg   (i) + zpplus  (i)) * moyhri
            zuvsavg  (i) = (    zuvsavg  (i) + uvs     (i)) * moyhri
            zuvsmax  (i) = max(zuvsmax  (i) , uvs     (i))
         enddo

         ! Compute conservation averages on request
         CONSERVE_RESIDUAL_AVERAGING: if (lcons) then
            do i = 1, ni
               dummy1Dni(i) = 0.
               zconedynm(i) = (zconedynm(i) + zconedyn(i)) * moyhri
               zconecndm(i) = (zconecndm(i) + zconecnd(i)) * moyhri
               zconedcm(i)  = (zconedcm(i)  + zconedc(i))  * moyhri
               zconemcm(i)  = (zconemcm(i)  + zconemc(i))  * moyhri
               zconemom(i)  = (zconemom(i)  + zconemo(i))  * moyhri
               zconepblm(i) = (zconepblm(i) + zconepbl(i)) * moyhri
               zconephym(i) = (zconephym(i) + zconephy(i)) * moyhri
               zconeradm(i) = (zconeradm(i) + zconerad(i)) * moyhri
               zconescm(i)  = (zconescm(i)  + zconesc(i))  * moyhri
               zconetotm(i) = (zconetotm(i) + zconetot(i)) * moyhri
               zconqdynm(i) = (zconqdynm(i) + zconqdyn(i)) * moyhri
               zconqcndm(i) = (zconqcndm(i) + zconqcnd(i)) * moyhri
               zconqdcm(i)  = (zconqdcm(i)  + zconqdc(i))  * moyhri
               zconqmcm(i)  = (zconqmcm(i)  + zconqmc(i))  * moyhri
               zconqmom(i)  = (zconqmom(i)  + zconqmo(i))  * moyhri
               zconqpblm(i) = (zconqpblm(i) + zconqpbl(i)) * moyhri
               zconqphym(i) = (zconqphym(i) + zconqphy(i)) * moyhri
               zconqradm(i) = (zconqradm(i) + zconqrad(i)) * moyhri
               zconqscm(i)  = (zconqscm(i)  + zconqsc(i))  * moyhri
               zconqtotm(i) = (zconqtotm(i) + zconqtot(i)) * moyhri
            enddo
         endif CONSERVE_RESIDUAL_AVERAGING

         !# minimum and maximum temperature and temp.tendencies
         do k = 1, nk
            do i = 1, ni
               zttmin  (i,k) = min(zttmin  (i,k), ztplus(i,k))
               zttmax  (i,k) = max(zttmax  (i,k), ztplus(i,k))
               ztadvmin(i,k) = min(ztadvmin(i,k), ztadv (i,k))
               ztadvmax(i,k) = max(ztadvmax(i,k), ztadv (i,k))
            end do
         enddo

         IF_LINOZ: if (llinoz .and. out_linoz) then

            do i = 1, ni
               zo3tcm  (i) = (zo3tcm  (i) + zo3tc  (i)) * moyhri
               zo3ctcm (i) = (zo3ctcm (i) + zo3ctc (i)) * moyhri
            end do

            do k = 1, nk-1
               do i = 1, ni
                  dummy2Dnk(i,k) = 0.
                  zo3avg  (i,k) = (zo3avg  (i,k) + zo3lplus (i,k)) * moyhri

                  zo3ccolm (i,k) = (zo3ccolm (i,k) + zo3ccol (i,k)) * moyhri
                  zo3colm  (i,k) = (zo3colm  (i,k) + zo3col  (i,k)) * moyhri

                  zo1chmtdm  (i,k) = ( zo1chmtdm (i,k) +  zo1chmtd (i,k)) * moyhri
                  zo4chmtdm  (i,k) = ( zo4chmtdm (i,k) +  zo4chmtd (i,k)) * moyhri
                  zo6chmtdm  (i,k) = ( zo6chmtdm (i,k) +  zo6chmtd (i,k)) * moyhri
                  zo7chmtdm  (i,k) = ( zo7chmtdm (i,k) +  zo7chmtd (i,k)) * moyhri
                  zo3chmtdm  (i,k) = ( zo3chmtdm (i,k) +  zo3chmtd (i,k)) * moyhri
               end do
            end do
            
         end if IF_LINOZ

         IF_LINGH: if (llingh .and. out_linoz) then
            ! 2D
            do i = 1, ni
               dummy1Dni(i) = 0.
               zch4tcm (i) = (zch4tcm (i) + zch4tc (i)) * moyhri
               zn2otcm (i) = (zn2otcm (i) + zn2otc (i)) * moyhri
               zf11tcm (i) = (zf11tcm (i) + zf11tc (i)) * moyhri
               zf12tcm (i) = (zf12tcm (i) + zf12tc (i)) * moyhri
            end do

            ! 3D
            do k = 1, nk-1
               do i = 1, ni
                  dummy2Dnk(i,k) = 0.
                  zch4avg (i,k) = (zch4avg (i,k) + zch4lplus (i,k)) * moyhri
                  zn2oavg (i,k) = (zn2oavg (i,k) + zn2olplus (i,k)) * moyhri
                  zf11avg (i,k) = (zf11avg (i,k) + zf11lplus (i,k)) * moyhri
                  zf12avg (i,k) = (zf12avg (i,k) + zf12lplus (i,k)) * moyhri

                  zch4colm (i,k) = (zch4colm (i,k) + zch4col (i,k)) * moyhri
                  zn2ocolm (i,k) = (zn2ocolm (i,k) + zn2ocol (i,k)) * moyhri
                  zf11colm (i,k) = (zf11colm (i,k) + zf11col (i,k)) * moyhri
                  zf12colm (i,k) = (zf12colm (i,k) + zf12col (i,k)) * moyhri

                  zch4chmtdm (i,k) = (zch4chmtdm (i,k) + zch4chmtd (i,k)) * moyhri
                  zn2ochmtdm (i,k) = (zn2ochmtdm (i,k) + zn2ochmtd (i,k)) * moyhri
                  zf11chmtdm (i,k) = (zf11chmtdm (i,k) + zf11chmtd (i,k)) * moyhri
                  zf12chmtdm (i,k) = (zf12chmtdm (i,k) + zf12chmtd (i,k)) * moyhri
               end do
            end do
         end if IF_LINGH

         !#TODO: split by scheme and avoid computation if scheme not active
         !#TODO: Make it conditional to reqested var
         do k = 1, nk-1
            do i = 1, ni
               dummy2Dnk(i,k) = 0.
               zccnm   (i,k) = (zccnm   (i,k) + zftot(i,k)) * moyhri
               ztim    (i,k) = (ztim    (i,k) + zti  (i,k)) * moyhri
               zt2m    (i,k) = (zt2m    (i,k) + zt2  (i,k)) * moyhri

               zudifvm (i,k) = (zudifvm (i,k) + zudifv(i,k)) * moyhri
               zvdifvm (i,k) = (zvdifvm (i,k) + zvdifv(i,k)) * moyhri
               ztdifvm (i,k) = (ztdifvm (i,k) + ztdifv(i,k)) * moyhri
               zqdifvm (i,k) = (zqdifvm (i,k) + zqdifv(i,k)) * moyhri
               
               ztadvm  (i,k) = (ztadvm  (i,k) + ztadv (i,k)) * moyhri
               zuadvm  (i,k) = (zuadvm  (i,k) + zuadv (i,k)) * moyhri
               zvadvm  (i,k) = (zvadvm  (i,k) + zvadv (i,k)) * moyhri
               zqadvm  (i,k) = (zqadvm  (i,k) + zqadv (i,k)) * moyhri
               
               zqmetoxm(i,k) = (zqmetoxm(i,k) + zqmetox(i,k)) * moyhri
               
               zhushalm(i,k) = (zhushalm(i,k) + zhushal(i,k)) * moyhri
               ztshalm (i,k) = (ztshalm (i,k) + ztshal(i,k)) * moyhri
               
               zlwcm   (i,k) = (zlwcm   (i,k) + zlwc(i,k)) * moyhri
               ziwcm   (i,k) = (ziwcm   (i,k) + ziwc(i,k)) * moyhri
               zlwcradm(i,k) = (zlwcradm(i,k) + zlwcrad(i,k)) * moyhri
               ziwcradm(i,k) = (ziwcradm(i,k) + ziwcrad(i,k)) * moyhri
               zcldradm(i,k) = (zcldradm(i,k) + zcldrad(i,k)) * moyhri
               
               ztphytdm(i,k) = (ztphytdm(i,k) + ztphytd(i,k)) * moyhri
               zhuphytdm(i,k)= (zhuphytdm(i,k)+ zhuphytd(i,k)) * moyhri
               zuphytdm(i,k) = (zuphytdm(i,k) + zuphytd(i,k)) * moyhri
               zvphytdm(i,k) = (zvphytdm(i,k) + zvphytd(i,k)) * moyhri

               !#deep
               zzctem(i,k)   = (zzctem(i,k)   + zcte(i,k)) * moyhri
               zzstem(i,k)   = (zzstem(i,k)   + zste(i,k)) * moyhri
               zzcqem(i,k)   = (zzcqem(i,k)   + zcqe(i,k)) * moyhri
               zzsqem(i,k)   = (zzsqem(i,k)   + zsqe(i,k)) * moyhri
               zzcqcem(i,k)  = (zzcqcem(i,k)  + zcqce(i,k)) * moyhri
               zzsqcem(i,k)  = (zzsqcem(i,k)  + zsqce(i,k)) * moyhri
            end do
         end do
         
         if (tofd /= 'NIL') then
            do k = 1, nk-1
               do i = 1, ni
                  zutofdm (i,k) = (zutofdm (i,k) + zutofd(i,k)) * moyhri
                  zvtofdm (i,k) = (zvtofdm (i,k) + zvtofd(i,k)) * moyhri
                  zttofdm (i,k) = (zttofdm (i,k) + zttofd(i,k)) * moyhri
               enddo
            enddo
         endif

         if (gwdrag /= 'NIL') then
            do k = 1, nk-1
               do i = 1, ni
                  ztgwdm (i,k)  = (ztgwdm (i,k)  + ztgwd(i,k)) * moyhri
                  zugwdm  (i,k) = (zugwdm  (i,k) + zugwd(i,k)) * moyhri
                  zvgwdm  (i,k) = (zvgwdm  (i,k) + zvgwd(i,k)) * moyhri
                  zugnom  (i,k) = (zugnom  (i,k) + zugno(i,k)) * moyhri
                  zvgnom  (i,k) = (zvgnom  (i,k) + zvgno(i,k)) * moyhri
                  ztgnom  (i,k) = (ztgnom  (i,k) + ztgno(i,k)) * moyhri
               enddo
            enddo
         endif

         if (mqem > 0) then
            do k = 1, nk-1
               do i = 1, ni
                  zmqem(i,k) = (zmqem(i,k) + zmqe(i,k) ) * moyhri
                  zmtem(i,k) = (zmtem(i,k) + zmte(i,k) ) * moyhri
                  zumim(i,k) = (zumim(i,k) + zumid(i,k) ) * moyhri
                  zvmim(i,k) = (zvmim(i,k) + zvmid(i,k) ) * moyhri
               enddo
            enddo
         endif

         !# Note: all convec str in the if below must have same len
         if (any(convec == (/ &
              'KFC     ', &
              'KFC2    ', &
              'BECHTOLD'  &
              /))) then
            do k = 1, nk-1
               do i = 1, ni
                  ztfcpm  (i,k) = (ztfcpm  (i,k) + ztfcp (i,k)) * moyhri
                  zhufcpm (i,k) = (zhufcpm (i,k) + zhufcp(i,k)) * moyhri
                  zqckfcm (i,k) = (zqckfcm (i,k) + zqckfc(i,k)) * moyhri
                  zumfkfcm(i,k) = (zumfkfcm(i,k) + zumfkfc(i,k)) * moyhri
                  zdmfkfcm(i,k) = (zdmfkfcm(i,k) + zdmfkfc(i,k)) * moyhri
                  zudcm(i,k)    = (zudcm(i,k)   + zufcp(i,k)) * moyhri
                  zvdcm(i,k)    = (zvdcm(i,k)   + zvfcp(i,k)) * moyhri
               enddo
            enddo
            if (cmt_comp_diag) then
               do k = 1, nk-1
                  do i = 1, ni
                     zufcp1m(i,k) = (zufcp1m(i,k) + zufcp1(i,k)) * moyhri
                     zufcp2m(i,k) = (zufcp2m(i,k) + zufcp2(i,k)) * moyhri
                     zufcp3m(i,k) = (zufcp3m(i,k) + zufcp3(i,k)) * moyhri
                     zsufcpm(i,k) = (zsufcpm(i,k) + zsufcp(i,k)) * moyhri
                     zvfcp1m(i,k) = (zvfcp1m(i,k) + zvfcp1(i,k)) * moyhri
                     zvfcp2m(i,k) = (zvfcp2m(i,k) + zvfcp2(i,k)) * moyhri
                     zvfcp3m(i,k) = (zvfcp3m(i,k) + zvfcp3(i,k)) * moyhri
                     zsvfcpm(i,k) = (zsvfcpm(i,k) + zsvfcp(i,k)) * moyhri
                  enddo
               enddo
            endif
            
            if (tusc > 0 .and. uscm > 0) &
                 zuscm(:,1:nk-1) = (zuscm(:,1:nk-1) + ztusc(:,1:nk-1)) * moyhri
            if (tvsc > 0 .and. vscm > 0) &
                 zvscm(:,1:nk-1) = (zvscm(:,1:nk-1) + ztvsc(:,1:nk-1)) * moyhri
            
            do i=1, ni
               dummy1Dni(i) = 0.
               zabekfcm  (i) = (zabekfcm(i) + zabekfc(i)) * moyhri
               zcapekfcm (i) = (zcapekfcm(i) + zcapekfc(i)) * moyhri
               zcinkfcm  (i) = (zcinkfcm(i) + zcinkfc(i)) * moyhri
               zwumkfcm  (i) = (zwumkfcm(i) + zwumaxkfc(i)) * moyhri
               zzbaskfcm (i) = (zzbaskfcm(i) + zzbasekfc(i)) * moyhri
               zztopkfcm (i) = (zztopkfcm(i) + zztopkfc(i)) * moyhri
               zkkfcm    (i) = (zkkfcm(i) + zkkfc(i)) * moyhri
               if (kmidm > 0) zkmidm(i) = (zkmidm(i) + zkmid(i)) * moyhri
            end do
         endif

         if (.not.is_fluvert_nil) then
            do i=1,ni
               dummy1Dni(i) = 0.
               zsumf(i) = (zsumf(i) + zustress(i)) * moyhri
               zsvmf(i) = (zsvmf(i) + zvstress(i)) * moyhri
               zfqm(i)  = (zfqm(i) + zfq(i)) * moyhri
            enddo
         endif

      endif IF_ACCUM_1


      !****************************************************************
      !     ACCUMULATORS

      !Set accumulators to zero at the beginning and after every acchr hours,
      !and by default (acchr=0) as the model step goes through 0.
      IF_RESET_ACCUMULATORS: if (lkount0 .or. lacchr) then
         do i = 1, ni
            zrainaf(i) = 0.
            zsnowaf(i) = 0.
            zeiaf(i) = 0.
            zevaf(i) = 0.
            zfiaf(i) = 0.
            zfsaf(i) = 0.
            zivaf(i) = 0.
            zntaf(i) = 0.
            zflusolaf(i) = 0.
            zclbaf(i) = 0.
            zcltaf(i) = 0.
            zcstaf(i) = 0.
            zcsbaf(i) = 0.
            zfsdaf(i) = 0.
            zfsfaf(i) = 0.
            zfsiaf(i) = 0.
            zfsvaf(i) = 0.
            zparraf(i) = 0.
            zsiaf(i) = 0.
            zfnsaf(i) = 0.
            zflaf(i) = 0.
            zfcaf(i) = 0.
            zfvaf(i) = 0.
            zafoudre(i) = 0.
         enddo
      endif IF_RESET_ACCUMULATORS

      IF_KOUNT_NOT_0b: if (.not.lkount0) then
!VDIR NODEP
         DO_NI_ACC: do i = 1,ni
            dummy1Dni(i) = 0.

            !# Accumulation of precipitation (in m)

            zrainaf(i) = zrainaf(i) + zrainrate(i)*dt
            zsnowaf(i) = zsnowaf(i) + zsnowrate(i)*dt

            if (is_radia .or. is_fluvert_sfc) then
               zeiaf    (i) = zeiaf (i) + zei  (i) * dt
               zevaf    (i) = zevaf (i) + zev  (i) * dt
               zfiaf    (i) = zfiaf (i) + zfdsi(i) * dt
               zfsaf    (i) = zfsaf (i) + zfdss(i) * dt
               zivaf    (i) = zivaf (i) + ziv  (i) * dt
               zntaf    (i) = zntaf (i) + znt  (i) * dt
               zflusolaf(i) = zflusolaf(i) + zflusolis(i) * dt
            endif

            !# Accumulation of sfc and toa net clear sky fluxes, available with cccmarad
            if (is_radia) then
               zclbaf    (i) = zclbaf (i) + zclb  (i) * dt
               zcltaf    (i) = zcltaf (i) + zclt  (i) * dt
               zcstaf    (i) = zcstaf (i) + zcstt (i) * dt
               zcsbaf    (i) = zcsbaf (i) + zcsb  (i) * dt
               zfsdaf    (i) = zfsdaf (i) + zfsd  (i) * dt
               zfsfaf    (i) = zfsfaf (i) + zfsf  (i) * dt
               zfsiaf    (i) = zfsiaf (i) + zfsi  (i) * dt
               zfsvaf    (i) = zfsvaf (i) + zfsv  (i) * dt
               zparraf   (i) = zparraf(i) + zparr (i) * dt
               zsiaf (i) = zsiaf (i) + zfnsi(i) * dt
               zfnsaf (i) = zfnsaf (i) + zfns(i) * dt
            endif

            if (.not.is_fluvert_nil) then
               zflaf (i) = zflaf (i) + zfl  (i) * dt
               zfcaf (i) = zfcaf (i) + zfc_ag(i) * dt
               zfvaf (i) = zfvaf (i) + zfv_ag(i) * dt
            endif

            !# Accumulation of lightning threat (in number of flashes/m2)
            if (llight) then
               zafoudre(i) = zafoudre(i) + zfoudre(i)*dt
            endif

         end do DO_NI_ACC
      endif IF_KOUNT_NOT_0B

      !# For output purpose, diag level values needs to be copied into nk level of corresponding dyn var
      do i = 1,ni
         zhuplus(i,nk) = zqdiag(i)
         ztplus(i,nk)  = ztdiag(i)
         zuplus(i,nk)  = zudiag(i)
         zvplus(i,nk)  = zvdiag(i)
      enddo
      
      if (timings_L) call timing_stop_omp(470)
      call msg_toall(MSG_DEBUG, 'calcdiag [END]')
      !----------------------------------------------------------------
      return
   end subroutine calcdiag1

end module calcdiag
