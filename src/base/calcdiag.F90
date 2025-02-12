
#define DO_AVG(A,V,moyhri,WRESET) if (associated(A)) then ; A = (A*WRESET + V) * moyhri ; endif
#define DO_ACC(A,V,dt,WRESET)     if (associated(A)) then ; A = A*WRESET + V*dt ; endif
#define DO_MIN(A,V,WRESET)        if (associated(A)) then ; A = min(A*WRESET + V*(1.-WRESET), V); endif
#define DO_MAX(A,V,WRESET)        if (associated(A)) then ; A = max(A*WRESET + V*(1.-WRESET), V); endif


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
           is_pcptype_nil, is_pcptype_sps, is_pcptype_b3d, is_fluvert_nil, is_fluvert_sfc, is_hrl, is_hrli
      integer :: i, k, moyhr_steps, istat, istat1, nkm1
      real :: moyhri, tempo, tempo2, sol_stra, sol_conv, liq_stra, liq_conv, sol_mid, liq_mid
      real :: w1, w2, w_p3v5, w_my2, w_reset
      real, dimension(ni) :: uvs, vmod, vdir, th_air, hblendm, ublend, &
           vblend, z0m_ec, z0t_ec, esdiagec, tsurfec, qsurfec, zrtrauw
      real(REAL64), dimension(ni) :: en0, pw0, en1, pw1
      real, dimension(ni,nk) :: presinv, t2inv, tiinv, lwcp, iwcp
      real, dimension(ni,nk-1) :: q_grpl, iiwc, prest

!!$      real, target :: dummy1Dni(ni), dummy2Dnk(ni,nk)

#define PHYPTRDCL
#include "calcdiag_ptr.hf"

      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'calcdiag [BEGIN]')
      if (timings_L) call timing_start_omp(470, 'calcdiag', 46)

      nkm1 = nk-1

#undef PHYPTRDCL
#include "calcdiag_ptr.hf"
      
      call init2nan(uvs, vmod, vdir, th_air, hblendm, ublend)
      call init2nan(vblend, z0m_ec, z0t_ec, esdiagec, tsurfec, qsurfec, zrtrauw)
      call init2nan(en0, pw0, en1, pw1)
      call init2nan(presinv, t2inv, tiinv, q_grpl, iiwc, lwcp, iwcp)

      w_p3v5 = 0.
      if (stcond == 'MP_P3') w_p3v5 = 1.
      w_my2 = 0.
      if (stcond(1:6) == 'MP_MY2') w_my2 = 1.
      is_mp = (stcond(1:3) == 'MP_')
      is_consun = (stcond == 'CONSUN' .or. stcond == 'S2')
      is_pcptype_nil = (pcptype(1:3)=='NIL')
      is_pcptype_sps = any(pcptype == (/'SPS_W19', 'SPS_FRC'/))
      is_pcptype_b3d = (pcptype == 'BOURGE3D')
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
      istat = blheight(zhpbl,ztplus,zhuplus,zuplus,zvplus,zgzmom,zgztherm,zsigt, &
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
      if (associated(zudis)) zudis(:,:) = zugwd(:,:) + zugno(:,:) + zufcp(:,:)
      if (associated(zvdis)) zvdis(:,:) = zvgwd(:,:) + zvgno(:,:) + zvfcp(:,:)

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
      w_reset = 1.
      if (lkount0 .or. lacchr) w_reset = 0.
      
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
            if (sw_dhmax > 0 .and. w_p3v5 > 0.) zsw_dhmax(i) = max(zsw_dhmax(i)*w_reset,a_diag_dhmax(i,nkm1))

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

            endif IF_MP

            !taux des precipitations, grid-scale condensation scheme
            zrr(i) = ztss(i) + ztls(i)

            !taux total
            zrt(i) = zrc(i) + zrr(i)

            if (is_pcptype_nil .or. is_pcptype_sps) then
               sol_stra = ztss(i)
               liq_stra = ztls(i)
               sol_conv = ztsc(i) + ztscs(i)
               liq_conv = ztlc(i) + ztlcm(i) + ztlcs(i)
            else
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
            zrn(i) = zrn(i)*w_reset
            zazr(i) = zazr(i)*w_reset
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
               zaip(i) = zaip(i)*w_reset                        &
                    + ztss_pe1(i)*dt                            &  !from microphysics
                    + zfip1d(i)*tempo*dt                           !from convective schemes
               !note: Hail from M-Y (tss_pe2) is not included in total ice pellets (aip)
               !add explicit + diagnostic portion of snow from convective schemes:
               zsn(i)  = zsn(i)*w_reset                         &
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
               zpl(i)  = zfip1d(i)*tempo*dt
               zaip(i) = zaip(i)*w_reset + zpl(i)
               !diagnostic snow:
               zsni(i) = (sol_stra+sol_conv)*dt
               zsn(i)  = zsn(i)*w_reset + zsni(i)

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
            en1(:) = dble(zcone1(:))
            pw1(:) = dble(zconq1(:))
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
      if (ISREQOUT('AG') .or. ISREQSTEP('FL')) then
         zfl(:) = zfns(:) - zfv_ag(:) - zfc_ag(:)
      else
         zfl = 0.
      endif

      !****************************************************************
      !     Moisture diagnostics
      !     --------------------

      is_hrl  = (hrl > 0  .and. ISREQSTEP('HRL'))
      is_hrli = (hrli > 0 .and. ISREQSTEP('HRLI'))
      if (is_hrl .or. is_hrli) then
         do k=1,nkm1
            prest(:,k) = zpplus(:) * zsigt(:,k) 
         enddo
         ! Calculate HRL (liquid), including diag level nk
         if (is_hrl) then
            call mfohr4(zhrl, zhuplus, ztplus, prest, ni, nkm1, ni, .false.)
            call mfohr4(zhrl(:,nk), zqdiag, ztdiag, zpplus, ni, 1, ni ,.false.)
         endif
         ! Calculate HRI or HRL (ice or liquid), including diag level nk
         if (is_hrli) then
            call mfohr4(zhrli, zhuplus, ztplus, prest, ni, nkm1, ni, .true.)
            call mfohr4(zhrli(:,nk), zqdiag, ztdiag, zpplus, ni, 1, ni ,.true.)
         endif
      endif

      !****************************************************************
      !     AVERAGES
      !     --------

      w_reset = 1.
      if (lkount0 .or. lreset) w_reset = 0.
         
      if (laccum) then
         uvs = sqrt(zudiag*zudiag + zvdiag*zvdiag)
         
         DO_MAX(zhrsmax, zrhdiag, w_reset)
         DO_MIN(zhrsmin, zrhdiag, w_reset)
         DO_MAX(zuvsmax, uvs, w_reset)

         DO_MIN(zttmin, ztplus, w_reset)
         DO_MAX(zttmax, ztplus, w_reset)
         DO_MIN(ztadvmin, ztadv, w_reset)
         DO_MAX(ztadvmax, ztadv, w_reset)
      endif
      
      if (laccum .and. .not. lkount0) then

         DO_AVG(ztccm, ztcc,  moyhri, w_reset)
         DO_AVG(ztslm, ztcsl, moyhri, w_reset)
         DO_AVG(ztsmm, ztcsm, moyhri, w_reset)
         DO_AVG(ztshm, ztcsh, moyhri, w_reset)
         DO_AVG(ztzlm, ztczl, moyhri, w_reset)
         DO_AVG(ztzmm, ztczm, moyhri, w_reset)
         DO_AVG(ztzhm, ztczh, moyhri, w_reset)

         DO_AVG(zhusavg,   zqdiag, moyhri, w_reset)
         DO_AVG(ztdiagavg, ztdiag, moyhri, w_reset)
         DO_AVG(zp0avg,    zpplus, moyhri, w_reset)
         DO_AVG(zuvsavg,   uvs,    moyhri, w_reset)

         DO_AVG(zflwm,   zflw,   moyhri, w_reset)
         DO_AVG(zfshm,   zfsh,   moyhri, w_reset)
         DO_AVG(zfcmy,   zfc_ag, moyhri, w_reset)
         DO_AVG(zfvmm,   zfvap_ag, moyhri, w_reset)
         DO_AVG(zfvmy,   zfv_ag, moyhri, w_reset)
         DO_AVG(zkshalm, zkshal, moyhri, w_reset)
         DO_AVG(ziwvm,   ziwv  , moyhri, w_reset)
         DO_AVG(zq1appm, zq1app, moyhri, w_reset)
         DO_AVG(zq2appm, zq2app, moyhri, w_reset)
         DO_AVG(ztlwpm,  ztlwp , moyhri, w_reset)
         DO_AVG(zt2im,   zt2i  , moyhri, w_reset)
         DO_AVG(ztiim,   ztii  , moyhri, w_reset)
         DO_AVG(ztiwpm,  ztiwp , moyhri, w_reset)

         ! Compute conservation averages on request
         DO_AVG(zconedynm, zconedyn, moyhri, w_reset)
         DO_AVG(zconecndm, zconecnd, moyhri, w_reset)
         DO_AVG(zconedcm,  zconedc,  moyhri, w_reset)
         DO_AVG(zconemcm,  zconemc,  moyhri, w_reset)
         DO_AVG(zconemom,  zconemo,  moyhri, w_reset)
         DO_AVG(zconepblm, zconepbl, moyhri, w_reset)
         DO_AVG(zconephym, zconephy, moyhri, w_reset)
         DO_AVG(zconeradm, zconerad, moyhri, w_reset)
         DO_AVG(zconescm,  zconesc,  moyhri, w_reset)
         DO_AVG(zconetotm, zconetot, moyhri, w_reset)
         DO_AVG(zconqdynm, zconqdyn, moyhri, w_reset)
         DO_AVG(zconqcndm, zconqcnd, moyhri, w_reset)
         DO_AVG(zconqdcm,  zconqdc,  moyhri, w_reset)
         DO_AVG(zconqmcm,  zconqmc,  moyhri, w_reset)
         DO_AVG(zconqmom,  zconqmo,  moyhri, w_reset)
         DO_AVG(zconqpblm, zconqpbl, moyhri, w_reset)
         DO_AVG(zconqphym, zconqphy, moyhri, w_reset)
         DO_AVG(zconqradm, zconqrad, moyhri, w_reset)
         DO_AVG(zconqscm,  zconqsc,  moyhri, w_reset)
         DO_AVG(zconqtotm, zconqtot, moyhri, w_reset)

         DO_AVG(zo3tcm,    zo3tc,    moyhri, w_reset)
         DO_AVG(zo3ctcm,   zo3ctc,   moyhri, w_reset)
         DO_AVG(zo3avg,    zo3lplus, moyhri, w_reset)
         DO_AVG(zo3ccolm,  zo3ccol,  moyhri, w_reset)
         DO_AVG(zo3colm,   zo3col,   moyhri, w_reset)
         DO_AVG(zo1chmtdm, zo1chmtd, moyhri, w_reset)
         DO_AVG(zo4chmtdm, zo4chmtd, moyhri, w_reset)
         DO_AVG(zo6chmtdm, zo6chmtd, moyhri, w_reset)
         DO_AVG(zo7chmtdm, zo7chmtd, moyhri, w_reset)
         DO_AVG(zo3chmtdm, zo3chmtd, moyhri, w_reset)

         DO_AVG(zch4tcm, zch4tc, moyhri, w_reset)
         DO_AVG(zn2otcm, zn2otc, moyhri, w_reset)
         DO_AVG(zf11tcm, zf11tc, moyhri, w_reset)
         DO_AVG(zf12tcm, zf12tc, moyhri, w_reset)

         DO_AVG(zch4avg,  zch4lplus, moyhri, w_reset)
         DO_AVG(zn2oavg,  zn2olplus, moyhri, w_reset)
         DO_AVG(zf11avg,  zf11lplus, moyhri, w_reset)
         DO_AVG(zf12avg,  zf12lplus, moyhri, w_reset)
         
         DO_AVG(zch4colm, zch4col, moyhri, w_reset)
         DO_AVG(zn2ocolm, zn2ocol, moyhri, w_reset)
         DO_AVG(zf11colm, zf11col, moyhri, w_reset)
         DO_AVG(zf12colm, zf12col, moyhri, w_reset)
         
         DO_AVG(zch4chmtdm, zch4chmtd, moyhri, w_reset)
         DO_AVG(zn2ochmtdm, zn2ochmtd, moyhri, w_reset)
         DO_AVG(zf11chmtdm, zf11chmtd, moyhri, w_reset)
         DO_AVG(zf12chmtdm, zf12chmtd, moyhri, w_reset)
         
         DO_AVG(zccnm, zftot, moyhri, w_reset)
         DO_AVG(ztim,  zti,   moyhri, w_reset)
         DO_AVG(zt2m,  zt2,   moyhri, w_reset)
         
         DO_AVG(zudifvm, zudifv, moyhri, w_reset)
         DO_AVG(zvdifvm, zvdifv, moyhri, w_reset)
         DO_AVG(ztdifvm, ztdifv, moyhri, w_reset)
         DO_AVG(zqdifvm, zqdifv, moyhri, w_reset)
         
         DO_AVG(ztadvm, ztadv, moyhri, w_reset)
         DO_AVG(zuadvm, zuadv, moyhri, w_reset)
         DO_AVG(zvadvm, zvadv, moyhri, w_reset)
         DO_AVG(zqadvm, zqadv, moyhri, w_reset)
         
         DO_AVG(zqmetoxm, zqmetox, moyhri, w_reset)
         
         DO_AVG(zhushalm, zhushal, moyhri, w_reset)
         DO_AVG(ztshalm,  ztshal,  moyhri, w_reset)
         
         DO_AVG(zlwcm,    zlwc,    moyhri, w_reset)
         DO_AVG(ziwcm,    ziwc,    moyhri, w_reset)
         DO_AVG(zlwcradm, zlwcrad, moyhri, w_reset)
         DO_AVG(ziwcradm, ziwcrad, moyhri, w_reset)
         DO_AVG(zcldradm, zcldrad, moyhri, w_reset)
         
         DO_AVG(ztphytdm,  ztphytd,  moyhri, w_reset)
         DO_AVG(zhuphytdm, zhuphytd, moyhri, w_reset)
         DO_AVG(zuphytdm,  zuphytd,  moyhri, w_reset)
         DO_AVG(zvphytdm,  zvphytd,  moyhri, w_reset)
         
         !#deep
         DO_AVG(zzctem,  zcte,  moyhri, w_reset)
         DO_AVG(zzstem,  zste,  moyhri, w_reset)
         DO_AVG(zzcqem,  zcqe,  moyhri, w_reset)
         DO_AVG(zzsqem,  zsqe,  moyhri, w_reset)
         DO_AVG(zzcqcem, zcqce, moyhri, w_reset)
         DO_AVG(zzsqcem, zsqce, moyhri, w_reset)
         
         DO_AVG(zutofdm, zutofd, moyhri, w_reset)
         DO_AVG(zvtofdm, zvtofd, moyhri, w_reset)
         DO_AVG(zttofdm, zttofd, moyhri, w_reset)
         
         DO_AVG(ztgwdm, ztgwd, moyhri, w_reset)
         DO_AVG(zugwdm, zugwd, moyhri, w_reset)
         DO_AVG(zvgwdm, zvgwd, moyhri, w_reset)
         DO_AVG(zugnom, zugno, moyhri, w_reset)
         DO_AVG(zvgnom, zvgno, moyhri, w_reset)
         DO_AVG(ztgnom, ztgno, moyhri, w_reset)
         
         DO_AVG(zmqem, zmqe, moyhri, w_reset)
         DO_AVG(zmtem, zmte, moyhri, w_reset)
         DO_AVG(zumim, zumid, moyhri, w_reset)
         DO_AVG(zvmim, zvmid, moyhri, w_reset)
         
         DO_AVG(ztfcpm,   ztfcp,   moyhri, w_reset)
         DO_AVG(zhufcpm,  zhufcp,  moyhri, w_reset)
         DO_AVG(zqckfcm,  zqckfc,  moyhri, w_reset)
         DO_AVG(zumfkfcm, zumfkfc, moyhri, w_reset)
         DO_AVG(zdmfkfcm, zdmfkfc, moyhri, w_reset)
         DO_AVG(zudcm,    zufcp,   moyhri, w_reset)
         DO_AVG(zvdcm,    zvfcp,   moyhri, w_reset)
         
         DO_AVG(zufcp1m, zufcp1, moyhri, w_reset)
         DO_AVG(zufcp2m, zufcp2, moyhri, w_reset)
         DO_AVG(zufcp3m, zufcp3, moyhri, w_reset)
         DO_AVG(zsufcpm, zsufcp, moyhri, w_reset)
         DO_AVG(zvfcp1m, zvfcp1, moyhri, w_reset)
         DO_AVG(zvfcp2m, zvfcp2, moyhri, w_reset)
         DO_AVG(zvfcp3m, zvfcp3, moyhri, w_reset)
         DO_AVG(zsvfcpm, zsvfcp, moyhri, w_reset)

         DO_AVG(zabekfcm,  zabekfc,   moyhri, w_reset)
         DO_AVG(zcapekfcm, zcapekfc,  moyhri, w_reset)
         DO_AVG(zcinkfcm,  zcinkfc,   moyhri, w_reset)
         DO_AVG(zwumkfcm,  zwumaxkfc, moyhri, w_reset)
         DO_AVG(zzbaskfcm, zzbasekfc, moyhri, w_reset)
         DO_AVG(zztopkfcm, zztopkfc,  moyhri, w_reset)
         DO_AVG(zkkfcm,    zkkfc,     moyhri, w_reset)
         DO_AVG(zkmidm,    zkmid,     moyhri, w_reset)

         DO_AVG(zsumf, zustress, moyhri, w_reset)
         DO_AVG(zsvmf, zvstress, moyhri, w_reset)
         DO_AVG(zfqm,  zfq,      moyhri, w_reset)

         DO_AVG(zuscm, ztusc, moyhri, w_reset)
         DO_AVG(zvscm, ztvsc, moyhri, w_reset)
      endif

      !****************************************************************
      !     ACCUMULATORS

      !Set accumulators to zero at the beginning and after every acchr hours,
      !and by default (acchr=0) as the model step goes through 0.
      w_reset = 1.
      if (lkount0 .or. lacchr) w_reset = 0.

      if (.not.lkount0) then
         
         DO_ACC(zals_rn1, ztls_rn1, dt, w_reset)
         DO_ACC(zals_rn2, ztls_rn2, dt, w_reset)
         DO_ACC(zals_fr1, ztls_fr1, dt, w_reset)
         DO_ACC(zals_fr2, ztls_fr2, dt, w_reset)
         DO_ACC(zass_sn1, ztss_sn1, dt, w_reset)
         DO_ACC(zass_sn2, ztss_sn2, dt, w_reset)
         if (stcond == 'MP_P3') then
            DO_ACC(zass_ws, ztss_ws, dt, w_reset)
         endif
         DO_ACC(zass_sn3, ztss_sn3, dt, w_reset)
         DO_ACC(zass_pe1, ztss_pe1, dt, w_reset)
         DO_ACC(zass_pe2, ztss_pe2, dt, w_reset)
         DO_ACC(zass_pe2l, ztss_pe2l, dt, w_reset)
         DO_ACC(zass_snd, ztss_snd, dt, w_reset)
         DO_ACC(zass_mx,  ztss_mx,  dt, w_reset)
         if (associated(zass_s2l)) &
              zass_s2l = zass_snd / max(zass_sn1, zass_sn2, zass_sn3, 1.e-18)

         DO_ACC(zasc,  ztsc,  dt, w_reset)
         DO_ACC(zascs, ztscs, dt, w_reset)
         DO_ACC(zalc,  ztlc , dt, w_reset)
         DO_ACC(zalcm, ztlcm, dt, w_reset)
         DO_ACC(zalcs, ztlcs, dt, w_reset)
         DO_ACC(zass,  ztss,  dt, w_reset)
         DO_ACC(zals,  ztls,  dt, w_reset)
         zpc  = zalc + zasc + zalcs + zascs + zalcm
         zpy  = zalc + zasc
         zpz  = zalcs + zascs
         zacm = zalcm
         zae  = zals + zass
         zpr  = zpc + zae
         
         !# Accumulation of precipitation (in m)
         DO_ACC(zrainaf, zrainrate, dt, w_reset)
         DO_ACC(zsnowaf, zsnowrate, dt, w_reset)

         DO_ACC(zeiaf, zei,   dt, w_reset)
         DO_ACC(zevaf, zev,   dt, w_reset)
         DO_ACC(zfiaf, zfdsi, dt, w_reset)
         DO_ACC(zfsaf, zfdss, dt, w_reset)
         DO_ACC(zivaf, ziv,   dt, w_reset)
         DO_ACC(zntaf, znt,   dt, w_reset)
         DO_ACC(zflusolaf, zflusolis, dt, w_reset)

         !# Accumu of sfc and toa net clear sky fluxes, avail with cccmarad
         DO_ACC(zclbaf,  zclb,  dt, w_reset)
         DO_ACC(zcltaf,  zclt,  dt, w_reset)
         DO_ACC(zcstaf,  zcstt, dt, w_reset)
         DO_ACC(zcsbaf,  zcsb,  dt, w_reset)
         DO_ACC(zfsdaf,  zfsd,  dt, w_reset)
         DO_ACC(zfsfaf,  zfsf,  dt, w_reset)
         DO_ACC(zfsiaf,  zfsi,  dt, w_reset)
         DO_ACC(zfsvaf,  zfsv,  dt, w_reset)
         DO_ACC(zparraf, zparr, dt, w_reset)
         DO_ACC(zsiaf ,  zfnsi, dt, w_reset)
         DO_ACC(zfnsaf,  zfns,  dt, w_reset)

         DO_ACC(zflaf, zfl,     dt, w_reset)
         DO_ACC(zfcaf, zfc_ag,  dt, w_reset)
         DO_ACC(zfvaf, zfv_ag,  dt, w_reset)
         
         !# Accumulation of lightning threat (in number of flashes/m2)
         DO_ACC(zafoudre, zfoudre, dt, w_reset)
      endif
      
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
