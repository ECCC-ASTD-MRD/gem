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
   subroutine calcdiag1(dbus, fbus, vbus, dt, kount, ni, nk)
      !@Object Calculates averages and accumulators of tendencies and diagnostics
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use tdpack_const, only: CHLC, CPD, GRAV, TCDK, RAUW, CAPPA
      use integrals, only: int_profile, INT_OK
      use pbl_height, only: pbl_height1
      use sfclayer, only: sl_prelim, sl_sfclayer, SL_OK
      use phy_options
      use phybus
      use phybudget, only: pb_compute, pb_residual
      use phy_status, only: PHY_OK
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      !          - Input/Output -
      ! dbus     dynamic bus
      ! fbus     permanent bus
      !          - input -
      ! vbus     volatile (output) bus
      !          - input -
      ! dsiz     dimension of d
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! kount    timestep number
      ! dt       length of timestep
      ! n        horizontal running length
      ! nk       vertical dimension

      integer, intent(in) :: kount, ni, nk
      real, intent(in) :: dt
      real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)

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

      logical :: lmoyhr, laccum, lreset, lavg, lkount0, lacchr
      integer :: i, k, moyhr_steps, istat, istat1, nkm1
      real :: moyhri, tempo, tempo2, sol_stra, sol_conv, liq_stra, liq_conv, sol_mid, liq_mid
      real, dimension(ni) :: uvs, vmod, vdir, th_air, hblendm, ublend, &
           vblend, z0m_ec, z0t_ec, esdiagec, tsurfec, qsurfec
      real(REAL64), dimension(ni) :: en0, pw0, en1, pw1
      real, dimension(ni,nk) :: presinv, t2inv, tiinv, lwcp, iwcp
      real, dimension(ni,nk-1) :: q_grpl, iiwc, prest

#define PHYPTRDCL
#include "calcdiag_ptr.hf"

      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'calcdiag [BEGIN]')
      if (timings_L) call timing_start_omp(470, 'calcdiag', 46)

      nkm1 = nk-1

#undef PHYPTRDCL
#include "calcdiag_ptr.hf"

      call init2nan(uvs, vmod, vdir, th_air, hblendm, ublend)
      call init2nan(vblend, z0m_ec, z0t_ec, esdiagec, tsurfec, qsurfec)
      call init2nan(en0, pw0, en1, pw1)
      call init2nan(presinv, t2inv, tiinv, q_grpl, iiwc, lwcp, iwcp)

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
         moyhri = 1./float(moyhr_steps)
      endif

      !# Final PBL height
      istat = pbl_height1(zh,ztplus,zhuplus,zuplus,zvplus,zgzmom,zgztherm,zsigt, &
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
      if (.not.(lkount0 .and. fluvert /= 'NIL').or.fluvert=='SURFACE') then

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
         esdiagec(:) = max(esdiagec(:), 0.)
         ztddiagec(:) = ztdiagec(:) - esdiagec(:)

      endif ECMWF_SCREEN
      
      !****************************************************************
      !     PRECIPITATION RATES AND ACCUMULATIONS
      !     -------------------------------------

      !# Set precipitation accumulators to zero at the beginning and after every
      !  acchr hours, and by default (acchr=0) as the model step goes through 0.
      IF_RESET_PRECIP: if (lkount0 .or. lacchr) then

         zasc(:)  = 0.
         zascs(:) = 0.
         zalc(:)  = 0.
         zalcm(:) = 0.
         zalcs(:) = 0.
         zass(:)  = 0.
         zals(:)  = 0.
         zpc(:)   = 0.
         zpy(:)   = 0.
         zpz(:)   = 0.
         zae(:)   = 0.
         zpr(:)   = 0.
         zazr(:)  = 0.
         zrn(:)   = 0.
         zaip(:)  = 0.
         zsn(:)   = 0.

         if (associated(zals_rn1)) zals_rn1(:)  = 0.
         if (associated(zals_rn2)) zals_rn2(:)  = 0.
         if (associated(zals_fr1)) zals_fr1(:)  = 0.
         if (associated(zals_fr2)) zals_fr2(:)  = 0.
         if (associated(zass_sn1)) zass_sn1(:)  = 0.
         if (associated(zass_sn2)) zass_sn2(:)  = 0.
         if (associated(zass_sn3)) zass_sn3(:)  = 0.
         if (associated(zass_pe1)) zass_pe1(:)  = 0.
         if (associated(zass_pe2)) zass_pe2(:)  = 0.
         if (associated(zass_pe2l)) zass_pe2l(:) = 0.
         if (associated(zass_snd)) zass_snd(:)  = 0.
         if (associated(zass_mx)) zass_mx(:)   = 0.
         if (associated(zass_s2l)) zass_s2l(:)  = 0.

      endif IF_RESET_PRECIP

      IF_KOUNT_NOT_0: if (.not.lkount0) then

         ! Diagnostics on precipitation type.
         if (pcptype == 'BOURGE' .or. &
              (stcond(1:3) == 'MP_' .and. pcptype == 'NIL')) then
            !Note: 'bourge2' is always called if microphysics scheme is used and pcptyp=nil,
            !       but it is only applied to implicit (convection) precipitation
            call bourge2(zfneige1d, zfip1d, ztplus, zsigw, zpplus, ni, nk-1)
         else if (pcptype == 'BOURGE3D') then
            call bourge1_3d(zfneige2d, zfip2d, ztplus, zsigw, zpplus, ni, nk-1)
            zfneige1d(1:ni) => zfneige2d(:,nk)
            zfip1d(1:ni)    => zfip2d(:,nk)
         endif

         IF_BOURG3D: if (pcptype == 'BOURGE3D' .and. stcond == 'CONSUN') then
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
                  if (ztplus(i,k) < TCDK) then
                     zzr3d(i,k) = (1.-zfip2d(i,k))*tempo*dt
                  else
                     zzr3d(i,k) = 0.
                  endif
                  zazr3d(i,k) = zazr3d(i,k) + zzr3d(i,k)
                  zpl3d(i,k) = zfip2d(i,k)*tempo*dt
                  zaip3d(i,k) = zaip3d(i,k) + zpl3d(i,k)
               enddo
            enddo
         endif IF_BOURG3D

!VDIR NODEP
         DO_NI: do i = 1,ni

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
            IF_MP: if (stcond(1:3) == 'MP_')  then

               !diagnose freezing drizzle/rain based on sfc temperature
               !---override rn/fr partition from my2 [to be removed]
               if (stcond(1:6) == 'MP_MY2') then
                  ztls_rn1(i) = ztls_rn1(i) + ztls_fr1(i)
                  ztls_rn2(i) = ztls_rn2(i) + ztls_fr2(i)
               endif

               if  (ztplus(i,nk) < TCDK) then
                  ztls_fr1(i) = ztls_rn1(i)
                  ztls_fr2(i) = ztls_rn2(i)
                  ztls_rn1(i) = 0.
                  ztls_rn2(i) = 0.
               endif

               !tls:  rate of liquid precipitation (sum of rain and drizzle)
               ztls(i) = ztls_rn1(i) + ztls_rn2(i) + ztls_fr1(i) + ztls_fr2(i)

               !tss:  rate of solid precipitation (sum of all fozen precipitation types)
               ztss(i) = ztss_sn1(i) + ztss_sn2(i) +ztss_sn3(i) + ztss_pe1(i) + ztss_pe2(i)

               !tss_mx:  rate of mixed precipitation (liquid and solid simultaneously [>0.01 mm/h each])
               if (ztls(i) > 2.78e-9 .and. ztss(i) > 2.78e-9) then   ![note: 2.78e-9 m/s = 0.01 mm/h]
                  ztss_mx(i) = ztls(i) + ztss(i)
               else
                  ztss_mx(i) = 0.
               endif

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

            if (stcond == 'CONSUN') then
               tempo    = ztls(i) + ztss(i)
               sol_stra = max(0., zfneige1d(i)*tempo)
               liq_stra = max(0., tempo-sol_stra)
               tempo    = ztlc(i) + ztlcm(i) + ztlcs(i) + ztsc(i) + ztscs(i)
               sol_conv = max(0., zfneige1d(i)*tempo)
               liq_conv = max(0., tempo-sol_conv)
            elseif (any(pcptype == (/ &
                 'NIL    ', &
                 'SPS_W19', &
                 'SPS_FRC'/))) then
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

            IF_MP_EXPL: if (stcond(1:3)=='MP_' .and. pcptype(1:3)=='NIL') then

               tempo    = ztlc(i) + ztlcm(i) + ztlcs(i) + ztsc(i) + ztscs(i)   !total convective (implicit)
               sol_conv = max(0., zfneige1d(i)*tempo)
               liq_conv = max(0., tempo-sol_conv)
               tempo2   = zazr(i)
               !add explicit + diagnostic portion of rain/freezing rain from convective schemes:
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
               !add explicit + diagnostic portion of ice pellets from convective schemes:
               zaip(i) = zaip(i)                                     &
                    + ztss_pe1(i)*dt                            &  !from microphysics
                    + zfip1d(i)*tempo*dt                           !from convective schemes
               !note: Hail from M-Y (tss_pe2) is not included in total ice pellets (aip)
               !add explicit + diagnostic portion of snow from convective schemes:
               zsn(i)  = zsn(i)                                      &
                    + (ztss_sn1(i)+ztss_sn2(i)+ztss_sn3(i))*dt  &  !from microphysics
                    + sol_conv*dt                                  !from convective schemes
               !note: contribution to SN from microphysics is from ice, snow,
               !      and graupel, instantaneous solid-to-liquid ratio for
               !      pcp rate total snow (i+s+g):
               tempo       =  max(ztss_sn1(i)+ztss_sn2(i)+ztss_sn3(i), 1.e-18)
               ztss_s2l(i) = ztss_snd(i)/tempo

            else

               !diagnostic partitioning of rain vs. freezing rain:
               tempo =  liq_stra + liq_conv
               if  (ztplus(i,nk) < TCDK) then
                  zzr(i) = (1.-zfip1d(i))*tempo*dt
                  zazr(i) = zazr(i) + zzr(i)
                  if (tempo > 0.) zzrflag(i) = 1.
               else
                  zrn(i) = zrn(i) + (1.-zfip1d(i))*tempo*dt
                  if (pcptype == 'BOURGE3D' .and. tempo > 0.) zzrflag(i) = 1.
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
            q_grpl(:,:) = a_qi_4(:,:) + a_qi_5(:,:)
            iiwc(:,:)   = a_qi_1(:,:) + a_qi_2(:,:) + a_qi_3(:,:) +     &
                 a_qi_4(:,:) + a_qi_5(:,:) + a_qi_6(:,:)
            call lightning2(zfoudre,zp0_plus,zsigm,ztplus,zwplus,q_grpl,iiwc,ni,nk)
         elseif (stcond(1:6)=='MP_MY2') then
            q_grpl(:,:) = zqgplus(:,:)
            iiwc(:,:)   = zqiplus(:,:) + zqnplus(:,:) + zqgplus(:,:)
            call lightning2(zfoudre,zp0_plus,zsigm,ztplus,zwplus,q_grpl,iiwc,ni,nk)
         endif
      endif

      !****************************************************************
      !     Energy budget diagnostics
      !     ---------------------------------

      ! Compute Q1, Q2 apparent heat source / moisture sink on request
      Q1_BUDGET: if (ebdiag) then
         do k=1,nk
            presinv(:,nk-(k-1)) = zpplus*zsigt(:,k)
            t2inv(:,nk-(k-1)) = zt2(:,k)
            tiinv(:,nk-(k-1)) = zti(:,k)
         enddo
         istat  = int_profile(zt2i,t2inv,presinv,zpplus*zsigt(:,1),zpplus*zsigt(:,nk))
         istat1 = int_profile(ztii,tiinv,presinv,zpplus*zsigt(:,1),zpplus*zsigt(:,nk))
         if (istat /= INT_OK .or. istat1 /= INT_OK) then
            call physeterror('calcdiag', 'Problem in radiative tendency integrals')
            return
         endif
         zq1app = CPD*(zt2i+ztii)/GRAV + CHLC*zrt*RAUW + zfc_ag
         zq2app = CHLC*(zrt*RAUW - zflw)
      endif Q1_BUDGET
      
      ! Compute physics budget residual
      en0(:) = 0.
      if (associated(zcone0)) en0(:) = dble(zcone0(:))
      pw0(:) = 0.
      if (associated(zconq0)) pw0(:) = dble(zconq0(:))
      if (pb_residual(zconephy, zconqphy, en0, pw0, dbus, fbus, vbus, &
           delt, nkm1, F_rain=zrt*RAUW, F_shf=zfc_ag, F_wvf=zflw, &
           F_rad=znetrad) /= PHY_OK) then
         call physeterror('calcdiag', &
              'Problem computing physics budget residual')
         return
      endif

      ! Compute full-model budget residual
      en1(:) = 0.
      if (associated(zcone1)) en1(:) = dble(zcone1(:))
      pw1(:) = 0.
      if (associated(zconq1)) pw1(:) = dble(zconq1(:))
      if (pb_residual(zconetot, zconqtot, en1, pw1, dbus, fbus, vbus, &
           delt, nkm1, F_rain=zrt*RAUW, F_shf=zfc_ag, F_wvf=zflw, &
           F_rad=znetrad) /= PHY_OK) then
         call physeterror('calcdiag', &
              'Problem computing full-model budget residual')
         return
      endif
      
      ! Compute post-physics budget state
      if (pb_compute(zcone1, zconq1, en1, pw1, dbus, fbus, vbus, nkm1) /= PHY_OK) then
         call physeterror('phystepinit', &
              'Problem computing post-physics budget state')
         return
      endif
      if (associated(zcone1)) zcone1(:) = real(en1(:))
      if (associated(zconq1)) zconq1(:) = real(pw1(:))

      ! Compute surface energy budget (
      zfl(:)   = zfns(:) - zfv_ag(:) - zfc_ag(:)

      !****************************************************************
      !     Moisture diagnostics
      !     --------------------
      do k=1,nkm1
         prest(:,k) = zpplus(:) * zsigt(:,k) 
      enddo
      ! Calculate HRL (liquid), including diag level nk
      call mfohr4(zhrl, zhuplus, ztplus, prest, ni, nkm1, ni, .false.)
      call mfohr4(zhrl(:,nk), zqdiag, ztdiag, zpplus, ni, 1, ni ,.false.)
      ! Calculate HRI or HRL (ice or liquid), including diag level nk
      call mfohr4(zhrli, zhuplus, ztplus, prest, ni, nkm1, ni, .true.)
      call mfohr4(zhrli(:,nk), zqdiag, ztdiag, zpplus, ni, 1, ni ,.true.)

      !****************************************************************
      !     AVERAGES
      !     --------

      !set averages to zero every moyhr hours
      IF_ACCUM_0: if (laccum) then

         !pre-calculate screen wind modulus
         uvs(1:ni) = sqrt(zudiag(1:ni)*zudiag(1:ni) + zvdiag(1:ni)*zvdiag(1:ni))

         IF_RESET_0: if (lkount0 .or. lreset) then

            zflwm(:)   = 0.
            zfshm(:)   = 0.
            zfcmy(:)   = 0.0
            zfvmm(:)   = 0.0
            zfvmy(:)   = 0.0
            zkshalm(:) = 0.0
            ziwvm(:)   = 0.0
            zq1appm(:) = 0.0
            zq2appm(:) = 0.0
            ztlwpm(:)  = 0.0
            zt2im(:)   = 0.0
            ztiim(:)   = 0.0
            ztiwpm(:)  = 0.0
            if (etccdiag) then
               ztccm(:)  = 0.0
               ztslm(:)  = 0.0
               ztsmm(:)  = 0.0
               ztshm(:)  = 0.0
               ztzlm(:)  = 0.0
               ztzmm(:)  = 0.0
               ztzhm(:)  = 0.0
            endif
            zhrsmax(:) = zrhdiag(:)
            zhrsmin(:) = zrhdiag(:)
            zhusavg(:) = 0.0
            ztdiagavg(:)= 0.0
            zp0avg(:)  = 0.0
            zuvsavg(:) = 0.0
            zuvsmax(:) = uvs(:)

            !minimum and maximum temperature and temp. tendencies
            zttmin(:,:) = ztplus(:,:)
            zttmax(:,:) = ztplus(:,:)
            ztadvmin(:,:) = ztadv(:,:)
            ztadvmax(:,:) = ztadv(:,:)

            if (llinoz .and. out_linoz) then
               ! 2D Ozone
               zo3tcm  (:) = 0.0
               zo3ctcm (:) = 0.0
               ! 3D Ozone
               zo3avg  (:,:) = 0.0
               zo3ccolm(:,:) = 0.0
               zo3colm (:,:) = 0.0
               zo1chmtdm(:,:) = 0.0
               zo4chmtdm(:,:) = 0.0
               zo6chmtdm(:,:) = 0.0
               zo7chmtdm(:,:) = 0.0
               zo3chmtdm(:,:) = 0.0
            end if

            if (llingh .and. out_linoz) then
               ! 2D GHG
               zch4tcm (:) = 0.0
               zn2otcm (:) = 0.0
               zf11tcm (:) = 0.0
               zf12tcm (:) = 0.0
               ! 3D GHG
               zch4avg (:,:) = 0.0
               zn2oavg (:,:) = 0.0
               zf11avg (:,:) = 0.0
               zf12avg (:,:) = 0.0
               zch4colm(:,:) = 0.0
               zn2ocolm(:,:) = 0.0
               zf11colm(:,:) = 0.0
               zf12colm(:,:) = 0.0
               zch4chmtdm(:,:) = 0.0
               zn2ochmtdm(:,:) = 0.0
               zf11chmtdm(:,:) = 0.0
               zf12chmtdm(:,:) = 0.0
            end if

            zccnm(:,1:nk-1)  = 0.0
            ztim(:,1:nk-1)   = 0.0
            zt2m(:,1:nk-1)   = 0.0
            ztgwdm(:,1:nk-1) = 0.0
            zutofdm(:,1:nk-1) = 0.0
            zvtofdm(:,1:nk-1) = 0.0
            zttofdm(:,1:nk-1) = 0.0
            zugwdm(:,1:nk-1) = 0.0
            zvgwdm(:,1:nk-1) = 0.0
            zugnom(:,1:nk-1) = 0.0
            zvgnom(:,1:nk-1) = 0.0
            ztgnom(:,1:nk-1) = 0.0
            zudifvm(:,1:nk-1) = 0.0
            zvdifvm(:,1:nk-1) = 0.0
            ztdifvm(:,1:nk-1) = 0.0
            zqdifvm(:,1:nk-1) = 0.0
            ztadvm(:,1:nk-1)  = 0.0
            zuadvm(:,1:nk-1)  = 0.0
            zvadvm(:,1:nk-1)  = 0.0
            zqadvm(:,1:nk-1)  = 0.0
            zqmetoxm(:,1:nk-1) = 0.0
            zhushalm(:,1:nk-1) = 0.0
            ztshalm (:,1:nk-1) = 0.0
            zlwcm(:,1:nk-1)    = 0.0
            ziwcm(:,1:nk-1)    = 0.0
            zlwcradm(:,1:nk-1) = 0.0
            ziwcradm(:,1:nk-1) = 0.0
            zcldradm(:,1:nk-1) = 0.0
            ztphytdm(:,1:nk-1) = 0.0
            zhuphytdm(:,1:nk-1)= 0.0
            zuphytdm(:,1:nk-1) = 0.0
            zvphytdm(:,1:nk-1) = 0.0
            zzctem(:,1:nk-1)  = 0.0
            zzstem(:,1:nk-1)  = 0.0
            zzcqem(:,1:nk-1)  = 0.0
            zzcqcem(:,1:nk-1) = 0.0
            zzsqem(:,1:nk-1)  = 0.0
            zzsqcem(:,1:nk-1) = 0.0
            if (associated(zmqem)) zmqem(:,1:nk-1) = 0.0
            if (associated(zmtem)) zmtem(:,1:nk-1) = 0.0
            if (associated(zumim)) zumim(:,1:nk-1) = 0.0
            if (associated(zvmim)) zvmim(:,1:nk-1) = 0.0
            if (associated(zkmidm)) zkmidm(:) = 0.0

            !#Note: all convec str in the if below must have same len
            if (any(convec == (/ &
                 'KFC     ', &
                 'KFC2    ', &
                 'BECHTOLD'  &
                 /))) then
               zabekfcm(:)  = 0.0
               zcapekfcm(:) = 0.0
               zcinkfcm(:)  = 0.0
               zwumkfcm(:)  = 0.0
               zzbaskfcm(:) = 0.0
               zztopkfcm(:) = 0.0
               zkkfcm(:)    = 0.0
               ztfcpm(:,:)   = 0.0
               zhufcpm(:,:)  = 0.0
               zqckfcm(:,:)  = 0.0
               zumfkfcm(:,:) = 0.0
               zdmfkfcm(:,:) = 0.0
               zudcm(:,:)   = 0.0
               zvdcm(:,:)   = 0.0

               if (cmt_comp_diag) then
                  zufcp1m(:,:) = 0.0
                  zufcp2m(:,:) = 0.0
                  zufcp3m(:,:) = 0.0
                  zsufcpm(:,:) = 0.0
                  zvfcp1m(:,:) = 0.0
                  zvfcp2m(:,:) = 0.0
                  zvfcp3m(:,:) = 0.0
                  zsvfcpm(:,:) = 0.0
               endif
            
               if (associated(zuscm)) zuscm(:,:) = 0.0
               if (associated(zvscm)) zvscm(:,:) = 0.0
            endif

            if (fluvert /= 'NIL') then
               zsumf(:) = 0.
               zsvmf(:) = 0.
               zfqm(:)  = 0.
            endif

            ! Reset conservation averages
            if (associated(zconecndm)) then
               zconedynm(:) = 0.
               zconecndm(:) = 0.
               zconedcm(:) = 0.
               zconemcm(:) = 0.
               zconemom(:) = 0.
               zconepblm(:) = 0.
               zconephym(:) = 0.
               zconeradm(:) = 0.
               zconescm(:) = 0.
               zconetotm(:) = 0.
               zconqdynm(:) = 0.
               zconqcndm(:) = 0.
               zconqdcm(:) = 0.
               zconqmcm(:) = 0.
               zconqmom(:) = 0.
               zconqpblm(:) = 0.
               zconqphym(:) = 0.
               zconqradm(:) = 0.
               zconqscm(:) = 0.
               zconqtotm(:) = 0.
            endif

         endif IF_RESET_0

      endif IF_ACCUM_0


      IF_ACCUM_1: if (laccum .and. .not.lkount0) then

         do i = 1, ni
            zflwm    (i) = zflwm   (i) + zflw(i) 
            zfshm    (i) = zfshm   (i) + zfsh(i) 
            zfcmy    (i) = zfcmy   (i) + zfc_ag(i)
            zfvmm    (i) = zfvmm   (i) + zfvap_ag(i)
            zfvmy    (i) = zfvmy   (i) + zfv_ag(i)
            zkshalm  (i) = zkshalm (i) + zkshal(i)
            ziwvm    (i) = ziwvm   (i) + ziwv  (i)
            zq1appm  (i) = zq1appm (i) + zq1app(i)
            zq2appm  (i) = zq2appm (i) + zq2app(i)
            ztlwpm   (i) = ztlwpm  (i) + ztlwp (i)
            zt2im    (i) = zt2im   (i) + zt2i  (i)
            ztiim    (i) = ztiim   (i) + ztii  (i)
            ztiwpm   (i) = ztiwpm  (i) + ztiwp (i)
            if (etccdiag) then
               ztccm(i) = ztccm(i) + ztcc(i)
               ztslm(i) = ztslm(i) + ztcsl(i)
               ztsmm(i) = ztsmm(i) + ztcsm(i)
               ztshm(i) = ztshm(i) + ztcsh(i)
               ztzlm(i) = ztzlm(i) + ztczl(i)
               ztzmm(i) = ztzmm(i) + ztczm(i)
               ztzhm(i) = ztzhm(i) + ztczh(i)
            endif
            zhrsmax  (i) = max(zhrsmax  (i) , zrhdiag (i))
            zhrsmin  (i) = min(zhrsmin  (i) , zrhdiag (i))
            zhusavg  (i) =     zhusavg  (i) + zqdiag  (i)
            ztdiagavg(i) =     ztdiagavg(i) + ztdiag  (i)
            zp0avg   (i) =     zp0avg   (i) + zpplus  (i)
            zuvsavg  (i) =     zuvsavg  (i) + uvs     (i)
            zuvsmax  (i) = max(zuvsmax  (i) , uvs     (i))
            if (lavg) then
               zflwm    (i) = zflwm    (i) * moyhri
               zfshm    (i) = zfshm    (i) * moyhri
               zfcmy    (i) = zfcmy    (i) * moyhri
               zfvmm    (i) = zfvmm    (i) * moyhri
               zfvmy    (i) = zfvmy    (i) * moyhri
               zkshalm  (i) = zkshalm  (i) * moyhri
               ziwvm    (i) = ziwvm    (i) * moyhri
               zq1appm  (i) = zq1appm  (i) * moyhri
               zq2appm  (i) = zq2appm  (i) * moyhri
               ztlwpm   (i) = ztlwpm   (i) * moyhri
               zt2im    (i) = zt2im    (i) * moyhri
               ztiim    (i) = ztiim    (i) * moyhri
               ztiwpm   (i) = ztiwpm   (i) * moyhri
               if (etccdiag) then
                  ztccm(i) = ztccm(i) * moyhri
                  ztslm(i) = ztslm(i) * moyhri
                  ztsmm(i) = ztsmm(i) * moyhri
                  ztshm(i) = ztshm(i) * moyhri
                  ztzlm(i) = ztzlm(i) * moyhri
                  ztzmm(i) = ztzmm(i) * moyhri
                  ztzhm(i) = ztzhm(i) * moyhri
               endif
               zhusavg  (i) = zhusavg  (i) * moyhri
               ztdiagavg(i) = ztdiagavg(i) * moyhri
               zp0avg   (i) = zp0avg   (i) * moyhri
               zuvsavg  (i) = zuvsavg  (i) * moyhri
            endif
         end do

         ! Compute conservation averages on request
         CONSERVE_RESIDUAL_AVERAGING: if (associated(zconecndm)) then
            zconedynm(:) = zconedynm(:) + zconedyn(:)
            zconecndm(:) = zconecndm(:) + zconecnd(:)
            zconedcm(:) = zconedcm(:) + zconedc(:)
            zconemcm(:) = zconemcm(:) + zconemc(:)
            zconemom(:) = zconemom(:) + zconemo(:)
            zconepblm(:) = zconepblm(:) + zconepbl(:)
            zconephym(:) = zconephym(:) + zconephy(:)
            zconeradm(:) = zconeradm(:) + zconerad(:)
            zconescm(:) = zconescm(:) + zconesc(:)
            zconetotm(:) = zconetotm(:) + zconetot(:)
            zconqdynm(:) = zconqdynm(:) + zconqdyn(:)
            zconqcndm(:) = zconqcndm(:) + zconqcnd(:)
            zconqdcm(:) = zconqdcm(:) + zconqdc(:)
            zconqmcm(:) = zconqmcm(:) + zconqmc(:)
            zconqmom(:) = zconqmom(:) + zconqmo(:)
            zconqpblm(:) = zconqpblm(:) + zconqpbl(:)
            zconqphym(:) = zconqphym(:) + zconqphy(:)
            zconqradm(:) = zconqradm(:) + zconqrad(:)
            zconqscm(:) = zconqscm(:) + zconqsc(:)
            zconqtotm(:) = zconqtotm(:) + zconqtot(:)
            if (lavg) then
               zconedynm(:) = zconedynm(:) * moyhri
               zconecndm(:) = zconecndm(:) * moyhri
               zconedcm(:) = zconedcm(:) * moyhri
               zconemcm(:) = zconemcm(:) * moyhri
               zconemom(:) = zconemom(:) * moyhri
               zconepblm(:) = zconepblm(:) * moyhri
               zconephym(:) = zconephym(:) * moyhri
               zconeradm(:) = zconeradm(:) * moyhri
               zconescm(:) = zconescm(:) * moyhri
               zconetotm(:) = zconetotm(:) * moyhri
               zconqdynm(:) = zconqdynm(:) * moyhri
               zconqcndm(:) = zconqcndm(:) * moyhri
               zconqdcm(:) = zconqdcm(:) * moyhri
               zconqmcm(:) = zconqmcm(:) * moyhri
               zconqmom(:) = zconqmom(:) * moyhri
               zconqpblm(:) = zconqpblm(:) * moyhri
               zconqphym(:) = zconqphym(:) * moyhri
               zconqradm(:) = zconqradm(:) * moyhri
               zconqscm(:) = zconqscm(:) * moyhri
               zconqtotm(:) = zconqtotm(:) * moyhri
            endif
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

               zo3tcm  (i) = zo3tcm  (i) + zo3tc  (i)
               zo3ctcm (i) = zo3ctcm (i) + zo3ctc (i)

               if (lavg) then
                  zo3tcm  (i) = zo3tcm  (i) * moyhri
                  zo3ctcm (i) = zo3ctcm (i) * moyhri
               end if

            end do

            do k = 1, nk-1
               do i = 1, ni

                  zo3avg  (i,k) = zo3avg  (i,k) + zo3lplus (i,k)

                  zo3ccolm (i,k) = zo3ccolm (i,k) + zo3ccol (i,k)
                  zo3colm  (i,k) = zo3colm  (i,k) + zo3col  (i,k)

                  zo1chmtdm  (i,k) =  zo1chmtdm (i,k) +  zo1chmtd (i,k)
                  zo4chmtdm  (i,k) =  zo4chmtdm (i,k) +  zo4chmtd (i,k)
                  zo6chmtdm  (i,k) =  zo6chmtdm (i,k) +  zo6chmtd (i,k)
                  zo7chmtdm  (i,k) =  zo7chmtdm (i,k) +  zo7chmtd (i,k)
                  zo3chmtdm  (i,k) =  zo3chmtdm (i,k) +  zo3chmtd (i,k)

                  if (lavg) then
                     zo3avg  (i,k) = zo3avg  (i,k) * moyhri

                     zo3ccolm (i,k) = zo3ccolm (i,k) * moyhri
                     zo3colm  (i,k) = zo3colm  (i,k) * moyhri

                     zo1chmtdm  (i,k) =  zo1chmtdm (i,k) * moyhri
                     zo4chmtdm  (i,k) =  zo4chmtdm (i,k) * moyhri
                     zo6chmtdm  (i,k) =  zo6chmtdm (i,k) * moyhri
                     zo7chmtdm  (i,k) =  zo7chmtdm (i,k) * moyhri
                     zo3chmtdm  (i,k) =  zo3chmtdm (i,k) * moyhri
                  endif

               end do
            end do
         end if IF_LINOZ

         IF_LINGH: if (llingh .and. out_linoz) then
            ! 2D
            do i = 1, ni

               zch4tcm (i) = zch4tcm (i) + zch4tc (i)
               zn2otcm (i) = zn2otcm (i) + zn2otc (i)
               zf11tcm (i) = zf11tcm (i) + zf11tc (i)
               zf12tcm (i) = zf12tcm (i) + zf12tc (i)

               if (lavg) then
                  zch4tcm (i) = zch4tcm (i) * moyhri
                  zn2otcm (i) = zn2otcm (i) * moyhri
                  zf11tcm (i) = zf11tcm (i) * moyhri
                  zf12tcm (i) = zf12tcm (i) * moyhri
               end if

            end do

            ! 3D
            do k = 1, nk-1
               do i = 1, ni
                  zch4avg (i,k) = zch4avg (i,k) + zch4lplus (i,k)
                  zn2oavg (i,k) = zn2oavg (i,k) + zn2olplus (i,k)
                  zf11avg (i,k) = zf11avg (i,k) + zf11lplus (i,k)
                  zf12avg (i,k) = zf12avg (i,k) + zf12lplus (i,k)

                  zch4colm (i,k) = zch4colm (i,k) + zch4col (i,k)
                  zn2ocolm (i,k) = zn2ocolm (i,k) + zn2ocol (i,k)
                  zf11colm (i,k) = zf11colm (i,k) + zf11col (i,k)
                  zf12colm (i,k) = zf12colm (i,k) + zf12col (i,k)

                  zch4chmtdm (i,k) = zch4chmtdm (i,k) + zch4chmtd (i,k)
                  zn2ochmtdm (i,k) = zn2ochmtdm (i,k) + zn2ochmtd (i,k)
                  zf11chmtdm (i,k) = zf11chmtdm (i,k) + zf11chmtd (i,k)
                  zf12chmtdm (i,k) = zf12chmtdm (i,k) + zf12chmtd (i,k)

                  if (lavg) then
                     zch4avg (i,k) = zch4avg (i,k) * moyhri
                     zn2oavg (i,k) = zn2oavg (i,k) * moyhri
                     zf11avg (i,k) = zf11avg (i,k) * moyhri
                     zf12avg (i,k) = zf12avg (i,k) * moyhri

                     zch4colm (i,k) = zch4colm (i,k) * moyhri
                     zn2ocolm (i,k) = zn2ocolm (i,k) * moyhri
                     zf11colm (i,k) = zf11colm (i,k) * moyhri
                     zf12colm (i,k) = zf12colm (i,k) * moyhri

                     zch4chmtdm (i,k) = zch4chmtdm (i,k) * moyhri
                     zn2ochmtdm (i,k) = zn2ochmtdm (i,k) * moyhri
                     zf11chmtdm (i,k) = zf11chmtdm (i,k) * moyhri
                     zf12chmtdm (i,k) = zf12chmtdm (i,k) * moyhri
                  endif
                  
               end do
            end do
         end if IF_LINGH

         do k = 1, nk-1
            do i = 1, ni

               zccnm   (i,k) = zccnm   (i,k) + zftot(i,k)
               ztim    (i,k) = ztim    (i,k) + zti  (i,k)
               zt2m    (i,k) = zt2m    (i,k) + zt2  (i,k)
               ztgwdm (i,k)  = ztgwdm (i,k)  + ztgwd(i,k)
               zutofdm (i,k) = zutofdm (i,k) + zutofd(i,k)
               zvtofdm (i,k) = zvtofdm (i,k) + zvtofd(i,k)
               zttofdm (i,k) = zttofdm (i,k) + zttofd(i,k)
               zugwdm  (i,k) = zugwdm  (i,k) + zugwd(i,k)
               zvgwdm  (i,k) = zvgwdm  (i,k) + zvgwd(i,k)
               zugnom  (i,k) = zugnom  (i,k) + zugno(i,k)
               zvgnom  (i,k) = zvgnom  (i,k) + zvgno(i,k)
               ztgnom  (i,k) = ztgnom  (i,k) + ztgno(i,k)
               zudifvm (i,k) = zudifvm (i,k) + zudifv(i,k)
               zvdifvm (i,k) = zvdifvm (i,k) + zvdifv(i,k)
               ztdifvm (i,k) = ztdifvm (i,k) + ztdifv(i,k)
               zqdifvm (i,k) = zqdifvm (i,k) + zqdifv(i,k)
               ztadvm  (i,k) = ztadvm  (i,k) + ztadv (i,k)
               zuadvm  (i,k) = zuadvm  (i,k) + zuadv (i,k)
               zvadvm  (i,k) = zvadvm  (i,k) + zvadv (i,k)
               zqadvm  (i,k) = zqadvm  (i,k) + zqadv (i,k)
               zqmetoxm(i,k) = zqmetoxm(i,k) + zqmetox(i,k)
               zhushalm(i,k) = zhushalm(i,k) + zhushal(i,k)
               ztshalm (i,k) = ztshalm (i,k) + ztshal(i,k)
               zlwcm   (i,k) = zlwcm   (i,k) + zlwc(i,k)
               ziwcm   (i,k) = ziwcm   (i,k) + ziwc(i,k)
               zlwcradm(i,k) = zlwcradm(i,k) + zlwcrad(i,k)
               ziwcradm(i,k) = ziwcradm(i,k) + ziwcrad(i,k)
               zcldradm(i,k) = zcldradm(i,k) + zcldrad(i,k)
               ztphytdm(i,k) = ztphytdm(i,k) + ztphytd(i,k)
               zhuphytdm(i,k)= zhuphytdm(i,k)+ zhuphytd(i,k)
               zuphytdm(i,k) = zuphytdm(i,k) + zuphytd(i,k)
               zvphytdm(i,k) = zvphytdm(i,k) + zvphytd(i,k)
               zzctem(i,k)   = zzctem(i,k)   + zcte(i,k)
               zzstem(i,k)   = zzstem(i,k)   + zste(i,k)
               zzcqem(i,k)   = zzcqem(i,k)   + zcqe(i,k)
               zzsqem(i,k)   = zzsqem(i,k)   + zsqe(i,k)
               zzcqcem(i,k)  = zzcqcem(i,k)  + zcqce(i,k)
               zzsqcem(i,k)  = zzsqcem(i,k)  + zsqce(i,k)

               if (associated(zmqem).and.associated(zmqe)) zmqem(i,k) = zmqem(i,k) + zmqe(i,k) 
               if (associated(zmtem).and.associated(zmte)) zmtem(i,k) = zmtem(i,k) + zmte(i,k) 
               if (associated(zumim).and.associated(zumid)) zumim(i,k) = zumim(i,k) + zumid(i,k) 
               if (associated(zvmim).and.associated(zvmid)) zvmim(i,k) = zvmim(i,k) + zvmid(i,k) 


               IF_AVG_1: if (lavg) then
                  zccnm   (i,k) = zccnm   (i,k) * moyhri
                  ztim    (i,k) = ztim    (i,k) * moyhri
                  zt2m    (i,k) = zt2m    (i,k) * moyhri
                  ztgwdm (i,k) = ztgwdm (i,k) * moyhri
                  zutofdm (i,k) = zutofdm (i,k) * moyhri
                  zvtofdm (i,k) = zvtofdm (i,k) * moyhri
                  zttofdm (i,k) = zttofdm (i,k) * moyhri
                  zugwdm  (i,k) = zugwdm  (i,k) * moyhri
                  zvgwdm  (i,k) = zvgwdm  (i,k) * moyhri
                  zugnom  (i,k) = zugnom  (i,k) * moyhri
                  zvgnom  (i,k) = zvgnom  (i,k) * moyhri
                  ztgnom  (i,k) = ztgnom  (i,k) * moyhri
                  zudifvm (i,k) = zudifvm (i,k) * moyhri
                  zvdifvm (i,k) = zvdifvm (i,k) * moyhri
                  ztdifvm (i,k) = ztdifvm (i,k) * moyhri
                  zqdifvm (i,k) = zqdifvm (i,k) * moyhri
                  ztadvm  (i,k) = ztadvm  (i,k) * moyhri
                  zuadvm  (i,k) = zuadvm  (i,k) * moyhri
                  zvadvm  (i,k) = zvadvm  (i,k) * moyhri
                  zqadvm  (i,k) = zqadvm  (i,k) * moyhri
                  zqmetoxm(i,k) = zqmetoxm(i,k) * moyhri
                  zzctem  (i,k) = zzctem  (i,k) * moyhri
                  zzcqem  (i,k) = zzcqem  (i,k) * moyhri
                  zzcqcem (i,k) = zzcqcem (i,k) * moyhri
                  zzstem  (i,k) = zzstem  (i,k) * moyhri
                  zzsqem  (i,k) = zzsqem  (i,k) * moyhri
                  zzsqcem (i,k) = zzsqcem (i,k) * moyhri
                  zhushalm(i,k) = zhushalm(i,k) * moyhri
                  ztshalm (i,k) = ztshalm (i,k) * moyhri
                  zlwcm   (i,k) = zlwcm   (i,k) * moyhri
                  ziwcm   (i,k) = ziwcm   (i,k) * moyhri
                  zlwcradm(i,k) = zlwcradm(i,k) * moyhri
                  ziwcradm(i,k) = ziwcradm(i,k) * moyhri
                  zcldradm(i,k) = zcldradm(i,k) * moyhri
                  ztphytdm(i,k) = ztphytdm(i,k) * moyhri
                  zhuphytdm(i,k)= zhuphytdm(i,k)* moyhri
                  zuphytdm(i,k) = zuphytdm(i,k) * moyhri
                  zvphytdm(i,k) = zvphytdm(i,k) * moyhri
                  if (associated(zmqem).and.associated(zmqe)) zmqem(i,k) = zmqem(i,k) * moyhri 
                  if (associated(zmtem).and.associated(zmte)) zmtem(i,k) = zmtem(i,k) * moyhri
                  if (associated(zumim).and.associated(zumid)) zumim(i,k) = zumim(i,k)* moyhri
                  if (associated(zvmim).and.associated(zvmid)) zvmim(i,k) = zvmim(i,k)* moyhri
               endif IF_AVG_1

            end do
         end do

         !# Note: all convec str in the if below must have same len
         if (any(convec == (/ &
              'KFC     ', &
              'KFC2    ', &
              'BECHTOLD'  &
              /))) then
            do k = 1, nk-1
               do i = 1, ni
                  ztfcpm  (i,k) = ztfcpm  (i,k) + ztfcp (i,k)
                  zhufcpm (i,k) = zhufcpm (i,k) + zhufcp(i,k)
                  zqckfcm (i,k) = zqckfcm (i,k) + zqckfc(i,k)
                  zumfkfcm(i,k) = zumfkfcm(i,k) + zumfkfc(i,k)
                  zdmfkfcm(i,k) = zdmfkfcm(i,k) + zdmfkfc(i,k)
               enddo
            enddo
            zudcm(:,1:nk-1)   = zudcm(:,1:nk-1)   + zufcp(:,1:nk-1)
            zvdcm(:,1:nk-1)   = zvdcm(:,1:nk-1)   + zvfcp(:,1:nk-1)
            if (cmt_comp_diag) then
               zufcp1m(:,1:nk-1) = zufcp1m(:,1:nk-1) + zufcp1(:,1:nk-1)
               zufcp2m(:,1:nk-1) = zufcp2m(:,1:nk-1) + zufcp2(:,1:nk-1)
               zufcp3m(:,1:nk-1) = zufcp3m(:,1:nk-1) + zufcp3(:,1:nk-1)
               zsufcpm(:,1:nk-1) = zsufcpm(:,1:nk-1) + zsufcp(:,1:nk-1)
               zvfcp1m(:,1:nk-1) = zvfcp1m(:,1:nk-1) + zvfcp1(:,1:nk-1)
               zvfcp2m(:,1:nk-1) = zvfcp2m(:,1:nk-1) + zvfcp2(:,1:nk-1)
               zvfcp3m(:,1:nk-1) = zvfcp3m(:,1:nk-1) + zvfcp3(:,1:nk-1)
               zsvfcpm(:,1:nk-1) = zsvfcpm(:,1:nk-1) + zsvfcp(:,1:nk-1)
            endif
            if (associated(ztusc) .and. associated(zuscm)) &
                 zuscm(:,1:nk-1) = zuscm(:,1:nk-1) + ztusc(:,1:nk-1)
            if (associated(ztvsc) .and. associated(zvscm)) &
                 zvscm(:,1:nk-1) = zvscm(:,1:nk-1) + ztvsc(:,1:nk-1)
            IF_AVG_2: if (lavg) then
               ztfcpm(:,1:nk-1)   = ztfcpm(:,1:nk-1)   * moyhri
               zhufcpm(:,1:nk-1)  = zhufcpm(:,1:nk-1)  * moyhri
               zqckfcm(:,1:nk-1)  = zqckfcm(:,1:nk-1)  * moyhri
               zumfkfcm(:,1:nk-1) = zumfkfcm(:,1:nk-1) * moyhri
               zdmfkfcm(:,1:nk-1) = zdmfkfcm(:,1:nk-1) * moyhri
               zudcm(:,1:nk-1)    = zudcm(:,1:nk-1) * moyhri
               zvdcm(:,1:nk-1)    = zvdcm(:,1:nk-1) * moyhri
               if (cmt_comp_diag) then
                  zufcp1m(:,1:nk-1) = zufcp1m(:,1:nk-1) * moyhri
                  zufcp2m(:,1:nk-1) = zufcp2m(:,1:nk-1) * moyhri
                  zufcp3m(:,1:nk-1) = zufcp3m(:,1:nk-1) * moyhri
                  zsufcpm(:,1:nk-1) = zsufcpm(:,1:nk-1) * moyhri
                  zvfcp1m(:,1:nk-1) = zvfcp1m(:,1:nk-1) * moyhri
                  zvfcp2m(:,1:nk-1) = zvfcp2m(:,1:nk-1) * moyhri
                  zvfcp3m(:,1:nk-1) = zvfcp3m(:,1:nk-1) * moyhri
                  zsvfcpm(:,1:nk-1) = zsvfcpm(:,1:nk-1) * moyhri
               endif
               if (associated(zuscm)) zuscm(:,1:nk-1) = zuscm(:,1:nk-1) * moyhri
               if (associated(zvscm)) zvscm(:,1:nk-1) = zvscm(:,1:nk-1) * moyhri
            endif IF_AVG_2

            do i=1, ni
               zabekfcm  (i) = zabekfcm(i) + zabekfc(i)
               zcapekfcm (i) = zcapekfcm(i) + zcapekfc(i)
               zcinkfcm (i) = zcinkfcm(i) + zcinkfc(i)
               zwumkfcm  (i) = zwumkfcm(i) + zwumaxkfc(i)
               zzbaskfcm (i) = zzbaskfcm(i) + zzbasekfc(i)
               zztopkfcm (i) = zztopkfcm(i) + zztopkfc(i)
               zkkfcm    (i) = zkkfcm(i) + zkkfc(i)
               if (associated(zkmidm)) zkmidm(i) = zkmidm(i) + zkmid(i)

               if (lavg) then
                  zabekfcm (i) = zabekfcm (i) * moyhri
                  zcapekfcm (i) = zcapekfcm (i) * moyhri
                  zcinkfcm (i) = zcinkfcm (i) * moyhri
                  zwumkfcm(i) = zwumkfcm(i) * moyhri
                  zzbaskfcm(i) = zzbaskfcm(i) * moyhri
                  zztopkfcm(i) = zztopkfcm(i) * moyhri
                  zkkfcm(i) = zkkfcm(i) * moyhri
                  if (associated(zkmidm)) zkmidm(i) = zkmidm(i) * moyhri
               endif

            end do
         endif

         if (fluvert /= 'NIL') then
            do i=1,ni
               zsumf(i) = zsumf(i) + zustress(i)
               zsvmf(i) = zsvmf(i) + zvstress(i)
               zfqm(i) = zfqm(i) + zfq(i)
            enddo
            if (lavg) then
               zsumf(:) = zsumf(:) * moyhri
               zsvmf(:) = zsvmf(:) * moyhri
               zfqm(:)  = zfqm(:)  * moyhri
            endif
         endif

      endif IF_ACCUM_1


      !****************************************************************
      !     ACCUMULATORS

      !Set accumulators to zero at the beginning and after every acchr hours,
      !and by default (acchr=0) as the model step goes through 0.
      IF_RESET_ACCUMULATORS: if (lkount0 .or. lacchr) then
         zrainaf(:) = 0.
         zsnowaf(:) = 0.
         if (radia /= 'NIL' .or. fluvert == 'SURFACE') then
            zeiaf(:) = 0.
            zevaf(:) = 0.
            zfiaf(:) = 0.
            zfsaf(:) = 0.
            zivaf(:) = 0.
            zntaf(:) = 0.
            zflusolaf(:) = 0.
         endif
         if (radia(1:8) == 'CCCMARAD') then
            zclbaf(:) = 0.
            zcltaf(:) = 0.
            zcstaf(:) = 0.
            zcsbaf(:) = 0.
            zfsdaf(:) = 0.
            zfsfaf(:) = 0.
            zfsiaf(:) = 0.
            zfsvaf(:) = 0.
            zparraf(:) = 0.
            zsiaf(:) = 0.
            zfnsaf(:) = 0.
         endif
         if (fluvert /= 'NIL') then
            zflaf(:) = 0.
            zfcaf(:) = 0.
            zfvaf(:) = 0.
         endif
         if (llight) then
            zafoudre(:) = 0.
         endif
      endif IF_RESET_ACCUMULATORS

      IF_KOUNT_NOT_0b: if (.not.lkount0) then
!VDIR NODEP
         DO_NI_ACC: do i = 1,ni

            !# Accumulation of precipitation (in m)

            zrainaf(i) = zrainaf(i) + zrainrate(i)*dt
            zsnowaf(i) = zsnowaf(i) + zsnowrate(i)*dt

            if (radia /= 'NIL' .or. fluvert == 'SURFACE') then
               zeiaf    (i) = zeiaf (i) + zei  (i) * dt
               zevaf    (i) = zevaf (i) + zev  (i) * dt
               zfiaf    (i) = zfiaf (i) + zfdsi(i) * dt
               zfsaf    (i) = zfsaf (i) + zfdss(i) * dt
               zivaf    (i) = zivaf (i) + ziv  (i) * dt
               zntaf    (i) = zntaf (i) + znt  (i) * dt
               zflusolaf(i) = zflusolaf(i) + zflusolis(i) * dt
            endif

            !# Accumulation of sfc and toa net clear sky fluxes, available with cccmarad
            if (radia(1:8) == 'CCCMARAD') then
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

            if (fluvert /= 'NIL') then
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

      !# For output purpose, diag level values needs to be copied into nk level of corresponding dynbus var
      zhuplus(:,nk) = zqdiag
      ztplus(:,nk)  = ztdiag
      zuplus(:,nk)  = zudiag
      zvplus(:,nk)  = zvdiag

      if (timings_L) call timing_stop_omp(470)
      call msg_toall(MSG_DEBUG, 'calcdiag [END]')
      !----------------------------------------------------------------
      return
   end subroutine calcdiag1

end module calcdiag
