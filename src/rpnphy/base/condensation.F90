!#-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ----------------------------

module condensation
   implicit none
   private
   public :: condensation3

contains

   !/@*
   subroutine condensation3(d, dsiz, f, fsiz, v, vsiz, &
        tplus0, t0, huplus0, q0, qc0, ilab, dbdt, &
        dt, ni, nk, kount, trnch)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use tdpack_const, only: GRAV
      use energy_budget, only: eb_en,eb_pw,eb_residual_en,eb_residual_pw,eb_conserve_en,eb_conserve_pw,EB_OK
      use module_mp_p3,  only: mp_p3_wrapper_gem,n_qiType
      use mp_my2_mod,    only: mp_my2_main
      use my_dmom_mod,   only: mydmom_main
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use tendency, only: apply_tendencies
      use water_integrated, only: water_integrated1
      use ens_perturb, only: ens_nc2d
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Object Interface to convection/condensation
      !@Arguments
      !          - Input -
      ! dsiz     dimension of dbus
      ! fsiz     dimension of fbus
      ! vsiz     dimension of vbus
      ! dt       timestep (sec.)
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! kount    timestep number
      ! trnch    slice number
      ! tplus0   temperature at t+dT at the beginning of the physics
      ! huplus0  humidity at t+dT at the beginning of the physics
      !
      !          - Input/Output -
      ! d        dynamics input field
      ! f        historic variables for the physics
      ! v        physics tendencies and other output fields from the physics
      !
      !          - Output -
      ! t0       initial temperature at t+dT
      ! q0       initial humidity humidity  at t+dT
      ! qc0      initial total condensate mixing ratio at t+dT
      ! ilab     flag array: an indication of convective activity from Kuo schemes
      ! dbdt     estimated averaged cloud fraction growth rate for kuostd

      integer, intent(in) :: fsiz,vsiz,dsiz,ni,nk,kount,trnch
      real,    intent(in) :: dt
      integer,dimension(ni,nk-1), intent(inout) :: ilab
      real,   dimension(ni), intent(inout)      :: dbdt
      real,   dimension(ni,nk-1), intent(inout) :: tplus0,t0,huplus0,q0,qc0
      real,   target, intent(inout)             :: f(fsiz), v(vsiz), d(dsiz)

      !@Author L.Spacek, November 2011
      !@Revisions
      ! 001      new arguments in call to mydmom_main
      ! 002      PV-nov2014: fix communication between deep convection and MY_Dm
      !*@/
#include <msg.h>
      include "surface.cdk"

      logical, parameter :: NK_BOTTOM = .true.  !(.T. for nk at bottom)
      integer, parameter :: N_DIAG_2D = 20      !number of diagnostic 2D fields
      integer, parameter :: N_DIAG_3D = 20      !number of diagnostic 3D fields

      integer :: nkm1, istat1, istat2, i
      real :: idt
      logical :: p3_comptend

      real, dimension(ni) :: tlcr,tscr,lv,ls

      real, dimension(ni,nk-1) :: zfm,zfm1,zcqer,zcqcer,zcter,iwc_total,lqip,lqrp,lqgp,lqnp,geop, lttp,lhup,ccf

      real, dimension(ni,N_DIAG_2D)      :: diag_2d    !diagnostic 2D fields
      real, dimension(ni,nk-1,N_DIAG_3D) :: diag_3d    !diagnostic 3D fields
      real, dimension(ni,nk-1,n_qiType)  :: qi_type    !diagnostic ice particle type  (mp_p3)

      real(REAL64), dimension(ni) :: l_en0,l_en,l_pw0,l_pw,l_enr,l_pwr

      real, dimension(:,:), pointer :: zttm,zhum

#define PHYPTRDCL
#include "condensation_ptr.hf"

      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'condensation [BEGIN]')

      nkm1 = nk -1
      idt  = 1./dt

#undef PHYPTRDCL
#include "condensation_ptr.hf"

      call init2nan(tlcr, tscr, lv, ls)
      call init2nan(zfm, zfm1, zcqer, zcqcer, zcter, iwc_total, lqip, lqrp, lqgp)
      call init2nan(lqnp, geop, lttp, lhup)
      call init2nan(diag_2d)
      call init2nan(diag_3d)
      call init2nan(l_en0, l_en, l_pw0, l_pw, l_enr, l_pwr)

      ! Startup operations
      if (kount == 0) then
         zhupostcnd = qqm
         zqcpostcnd = qcm
         ztpostcnd = ttm
         if (associated(zhums) .and. associated(zttms)) then
            zhums = zhupostcnd
            zttms = ztpostcnd
         endif
      endif

      ! Local initializations
      zste = 0.
      zsqe = 0.
      zsqce = 0.
      zsqre = 0.
      ccf = zfdc + zfmc  !sum convective cloud fractions (reasonable b/c of vertical correlation)

      ! Compute tendencies from P3 microphysics only on demand
      p3_comptend = .false.
      i = 1
      do while (.not.p3_comptend .and. i <= nphyoutlist)
         if (any(phyoutlist_S(i) == (/ &
              'ste ','sqe ','sqce','sqre','mtqi', &
              'w7  ','w9  ','w5  ','ta  ' &
              /))) p3_comptend = .true.
         i = i+1
      enddo

      ! Run selected gridscale condensation scheme
      GRIDSCALE_SCHEME: select case(stcond)

      case('NEWSUND')

         !# sundqvist (deuxieme version) :
         call skocon5(zste, zsqe, zsqce, a_tls, &
              a_tss, a_fxp, ttp, ttm, qqp, &
              qqm, qcp, qcm, psp, &
              psm, ilab, sigma, ni, nkm1, &
              dt, satuco, convec, zrnflx, zsnoflx)

      case('CONSUN')

         tlcr = 0.
         tscr = 0.
         zcter = 0.
         zcqer = 0.
         zcqcer = 0.
         if (convec=='KUOSTD') then
            !# transvider les tendances convectives pour kuostd
            !# par contre, on ne veut pas d'interaction
            !# entre les schemas kfc et consun.
            zcter = zcte
            zcqer = zcqe
         endif

         ! Apply physics tendencies to smoothed thermodynamic fields if available
         if (associated(zttps) .and. associated(zhups) .and. associated(zttms) .and. associated(zhums)) then
            zttm => zttms
            zhum => zhums
            lttp = zttps + (ttp-tplus0)
            lhup = zhups + (qqp-huplus0)
         else
            zttm => ztpostcnd
            zhum => zhupostcnd
            lttp = ttp
            lhup = qqp
         endif
         zfm1 = zqcpostcnd
         zfm  = qcp
         call consun3(zste , zsqe , zsqce , a_tls, a_tss, a_fxp, &
              zcter, zcqer, zcqcer, tlcr  , tscr  , ccf, &
              lttp    , zttm   , lhup     , zhum    , zfm   , zfm1  , &
              psp , psm  , ilab  , dbdt  , sigma, dt  , &
              zrnflx, zsnoflx, zf12 , zfevp  , &
              zfice, zmrk2, ni , nkm1)

         ! Adjust tendencies to impose conservation of total water and liquid
         ! water static energy on request.
         ! Apply humidity tendency correction for total water conservation
         ! and of liquid water static energy
         CONSUN_CONSERVATION: if (cond_conserve == 'TEND') then
            istat1 = eb_conserve_pw(zsqe,zsqe,ttp,qqp,sigma,psp,nkm1,F_dqc=zsqce,F_rain=a_tls,F_snow=a_tss)
            istat2 = eb_conserve_en(zste,zste,zsqe,ttp,qqp,qcp,sigma,psp,nkm1,F_dqc=zsqce,F_rain=a_tls,F_snow=a_tss)
            if (istat1 /= EB_OK .or. istat2 /= EB_OK) then
               call physeterror('condensation', 'Problem correcting for liquid water static energy conservation in condensation')
               return
            endif
         endif CONSUN_CONSERVATION

         !# transvider les tendances convectives et les taux
         !# des precipitations pour kuostd
         !# transvider les tendances de t et hu ainsi que la fraction nuageuse.
         !# amalgamer les champs de sortie de kuosym et de fcp.
         if (convec == 'KUOSTD') then
            ztsc(:) =  tscr(:)
            ztlc(:) =  tlcr(:)
            zcqce = zcqcer
            zcte  = zcter
            zcqe  = zcqer
         endif

      case('MP_MY2_OLD')

         !# "old" Milbrandt-Yau 2-moment microphysics (as in HRDPS v4.0.0 - v4.2.0)
         geop = zgztherm*GRAV

         call mydmom_main(ww,&
              ttp,qqp,qcp,qrp,qip,qnp,qgp,qhp,ncp,nrp,nip,nnp,ngp,nhp,psp, &
              ttm,qqm,qcm,qrm,qim,qnm,qgm,qhm,ncm,nrm,nim,nnm,ngm,nhm,psm,sigma,&
              a_tls_rn1, a_tls_rn2, a_tls_fr1, a_tls_fr2,&
              a_tss_sn1, a_tss_sn2, a_tss_sn3,&
              a_tss_pe1, a_tss_pe2, a_tss_pe2l, a_tss_snd,geop,&
              zste, zsqe, zsqce,zsqre,&
              qitend, qntend, qgtend, qhtend, nctend,&
              nrtend, nitend, nntend, ngtend, nhtend,&
              dt, ni, ni, nkm1, trnch, kount, my_ccntype, &
              my_diagon, my_sedion, my_warmon, &
              my_rainon, my_iceon,  my_snowon, my_initn,&
              my_dblmom_c, my_dblmom_r, my_dblmom_i, &
              my_dblmom_s, my_dblmom_g, my_dblmom_h,&
              a_dm_c, a_dm_r, a_dm_i, a_dm_s, a_dm_g, a_dm_h,&
              a_zet,  a_zec,  a_slw,  a_vis,  a_vis1, a_vis2, a_vis3,&
              a_h_cb, a_h_ml, a_h_m2, a_h_sn, &
              a_ss01, a_ss02, a_ss03, a_ss04, a_ss05, a_ss06, a_ss07,&
              a_ss08, a_ss09, a_ss10, a_ss11, a_ss12, a_ss13, a_ss14,&
              a_ss15, a_ss16, a_ss17, a_ss18, a_ss19, a_ss20, my_tc3comp,  &
              zrnflx,zsnoflx,zf12,zfevp,a_fxp,a_effradc,&
              a_effradi1, a_effradi2, a_effradi3, a_effradi4)

      case('MP_MY2')

         !#  modified (from HRDPS_4.2.0) Milbrandt-Yau 2-moment microphysics (MY2, v2.25.2)
         call mp_my2_main(ww,ttp,qqp,qcp,qrp,qip,qnp,qgp,qhp,ncp,nrp,nip,nnp,ngp,nhp,a_nwfa,     &
              psp, sigma, a_tls_rn1, a_tls_rn2, a_tls_fr1, a_tls_fr2, a_tss_sn1, a_tss_sn2,      &
              a_tss_sn3, a_tss_pe1, a_tss_pe2, a_tss_pe2l, a_tss_snd,dt, ni, nkm1, 1, kount,     &
              mp_aeroact, my_ccntype, my_diagon, my_sedion, my_warmon, my_rainon, my_iceon,      &
              my_snowon, a_dm_c,a_dm_r,a_dm_i,a_dm_s,a_dm_g,a_dm_h, a_zet, a_zec,diag_3d,        &
              a_effradc, a_effradi1, a_effradi2, a_effradi3, a_effradi4, a_fxp, NK_BOTTOM)
         if (phy_error_L) return

         !-- temporary:  (until RN/FR separation gets removed in MY2)
         a_tls_rn1 = a_tls_rn1 + a_tls_fr1
         a_tls_rn2 = a_tls_rn2 + a_tls_fr2
         a_tls_fr1 = 0.
         a_tls_fr2 = 0.

         a_ss01 = diag_3d(:,:,1);   a_ss08 = diag_3d(:,:, 8);   a_ss15 = diag_3d(:,:,15)
         a_ss02 = diag_3d(:,:,2);   a_ss09 = diag_3d(:,:, 9);   a_ss16 = diag_3d(:,:,16)
         a_ss03 = diag_3d(:,:,3);   a_ss10 = diag_3d(:,:,10);   a_ss17 = diag_3d(:,:,17)
         a_ss04 = diag_3d(:,:,4);   a_ss11 = diag_3d(:,:,11);   a_ss18 = diag_3d(:,:,18)
         a_ss05 = diag_3d(:,:,5);   a_ss12 = diag_3d(:,:,12);   a_ss19 = diag_3d(:,:,19)
         a_ss06 = diag_3d(:,:,6);   a_ss13 = diag_3d(:,:,13);   a_ss20 = diag_3d(:,:,20)
         a_ss07 = diag_3d(:,:,7);   a_ss14 = diag_3d(:,:,14)

      case('MP_P3')

        !save values before call, for computation of tendencies immediately after
        ! note: qitend, this is done below since it is p3_ncat-dependent
         if (p3_comptend) then
            zste(:,:)   = ttp(:,:)
            zsqe(:,:)   = qqp(:,:)
            zsqce(:,:)  = qcp(:,:)
            zsqre(:,:)  = qrp(:,:)
         endif

         !#  Predicted Particle Properties (P3) microphysics (v3.1.4)
         if (p3_ncat == 1) then

            if (p3_comptend) qitend(:,:) = qti1p(:,:)

            istat1 = mp_p3_wrapper_gem(qqm,qqp,ttm,ttp,dt,p3_dtmax,ww,psp,zgztherm,sigma,   &
                    kount,trnch,ni,nkm1,a_tls,a_tss,a_tls_rn1,a_tls_rn2,a_tss_sn1,          &
                    a_tss_sn2,a_tss_sn3,a_tss_pe1,a_tss_pe2,a_tss_snd,a_zet,a_zec,          &
                    a_effradc,qcp,ncp,qrp,nrp,p3_ncat,N_DIAG_2D,diag_2d,N_DIAG_3D,diag_3d,  &
                    qi_type,p3_depfact,p3_subfact,p3_debug,a_h_cb,a_h_sn,a_vis,a_vis1,      &
                    a_vis2,a_vis3,a_slw,p3_scpf_on,p3_pfrac,p3_resfact,a_fxp,p3_clip_qv,    &
                    qti1p,qmi1p,nti1p,bmi1p,a_effradi1)
            if (istat1 >= 0) then
               iwc_total = qti1p
               where (qti1p(:,1:nkm1)<1.e-14) a_effradi1 = 0.
            endif

         elseif (p3_ncat == 2) then

            if (p3_comptend) qitend(:,:) = qti1p(:,:) + qti2p(:,:)

            istat1 = mp_p3_wrapper_gem(qqm,qqp,ttm,ttp,dt,p3_dtmax,ww,psp,zgztherm,sigma,   &
                    kount,trnch,ni,nkm1,a_tls,a_tss,a_tls_rn1,a_tls_rn2,a_tss_sn1,          &
                    a_tss_sn2,a_tss_sn3,a_tss_pe1,a_tss_pe2,a_tss_snd,a_zet,a_zec,          &
                    a_effradc,qcp,ncp,qrp,nrp,p3_ncat,N_DIAG_2D,diag_2d,N_DIAG_3D,diag_3d,  &
                    qi_type,p3_depfact,p3_subfact,p3_debug,a_h_cb,a_h_sn,a_vis,a_vis1,      &
                    a_vis2,a_vis3,a_slw,p3_scpf_on,p3_pfrac,p3_resfact,a_fxp,p3_clip_qv,    &
                    qti1p,qmi1p,nti1p,bmi1p,a_effradi1,                                     &
                    qti2p,qmi2p,nti2p,bmi2p,a_effradi2)
            if (istat1 >= 0) then
               iwc_total = qti1p + qti2p
               where (qti1p(:,1:nkm1)<1.e-14) a_effradi1 = 0.
               where (qti2p(:,1:nkm1)<1.e-14) a_effradi2 = 0.
            endif

         elseif (p3_ncat == 3) then

            if (p3_comptend) qitend(:,:) = qti1p(:,:) + qti2p(:,:) + qti3p(:,:)

            istat1 = mp_p3_wrapper_gem(qqm,qqp,ttm,ttp,dt,p3_dtmax,ww,psp,zgztherm,sigma,   &
                    kount,trnch,ni,nkm1,a_tls,a_tss,a_tls_rn1,a_tls_rn2,a_tss_sn1,          &
                    a_tss_sn2,a_tss_sn3,a_tss_pe1,a_tss_pe2,a_tss_snd,a_zet,a_zec,          &
                    a_effradc,qcp,ncp,qrp,nrp,p3_ncat,N_DIAG_2D,diag_2d,N_DIAG_3D,diag_3d,  &
                    qi_type,p3_depfact,p3_subfact,p3_debug,a_h_cb,a_h_sn,a_vis,a_vis1,      &
                    a_vis2,a_vis3,a_slw,p3_scpf_on,p3_pfrac,p3_resfact,a_fxp,p3_clip_qv,    &
                    qti1p,qmi1p,nti1p,bmi1p,a_effradi1,                                     &
                    qti2p,qmi2p,nti2p,bmi2p,a_effradi2,                                     &
                    qti3p,qmi3p,nti3p,bmi3p,a_effradi3)
            if (istat1 >= 0) then
               iwc_total = qti1p + qti2p + qti3p
               where (qti1p(:,1:nkm1)<1.e-14) a_effradi1 = 0.
               where (qti2p(:,1:nkm1)<1.e-14) a_effradi2 = 0.
               where (qti3p(:,1:nkm1)<1.e-14) a_effradi3 = 0.
            endif

         elseif (p3_ncat == 4) then

            if (p3_comptend) qitend(:,:) = qti1p(:,:) + qti2p(:,:) + qti3p(:,:) + qti4p(:,:)

            istat1 = mp_p3_wrapper_gem(qqm,qqp,ttm,ttp,dt,p3_dtmax,ww,psp,zgztherm,sigma,   &
                    kount,trnch,ni,nkm1,a_tls,a_tss,a_tls_rn1,a_tls_rn2,a_tss_sn1,          &
                    a_tss_sn2,a_tss_sn3,a_tss_pe1,a_tss_pe2,a_tss_snd,a_zet,a_zec,          &
                    a_effradc,qcp,ncp,qrp,nrp,p3_ncat,N_DIAG_2D,diag_2d,N_DIAG_3D,diag_3d,  &
                    qi_type,p3_depfact,p3_subfact,p3_debug,a_h_cb,a_h_sn,a_vis,a_vis1,      &
                    a_vis2,a_vis3,a_slw,p3_scpf_on,p3_pfrac,p3_resfact,a_fxp,p3_clip_qv,    &
                    qti1p,qmi1p,nti1p,bmi1p,a_effradi1,                                     &
                    qti2p,qmi2p,nti2p,bmi2p,a_effradi2,                                     &
                    qti3p,qmi3p,nti3p,bmi3p,a_effradi3,                                     &
                    qti4p,qmi4p,nti4p,bmi4p,a_effradi4)
            if (istat1 >= 0) then
               iwc_total = qti1p + qti2p + qti3p + qti4p
               where (qti1p(:,1:nkm1)<1.e-14) a_effradi1 = 0.
               where (qti2p(:,1:nkm1)<1.e-14) a_effradi2 = 0.
               where (qti3p(:,1:nkm1)<1.e-14) a_effradi3 = 0.
               where (qti4p(:,1:nkm1)<1.e-14) a_effradi4 = 0.
            endif

         endif  !p3_ncat

         if (istat1 < 0) then
            call physeterror('condensation', 'Problem in p3')
            return
         endif

        !compute tendencies from microphysics:
         if (p3_comptend) then
            zste(:,:)  = (ttp(:,:)-zste(:,:) )*idt
            zsqe(:,:)  = (qqp(:,:)-zsqe(:,:) )*idt
            zsqce(:,:) = (qcp(:,:)-zsqce(:,:))*idt
            zsqre(:,:) = (qrp(:,:)-zsqre(:,:))*idt
            if (p3_ncat==1) then
               qitend(:,:) = (qti1p(:,:)-qitend(:,:))*idt
            elseif (p3_ncat==2) then
               qitend(:,:) = ((qti1p(:,:)+qti2p(:,:))-qitend(:,:))*idt
            elseif (p3_ncat==3) then
               qitend(:,:) = ((qti1p(:,:)+qti2p(:,:)+qti3p(:,:))-qitend(:,:))*idt
            elseif (p3_ncat==4) then
               qitend(:,:) = ((qti1p(:,:)+qti2p(:,:)+qti3p(:,:)+qti4p(:,:))-qitend(:,:))*idt
            endif
         endif

        !temporary; rn/fr (rate) partition should be done in s/r 'calcdiag'
        !(but currently it is still done inside microphyics scheme for MY2)
         a_tls_fr1 = 0.
         a_tls_fr2 = 0.

        !diagnostic ice particle types:
         a_qi_1 = qi_type(:,:,1)  !small ice crystals
         a_qi_2 = qi_type(:,:,2)  !unrimed snow crystals
         a_qi_3 = qi_type(:,:,3)  !lightly rimed snow
         a_qi_4 = qi_type(:,:,4)  !graupel
         a_qi_5 = qi_type(:,:,5)  !hail
         a_qi_6 = qi_type(:,:,6)  !ice pellets

         if (.true.) then  ! namelist switch to be added
            a_d2d01 = diag_2d(:,1);    a_d2d08 = diag_2d(:, 8);    a_d2d15 = diag_2d(:,15)
            a_d2d02 = diag_2d(:,2);    a_d2d09 = diag_2d(:, 9);    a_d2d16 = diag_2d(:,16)
            a_d2d03 = diag_2d(:,3);    a_d2d10 = diag_2d(:,10);    a_d2d17 = diag_2d(:,17)
            a_d2d04 = diag_2d(:,4);    a_d2d11 = diag_2d(:,11);    a_d2d18 = diag_2d(:,18)
            a_d2d05 = diag_2d(:,5);    a_d2d12 = diag_2d(:,12);    a_d2d19 = diag_2d(:,19)
            a_d2d06 = diag_2d(:,6);    a_d2d13 = diag_2d(:,13);    a_d2d20 = diag_2d(:,20)
            a_d2d07 = diag_2d(:,7);    a_d2d14 = diag_2d(:,14)

            a_ss01 = diag_3d(:,:,1);   a_ss08 = diag_3d(:,:, 8);   a_ss15 = diag_3d(:,:,15)
            a_ss02 = diag_3d(:,:,2);   a_ss09 = diag_3d(:,:, 9);   a_ss16 = diag_3d(:,:,16)
            a_ss03 = diag_3d(:,:,3);   a_ss10 = diag_3d(:,:,10);   a_ss17 = diag_3d(:,:,17)
            a_ss04 = diag_3d(:,:,4);   a_ss11 = diag_3d(:,:,11);   a_ss18 = diag_3d(:,:,18)
            a_ss05 = diag_3d(:,:,5);   a_ss12 = diag_3d(:,:,12);   a_ss19 = diag_3d(:,:,19)
            a_ss06 = diag_3d(:,:,6);   a_ss13 = diag_3d(:,:,13);   a_ss20 = diag_3d(:,:,20)
            a_ss07 = diag_3d(:,:,7);   a_ss14 = diag_3d(:,:,14)
         endif

      end select GRIDSCALE_SCHEME

      !# application des tendances convectives de qc (pour consun)
      if (any(stcond == (/'CONSUN','KUOSTD'/))) then
         call apply_tendencies(qcp, zcqce, ztdmask, ni, nk, nkm1)
      endif

      ! Pre-scheme state for energy budget
      if (associated(zconecnd)) then
         istat1 = eb_en(l_en0,ttp,qqp,qcp,sigma,psp,nkm1)
         istat2 = eb_pw(l_pw0,qqp,qcp,sigma,psp,nkm1)
         if (istat1 /= EB_OK .or. istat2 /= EB_OK) then
            call physeterror('condensation', 'Problem computing preliminary energy budget for '//trim(stcond))
            return
         endif
      endif

      !Application of general tendencies provided by most explicit schemes
      if (.not.any(stcond == (/'MP_MY2','MP_P3 '/))) then
         call apply_tendencies(ttp, zste,  ztdmask, ni, nk, nkm1)
         call apply_tendencies(qqp, zsqe,  ztdmask, ni, nk, nkm1)
         call apply_tendencies(qcp, zsqce, ztdmask, ni, nk, nkm1)
      endif

      ! Apply tendencies from microphysics scheme (my2_old only)
      if (stcond == 'MP_MY2_OLD') then
         call apply_tendencies(qrp, zsqre,  ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(qip, qitend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(qgp, qgtend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(qnp, qntend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(qhp, qhtend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(ncp, nctend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(nrp, nrtend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(nip, nitend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(ngp, ngtend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(nnp, ngtend, ztdmask, ni, nk, nkm1, F_minval=0.)
         call apply_tendencies(nhp, nhtend, ztdmask, ni, nk, nkm1, F_minval=0.)
      endif

      if (stcond(1:6) == 'MP_MY2') then
         lqrp = qrp
         lqip = qip
         lqgp = qgp
         lqnp = qnp
      elseif (stcond == 'MP_P3') then
         lqrp = qrp
         lqip = iwc_total
         lqgp = 0.
         lqnp = 0.
      else
         lqrp = 0.
         lqip = 0.
         lqgp = 0.
         lqnp = 0.
      endif

      ! Post-scheme energy budget analysis: post-scheme state and residuals
      if (associated(zconecnd)) then
         istat1 = eb_en(l_en,ttp,qqp,qcp,sigma,psp,nkm1)
         istat2 = eb_pw(l_pw,qqp,qcp,sigma,psp,nkm1)
         if (istat1 == EB_OK .and. istat2 == EB_OK) then
            istat1 = eb_residual_en(l_enr,l_en0,l_en,ttp,qqp,qcp,delt,nkm1,F_rain=a_tls,F_snow=a_tss)
            istat2 = eb_residual_pw(l_pwr,l_pw0,l_pw,ttp,delt,nkm1,F_rain=a_tls,F_snow=a_tss)
         endif
         if (istat1 /= EB_OK .or. istat2 /= EB_OK) then
            call physeterror('condensation', 'Problem computing final energy budget for '//trim(stcond))
            return
         endif
         zconecnd(:) = real(l_enr)
         zconqcnd(:) = real(l_pwr)
      endif

      ! Local copies of the MY2 masses are used to avoid passing a nullified
      ! pointer through the interface.
      call water_integrated1(ttp,qqp,qcp,lqip,lqnp,sigma,psp, &
           zicw,ziwv,ziwv700,ziwp,zlwp2,zslwp,zslwp2,zslwp3,zslwp4,ni,nkm1)

      ! Store post-scheme state information for "accession" calculations on next step
      zhupostcnd = qqp
      zqcpostcnd = qcp
      ztpostcnd = ttp

      ! Reset state variables for tendency application unless the microphysics
      ! scheme does not update state variables
      if (.not.any(stcond == (/'MP_MY2','MP_P3 '/))) then
         ttp(:,1:nkm1) = t0(:,1:nkm1)
         qqp(:,1:nkm1) = q0(:,1:nkm1)
         qcp(:,1:nkm1) = qc0(:,1:nkm1)
      endif

      ! Clip negative cloud water
      if (stcond == 'MP_MY2_OLD') qcp = max(qcp, 0.)

      call msg_toall(MSG_DEBUG, 'condensation [END]')
      !----------------------------------------------------------------
      return
   end subroutine condensation3

end module condensation
