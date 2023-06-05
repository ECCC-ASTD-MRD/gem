!-------------------------------------- LICENCE BEGIN --------------------------
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

module difver
   implicit none
   private
   public :: difver8

contains

   !/@*
   subroutine difver8(dbus, fbus, vbus, seloc, tau, &
        kount, trnch, ni, nk, nkm1)
      use debug_mod, only: init2nan
      use tdpack_const, only: CAPPA, CHLC, CHLF, CPD, DELTA, GRAV, RGASD
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use series_mod, only: series_xst
      use ens_perturb, only: ens_nc2d, ens_spp_get
      implicit none
!!!#include <arch_specific.hf>
      !@Object to perform the implicit vertical diffusion
      !@Arguments
      ! 002      A. Zadra (Oct 2015) - add land-water mask (MG) to input
      !                 list of baktotq4
      ! 003      A. Zadra / R. McT-C (Sep 2016) - added nonlocal scaling option
      !                 based on J. Mailhot/A. Lock (Aug 2012)

      !          - Input/Output -
      ! dbus     dynamic bus
      ! fbus     field for permanent physics variables
      ! vbus     volatile bus
      !          - output -
      ! tu       u  tendency
      ! tv       v  tendency
      ! tt       t  tendency
      ! tq       q  tendency
      ! tl       l tendency
      ! t        temperature
      ! uu       x component of the wind
      ! vv       y component of the wind
      ! q        specific humidity
      ! ql       liquid water
      ! ps       surface pressure
      ! tm       temperature at time t-dt
      !          - input -
      ! sg       sigma levels
      ! seloc    staggered sigma levels
      ! tau      timestep
      !          see common block "options"
      ! kount    timestep number
      ! trnch    row number
      ! ni       horizontal dimension
      ! nk       vertical dimension
      ! nkm1     vertical operator scope

      integer, intent(in) :: kount, trnch, ni, nk, nkm1
      real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)
      real, target :: seloc(ni,nkm1)
      real, intent(in) :: tau

      !@Author J. Cote (Oct 1984)
      !*@/
      real, parameter :: WEIGHT_BM=0.5

      integer j,k,typet,typem
      real rhortvsg,mrhocmu
      real tplusnk,qplusnk,uplusnk,vplusnk
      real rsg,du,dv,dsig, rgam
      logical :: compute_tend
#include "phymkptr.hf"
      include "surface.cdk"

      real, dimension(ni) :: bmsg, btsg, aq, bq, sfc_density, fsloflx, fm_mult, fh_mult
      real, dimension(nkm1) :: se, sig
      real, dimension(ni,nkm1) ::  kmsg, ktsg, wthl_ng, wqw_ng, uw_ng, vw_ng, &
           c, d, unused, zero, qclocal, ficelocal
      real, dimension(ni,nk) :: gam0
      real, dimension(ni,nkm1), target :: thl, qw, tthl, tqw, dkem, dket, tu2m, tu2t

      !     Pointeurs pour champs deja definis dans les bus
      real, pointer, dimension(:) :: zbt_ag, zbt_ag0, zfc_ag, zfv_ag, zqsurf_ag, zfvap_ag  !#TODO: should be contiguous
      real, pointer, dimension(:), contiguous   :: ps
      real, pointer, dimension(:), contiguous   :: zalfat, zalfat0, zalfaq, &
           zalfaq0, zbm, zbm0, zfdsi, zfdss, zfq, zmg, ztsrad, &
           zustress, zvstress, zue, zh
      real, pointer, dimension(:), contiguous :: zflw, zfsh
      real, pointer, dimension(:,:), contiguous :: tu, tv, tw, tt, tq, tl, uu, vv, w, &
           t, q, sg, zsigw, zsigt, zsigm, tm, &
           sigef, sigex, conserv_t,conserv_q,tconserv_t,tconserv_q, zgztherm, &
           zpblsigs,zpblq1, zqcplus, tqc
      real, pointer, dimension(:,:), contiguous :: zgq, zgql, zgte, zkm, zkt, zqtbl, ztve, &
           zwtng, zwqng, zuwng, zvwng, zfbl, zfblgauss, zfblnonloc, zfnn, &
           zzd, zzn, zfc, zfv, zc1pbl, zturbqf, zturbtf, zturbuf, zturbvf, &
           zturbuvf, zmrk2
      real, pointer, dimension(:,:,:), contiguous :: zvcoef

      !---------------------------------------------------------------------

      call init2nan(bmsg, btsg, aq, bq, sfc_density, fsloflx, se, sig)
      call init2nan(kmsg, ktsg, wthl_ng, wqw_ng, uw_ng, vw_ng, c, d, unused)
      call init2nan(zero, qclocal, ficelocal, gam0)
      call init2nan(thl, qw, tthl, tqw, dkem, dket)

      ! Pointer assignments
      MKPTR1D(ps, pmoins, fbus)

      if (sfcflx_filter_order <= 0) then
         MKPTR1D(zalfaq, alfaq, vbus)
         MKPTR1D(zalfat, alfat, vbus)
         MKPTR1DK(zbt_ag, bt, indx_agrege, vbus)
      else
         if (kount == 0) then
            MKPTR1D(zalfaq0, alfaq, vbus)
            MKPTR1D(zalfat0, alfat, vbus)
            MKPTR1DK(zbt_ag0, bt, indx_agrege, vbus)
         endif
         MKPTR1D(zalfaq, falfaq, fbus)
         MKPTR1D(zalfat, falfat, fbus)
         MKPTR1DK(zbt_ag, fbt, indx_agrege, fbus)
      endif
      
      MKPTR1D(zbm, bm, vbus)

      MKPTR1D(zbm0, bm0, fbus)
      MKPTR1D(zfdsi, fdsi, fbus)
      MKPTR1D(zfdss, fdss, fbus)
      MKPTR1D(zflw, flw, vbus)
      MKPTR1D(zfq, fq, fbus)
      MKPTR1D(zfsh, fsh, vbus)
      MKPTR1D(zh, h, fbus)
      MKPTR1D(zmg, mg, fbus)
      MKPTR1D(ztsrad, tsrad, fbus)
      MKPTR1D(zue, ue, vbus)
      MKPTR1D(zustress, ustress, vbus)
      MKPTR1D(zvstress, vstress, vbus)

      MKPTR1DK(zfc_ag, fc, indx_agrege, vbus)
      MKPTR1DK(zfv_ag, fv, indx_agrege, vbus)
      MKPTR1DK(zfvap_ag, fvap, indx_agrege, fbus)
      MKPTR1DK(zqsurf_ag, qsurf, indx_agrege, fbus)

      MKPTR2DN(zfc, fc, ni, nagrege, vbus)
      MKPTR2DN(zfv, fv, ni, nagrege, vbus)
      MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, fbus)

      MKPTR2D(zturbqf, turbqf, vbus)
      MKPTR2D(zturbtf, turbtf, vbus)
      MKPTR2D(zturbuf, turbuf, vbus)
      MKPTR2D(zturbuvf, turbuvf, vbus)
      MKPTR2D(zturbvf, turbvf, vbus)

      MKPTR2Dm1(q, huplus, dbus)
      MKPTR2Dm1(zqcplus, qcplus, dbus)
      MKPTR2Dm1(sg, sigm, dbus)
      MKPTR2Dm1(zsigm, sigm, dbus)
      MKPTR2Dm1(zsigt, sigt, dbus)
      MKPTR2Dm1(zsigw, sigw, dbus)
      MKPTR2Dm1(t, tplus, dbus)
      MKPTR2Dm1(tm, tmoins, dbus)
      MKPTR2Dm1(tu, udifv, vbus)
      MKPTR2Dm1(tv, vdifv, vbus)
      MKPTR2Dm1(tt, tdifv, vbus)
      MKPTR2Dm1(tq, qdifv, vbus)
      MKPTR2Dm1(tqc, qcdifv, vbus)
      MKPTR2Dm1(uu, uplus, dbus)
      MKPTR2Dm1(vv, vplus, dbus)
      MKPTR2Dm1(w, wplus, dbus)
      MKPTR2Dm1(zc1pbl, c1pbl, vbus)
      MKPTR2Dm1(zfbl, fbl, fbus)
      MKPTR2Dm1(zfblgauss, fblgauss, fbus)
      MKPTR2Dm1(zfblnonloc, fblnonloc, vbus)
      MKPTR2Dm1(zfnn, fnn, fbus)
      MKPTR2Dm1(zgq, gq, vbus)
      MKPTR2Dm1(zgql, gql, vbus)
      MKPTR2Dm1(zgte, gte, vbus)
      MKPTR2Dm1(zgztherm, gztherm, vbus)
      MKPTR2Dm1(zkm, km, vbus)
      MKPTR2Dm1(zkt, kt, vbus)
      MKPTR2Dm1(zqtbl, qtbl, fbus)
      MKPTR2Dm1(ztve, tve, vbus)
      MKPTR2Dm1(zuwng, uwng, vbus)
      MKPTR2Dm1(zvwng, vwng, vbus)
      MKPTR2Dm1(zwqng, wqng, vbus)
      MKPTR2Dm1(zwtng, wtng, vbus)
      MKPTR2Dm1(zzd, zd, fbus)
      MKPTR2Dm1(zzn, zn, fbus)

      MKPTR2Dm1(tw, wdifv, vbus)
      MKPTR2Dm1(tl, ldifv, vbus)
      MKPTR2Dm1(zpblsigs, pblsigs, vbus)
      MKPTR2Dm1(zpblq1, pblq1, vbus)

      MKPTR3D(zvcoef, vcoef, 2, vbus)

      ! Decide whether to compute PBL tendencies
      compute_tend = .true.
      if (fluvert == 'YSU') compute_tend = pbl_ysu_rpnsolve
      
      if (kount == 0 .and. sfcflx_filter_order > 0) then
         do j=1,ni
            zalfaq(j)=zalfaq0(j)
            zalfat(j)=zalfat0(j)
            zbt_ag(j)=zbt_ag0(j)
         enddo
      endif

      if(zsigt(1,1) > 0) then
         sigef(1:,1:) => zsigt(1:ni,1:nkm1)
         sigex(1:,1:) => zsigt(1:ni,1:nkm1)
      else
         sigef(1:,1:) => seloc(1:ni,1:nkm1)
         sigex(1:,1:) => zsigm(1:ni,1:nkm1)
      endif

      if (.not.diffuw) nullify(tw)

      ! Retrieve stochastic parameter information on request
      fm_mult(:) = ens_spp_get('fm_mult', zmrk2, default=1.)
      fh_mult(:) = ens_spp_get('fh_mult', zmrk2, default=1.)

      ! Choose diffusion operators based on vertical grid structure
      typem = 1
      if(zsigt(1,1) > 0) then
         if (tlift.eq.1) then
            typet = 6
         else
            typet = 5
         endif
      endif

      ! Conversion precalculations for sigma coordinates
      zero = 0.
      rsg = (GRAV/RGASD)
      do k=1,nkm1
!vdir nodep
         do j=1,ni
            gam0(j,k) = rsg*seloc(j,k)/ztve(j,k)
            kmsg(j,k) = zkm(j,k)*gam0(j,k)**2
            ktsg(j,k) = zkt(j,k)*gam0(j,k)**2
            NONLOCAL: if (pbl_nonloc == 'LOCK06') then
               ! Normalize the non-gradient flux terms (added to implicit solver)
               wthl_ng(j,k) = zwtng(j,k)*gam0(j,k)
               wqw_ng(j,k)  = zwqng(j,k)*gam0(j,k)
               uw_ng(j,k)   = zuwng(j,k)*gam0(j,k)
               vw_ng(j,k)   = zvwng(j,k)*gam0(j,k)
            else
               wthl_ng(j,k) = 0.0
               wqw_ng(j,k)  = 0.0
               uw_ng(j,k)   = 0.0
               vw_ng(j,k)   = 0.0
            endif NONLOCAL
            rgam = 1./gam0(j,k)
            zgte(j,k) = zgte(j,k) * rgam
            zgq(j,k)  = zgq(j,k)  * rgam
            zgql(j,k) = zgql(j,k) * rgam
         end do
      end do
      se = seloc(1,:)
      sig = sg(1,:)

      ! Start land surface fluxes slowly on user request over nsloflux time steps
      fsloflx = 1.
      if (nsloflux.gt.0) then
         do j=1,ni
            if (zmg(j)>0.5) then
               fsloflx(j) = (float(kount-1))/max(float(nsloflux),1.)
               if (kount.eq.0) fsloflx(j) = 0.0
            endif
         end do
      endif

      ! Compute turbulent flux surface boundary conditions
!vdir nodep
      do j=1,ni
         aq(j)=-rsg/(t(j,nkm1)*(1. + DELTA * zqsurf_ag(j)))
         bmsg(j)   = fm_mult(j) * zbm(j) * aq(j)
         btsg(j)   = zbt_ag(j) * aq(j) * fsloflx(j)
         zalfat(j) = zalfat(j)  * aq(j) * fsloflx(j)
         zalfaq(j) = zalfaq(j)  * aq(j) * fsloflx(j)
         sfc_density(j) = -aq(j) * ps(j)/GRAV
      end do
      if (pbl_cmu_timeavg .and. kount > 0) then
         bmsg(:) = fm_mult(:) * ((zbm(:)*WEIGHT_BM + zbm0(:)*(1.-WEIGHT_BM) )*aq(:))
         zbm0(:) = zbm(:)*WEIGHT_BM + zbm0(:)*(1.-WEIGHT_BM)
      endif
      gam0 = 0.

      ! ** U-Component Wind Diffusion **

      ! Compute tendency for u-component wind diffusion
      if (compute_tend) then
         call difuvdfj1(tu,uu,kmsg,zero,uw_ng,zero,zero,bmsg,sg,seloc,tau, &
              typem,1.,ni,ni,ni,nkm1)
         if (phy_error_L) return
      endif

      ! Diagnose atmospheric fluxes for u-component wind
      call atmflux3(zturbuf,uu,tu,kmsg,gam0,uw_ng,zero(1,nkm1), &
           bmsg,ps,t,q,tau,sg,seloc,zvcoef,c,d,0,ni,nkm1,trnch)

      ! ** V-Component Wind Diffusion **

      ! Compute tendency for v-component wind diffusion
      if (compute_tend) then
         call difuvdfj1(tv,vv,kmsg,zero,vw_ng,zero,zero,bmsg,sg,seloc,tau, &
              typem,1.,ni,ni,ni,nkm1)
         if (phy_error_L) return
      endif

      ! Diagnose atmospheric fluxes for v-component wind
      call atmflux3(zturbvf,vv,tv,kmsg,gam0,vw_ng,zero(1,nkm1), &
           bmsg,ps,t,q,tau,sg,seloc,zvcoef,c,d,1,ni,nkm1,trnch)

      ! Diagnose total turbulent momentum flux
      zturbuvf(:,:) = sqrt(zturbuf(:,:)**2+zturbvf(:,:)**2)

      ! ** Vertical Motion Diffusion **

      ! Compute tendency for vertical motion diffusion on user request
      if (diffuw .and. compute_tend) then
         call difuvdfj1(tw,w,kmsg,zero,zero,zero,zero,zero,sg,sigef,tau, &
              typet,1.,ni,ni,ni,nkm1)
         if (phy_error_L) return
      endif

      ! ** Moisture Diffusion **

      ! Set surface flux boundary terms
      aq(:) = fh_mult(:) * zalfaq(:)
      bq(:) = fh_mult(:) * btsg(:)

      ! Compute conserved moisture variable for vertical diffusion
      if (fluvert == 'MOISTKE') then
         ! Total water (qw)
         do k=1,nkm1
            do j=1,ni
               qclocal(j,k) = max( 0.0 , zqtbl(j,k))
               qw(j,k) = q(j,k) + max( 0.0 , qclocal(j,k) )
            enddo
         enddo
         conserv_q => qw
         tconserv_q => tqw
      else
         conserv_q => q
         tconserv_q => tq
      endif

      ! Compute tendency for moisture diffusion
      if (compute_tend) then
         call difuvdfj1(tconserv_q,conserv_q,ktsg,zgq,wqw_ng,zero,aq,bq,sg,sigef,tau, &
              typet,1.,ni,ni,ni,nkm1)
         if (phy_error_L) return
      endif

      ! Diagnose atmospheric fluxes for moisture
      call atmflux3(zturbqf,conserv_q,tconserv_q,ktsg,gam0,wqw_ng, &
           aq,bq,ps,t,q,tau,sg,sigef,zvcoef,c,d,2,ni,nkm1,trnch)

      ! ** Temperature Diffusion **

      ! Set surface flux boundary terms
      aq(:) = fh_mult(:) * zalfat(:)
      bq(:) = fh_mult(:) * btsg(:)

      ! Compute conserved temperature variable for vertical diffusion
      if (fluvert == 'MOISTKE') then
         ! Liquid water potential temperature (thl)
         call ficemxp2(ficelocal,unused,unused,tm,ni,nkm1)
         do k=1,nkm1
            do j=1,ni
               thl(j,k) = t(j,k) &
                    - ((CHLC+ficelocal(j,k)*CHLF)/CPD) &
                    *max( 0.0 , qclocal(j,k) )
               thl(j,k) = thl(j,k) * sigex(j,k)**(-CAPPA)
            end do
         end do
         conserv_t => thl
         tconserv_t => tthl
      else
         conserv_t => t
         tconserv_t => tt
      endif

      ! Compute tendency for temperature diffusion
      if (compute_tend) then
         call difuvdfj1(tconserv_t,conserv_t,ktsg,zgte,wthl_ng,zero,aq,bq,sg,sigef,tau, &
              typet,1.,ni,ni,ni,nkm1)
         if (phy_error_L) return
      endif
         
      ! Convert conserved variable tendencies to state (T,Q) variable tendencies
      if (fluvert == 'MOISTKE') then
         call baktotq8(tt,tq,tl,thl,qw,tthl,tqw,qclocal,sg,zsigw, &
              ps,zgztherm,tm,ficelocal,ztve,zh,zfc_ag,zqtbl,zfnn,zfbl,zfblgauss, &
              zfblnonloc,zc1pbl,zzn,zzd,zmg,zmrk2,zvcoef,zpblsigs,zpblq1,tau,ni,nkm1)
      endif

      ! Compute tendency for condensate diffusion
      if (pbl_diff_condens .and. compute_tend) then
         aq = 0.
         bq = 0.
         call difuvdfj1(tqc,zqcplus,ktsg,zgq,wqw_ng,zero,aq,bq,sg,sigef,tau, &
              typet,1.,ni,ni,ni,nkm1)
         if (phy_error_L) return
      else
         tqc = 0.
      endif

      ! Dissipative heating
      select case (pbl_dissheat)
      case ('LOCAL_TEND')
         dkem = (uu(:,1:nkm1)+0.5*tau*tu(:,1:nkm1)) * tu(:,1:nkm1) + &
              (vv(:,1:nkm1)+0.5*tau*tv(:,1:nkm1)) * tv(:,1:nkm1)
         call vint_mom2thermo(dket, dkem, zvcoef, ni, nkm1)
         do j=1,ni
            dsig = (zsigt(j,nkm1) - 1.)/(zsigt(j,nkm1-1) - 1.)
            dket(j,nkm1) = dsig * dket(j,nkm1-1)
         end do
      case ('LOCAL_K')
         tu2m = tu(:,1:nkm1)**2 + tv(:,1:nkm1)**2
         call vint_mom2thermo(tu2t, tu2m, zvcoef, ni, nkm1)
         do k=1,(nkm1-1)
            do j=1,ni
               du = uu(j,k) - uu(j,k+1) + tau*(tu(j,k) - tu(j,k+1))
               dv = vv(j,k) - vv(j,k+1) + tau*(tv(j,k) - tv(j,k+1))
               dsig = zsigm(j,k) - zsigm(j,k+1)
               dket(j,k) = - kmsg(j,k)*( (du/dsig)**2 + (dv/dsig)**2 ) &
                           - 0.5*tau*tu2t(j,k)
            end do
         end do
         do j=1,ni
            dsig = (zsigt(j,nkm1) - 1.)/(zsigt(j,nkm1-1) - 1.)
            dket(j,nkm1) = dsig * dket(j,nkm1-1)
         end do
      case DEFAULT
         dket = 0.
      end select
      tt(:,1:nkm1) = tt(:,1:nkm1) - (1./CPD) * dket(:,1:nkm1)

      ! Diagnose atmospheric fluxes for temperature
      gam0 = -GRAV/CPD !temperature is used in difuvdf instead of potential temperature
      call atmflux3(zturbtf,t,tt,ktsg,gam0,wthl_ng,aq,bq, &
           ps,t,q,tau,sg,sigef,zvcoef,c,d,3,ni,nkm1,trnch)

!vdir nodep
      do j = 1, ni
         rhortvsg = ps(j)/GRAV
         mrhocmu  = ps(j)/GRAV*bmsg(j)
         tplusnk  = t (j,nkm1)+tau*tt(j,nkm1)
         qplusnk  = q (j,nkm1)+tau*tq(j,nkm1)
         uplusnk  = uu(j,nkm1)+tau*tu(j,nkm1)
         vplusnk  = vv(j,nkm1)+tau*tv(j,nkm1)

         ! Recalculer les flux partout

         ! ustress et vstress sont calcules apres diffusion car
         ! on utilise toujours une formulation implicite pour
         ! la condition aux limites de surface pour les vents.
         ! Par contre, la formulation est explicite pour
         ! la temperature et l'humidite.
         ! A noter que, puisque la formulation est explicite,
         ! on agrege les flux fc et fv dans le sous-programme
         ! AGREGE; on pourrait aussi les calculer ici en tout
         ! temps, mais on ne le fait que pendant le "depart lent".
         ! Si on utilisait une formulation implicite pour
         ! fc et fv, il faudrait que ces derniers soient
         ! toujours calcules ici.

         if (nsloflux.gt.0.and.kount.le.nsloflux) then
            zfc_ag(j) = CPD * rhortvsg * (zalfat(j)+btsg(j)*tplusnk)
            zfv_ag(j) = CHLC * rhortvsg * (zalfaq(j)+btsg(j)*qplusnk)
         endif

         zustress(j) = -mrhocmu*uplusnk
         zvstress(j) = -mrhocmu*vplusnk
         zfq(j)      = -mrhocmu*sqrt(uplusnk**2 + vplusnk**2)
         zue(j)      = sqrt(zfq(j)/sfc_density(j))

         zfsh(j)    = zalfaq(j) + btsg(j)*qplusnk  ! units of s-1
         zflw(j)    = rhortvsg*zfsh(j)             ! water density flux - units of kg m-2 s-1
         zfsh(j)    = zflw(j)/sfc_density(j)       ! spec hum flux - units of ms-1
      enddo

      ! Time series diagnostics
      call series_xst(tu, 'tu', trnch)
      call series_xst(tv, 'tv', trnch)
      call series_xst(tt, 'tf', trnch)
      call series_xst(tq, 'qf', trnch)
      if (associated(tl)) &
           call series_xst(tl, 'lf', trnch)
      call series_xst(zfq, 'fq', trnch)
      call series_xst(zue, 'ue', trnch)
      !#TODO: check the following subarray calls
      call series_xst(zfc(:,indx_soil),    'f4', trnch)
      call series_xst(zfc(:,indx_glacier), 'f5', trnch)
      call series_xst(zfc(:,indx_water),   'f6', trnch)
      call series_xst(zfc(:,indx_ice),     'f7', trnch)
      call series_xst(zfc(:,indx_agrege),  'fc', trnch)
      call series_xst(zfv(:,indx_soil),    'h4', trnch)
      call series_xst(zfv(:,indx_glacier), 'h5', trnch)
      call series_xst(zfv(:,indx_water),   'h6', trnch)
      call series_xst(zfv(:,indx_ice),     'h7', trnch)
      call series_xst(zfv(:,indx_agrege),  'fv', trnch)
      call series_xst(zkm, 'km', trnch)
      call series_xst(zkt, 'kt', trnch)

      return
   end subroutine difver8

end module difver
