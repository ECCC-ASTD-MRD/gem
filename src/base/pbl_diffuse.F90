module pbl_diffuse
  implicit none
  private

  public :: difver             !Manage diffusion of all variables
  public :: difuvdfj           !Diffusion operator

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine difver(pvars, tau, kount, ni, nk, nkm1, trnch)
    use debug_mod, only: init2nan
    use tdpack_const, only: CAPPA, CHLC, CHLF, CPD, DELTA, GRAV, RGASD
    use phy_options, except=>pbl_dissheat
    use phy_status, only: phy_error_L
    use phybusidx
    use phymem, only: phyvar
    use series_mod, only: series_xst
    use ens_perturb, only: ens_nc2d, ens_spp_get
    use atmflux, only: atmflux4
    use pbl_utils, only: ficemxp, dissheat, sfcflux
    use pbl_mtke_utils, only: baktotq
    use integrals, only: INT_TYPE_LINEAR
    implicit none
!!!#include <arch_specific.hf>
    !@Object to perform the implicit vertical diffusion
    !@Arguments
    ! 002      A. Zadra (Oct 2015) - add land-water mask (MG) to input
    !                 list of baktotq4
    ! 003      A. Zadra / R. McT-C (Sep 2016) - added nonlocal scaling option
    !                 based on J. Mailhot/A. Lock (Aug 2012)

    !          - Input/Output -
    ! pvars    list of all phy vars (meta + slab data)
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
    ! tau      timestep
    !          see common block "options"
    ! kount    timestep number
    ! trnch    row number
    ! ni       horizontal dimension
    ! nk       vertical dimension
    ! nkm1     vertical operator scope

    type(phyvar), pointer, contiguous :: pvars(:)
    integer, intent(in) :: kount, trnch, ni, nk, nkm1
    real, intent(in) :: tau

    !@Author J. Cote (Oct 1984)
    !*@/
    real, parameter :: WEIGHT_BM=0.5
    integer, parameter :: FLUX_INTTYPE=INT_TYPE_LINEAR

    integer i, j, k, typet, typem
    real rsg, rgam
    logical :: compute_tend
#include "phymkptr.hf"
    include "surface.cdk"

    real, dimension(ni) :: bmsg, btsg, aq, bq, sfc_density, fsloflx, fm_mult, fh_mult
    real, dimension(nkm1) :: se, sig
    real, dimension(ni,nkm1) ::  kmsg, ktsg, wthl_ng, wqw_ng, uw_ng, vw_ng, &
         c, d, unused, zero, qclocal, ficelocal
    real, dimension(ni,nk) :: gam0
    real, dimension(ni,nkm1), target :: thl, qw, tthl, tqw, dkem, dket

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
         zsigmas,zpblq1, zqcplus, tqc
    real, pointer, dimension(:,:), contiguous :: zgq, zgql, zgte, zkm, zkt, zqtbl, ztve, &
         zwtng, zwqng, zuwng, zvwng, zfbl, zfblgauss, zfblnonloc, zfnn, &
         zzd, zzn, zfc, zfv, zc1pbl, zturbqf, zturbtf, zturbuf, zturbvf, &
         zturbuvf, zmrk2, zsige
    real, pointer, dimension(:,:,:), contiguous :: zvcoef
    !---------------------------------------------------------------------

    call init2nan(bmsg, btsg, aq, bq, sfc_density, fsloflx, se, sig)
    call init2nan(kmsg, ktsg, wthl_ng, wqw_ng, uw_ng, vw_ng, c, d, unused)
    call init2nan(zero, qclocal, ficelocal, gam0)
    call init2nan(thl, qw, tthl, tqw, dkem, dket)

    ! Pointer assignments
    MKPTR1D(ps, pmoins, pvars)

    if (sfcflx_filter_order <= 0) then
       MKPTR1D(zalfaq, alfaq, pvars)
       MKPTR1D(zalfat, alfat, pvars)
       MKPTR1DK(zbt_ag, bt, indx_agrege, pvars)
    else
       if (kount == 0) then
          MKPTR1D(zalfaq0, alfaq, pvars)
          MKPTR1D(zalfat0, alfat, pvars)
          MKPTR1DK(zbt_ag0, bt, indx_agrege, pvars)
       endif
       MKPTR1D(zalfaq, falfaq, pvars)
       MKPTR1D(zalfat, falfat, pvars)
       MKPTR1DK(zbt_ag, fbt, indx_agrege, pvars)
    endif

    MKPTR1D(zbm, bm, pvars)

    MKPTR1D(zbm0, bm0, pvars)
    MKPTR1D(zfdsi, fdsi, pvars)
    MKPTR1D(zfdss, fdss, pvars)
    MKPTR1D(zflw, flw, pvars)
    MKPTR1D(zfq, fq, pvars)
    MKPTR1D(zfsh, fsh, pvars)
    MKPTR1D(zh, h, pvars)
    MKPTR1D(zmg, mg, pvars)
    MKPTR1D(ztsrad, tsrad, pvars)
    MKPTR1D(zue, ue, pvars)
    MKPTR1D(zustress, ustress, pvars)
    MKPTR1D(zvstress, vstress, pvars)

    MKPTR1DK(zfc_ag, fc, indx_agrege, pvars)
    MKPTR1DK(zfv_ag, fv, indx_agrege, pvars)
    MKPTR1DK(zfvap_ag, fvap, indx_agrege, pvars)
    MKPTR1DK(zqsurf_ag, qsurf, indx_agrege, pvars)

    MKPTR2DN(zfc, fc, ni, nagrege, pvars)
    MKPTR2DN(zfv, fv, ni, nagrege, pvars)
    MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, pvars)

    MKPTR2D(zturbqf, turbqf, pvars)
    MKPTR2D(zturbtf, turbtf, pvars)
    MKPTR2D(zturbuf, turbuf, pvars)
    MKPTR2D(zturbuvf, turbuvf, pvars)
    MKPTR2D(zturbvf, turbvf, pvars)

    MKPTR2Dm1(q, huplus, pvars)
    MKPTR2Dm1(zqcplus, qcplus, pvars)
    MKPTR2Dm1(sg, sigm, pvars)
    MKPTR2Dm1(zsige, sige, pvars)
    MKPTR2Dm1(zsigm, sigm, pvars)
    MKPTR2Dm1(zsigt, sigt, pvars)
    MKPTR2Dm1(zsigw, sigw, pvars)
    MKPTR2Dm1(t, tplus, pvars)
    MKPTR2Dm1(tm, tmoins, pvars)
    MKPTR2Dm1(tu, udifv, pvars)
    MKPTR2Dm1(tv, vdifv, pvars)
    MKPTR2Dm1(tt, tdifv, pvars)
    MKPTR2Dm1(tq, qdifv, pvars)
    MKPTR2Dm1(tqc, qcdifv, pvars)
    MKPTR2Dm1(uu, uplus, pvars)
    MKPTR2Dm1(vv, vplus, pvars)
    MKPTR2Dm1(w, wplus, pvars)
    MKPTR2Dm1(zc1pbl, c1pbl, pvars)
    MKPTR2Dm1(zfbl, fbl, pvars)
    MKPTR2Dm1(zfblgauss, fblgauss, pvars)
    MKPTR2Dm1(zfblnonloc, fblnonloc, pvars)
    MKPTR2Dm1(zfnn, fnn, pvars)
    MKPTR2Dm1(zgq, gq, pvars)
    MKPTR2Dm1(zgql, gql, pvars)
    MKPTR2Dm1(zgte, gte, pvars)
    MKPTR2Dm1(zgztherm, gztherm, pvars)
    MKPTR2Dm1(zkm, km, pvars)
    MKPTR2Dm1(zkt, kt, pvars)
    MKPTR2Dm1(zqtbl, qtbl, pvars)
    MKPTR2Dm1(ztve, tve, pvars)
    MKPTR2Dm1(zuwng, uwng, pvars)
    MKPTR2Dm1(zvwng, vwng, pvars)
    MKPTR2Dm1(zwqng, wqng, pvars)
    MKPTR2Dm1(zwtng, wtng, pvars)
    MKPTR2Dm1(zzd, zd, pvars)
    MKPTR2Dm1(zzn, zn, pvars)

    MKPTR2Dm1(tw, wdifv, pvars)
    MKPTR2Dm1(tl, qcdifv, pvars)
    MKPTR2Dm1(zsigmas, sigmas, pvars)
    MKPTR2Dm1(zpblq1, pblq1, pvars)

    MKPTR3D(zvcoef, vcoef, pvars)

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
       sigef(1:,1:) => zsige(1:ni,1:nkm1)
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
          gam0(j,k) = rsg*zsige(j,k)/ztve(j,k)
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
    se = zsige(1,:)
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
       !#TODO: check - Looks like zbm0 was never init in the code (only =0 at bus init)...
       bmsg(:) = fm_mult(:) * ((zbm(:)*WEIGHT_BM + zbm0(:)*(1.-WEIGHT_BM) )*aq(:))
       zbm0(:) = zbm(:)*WEIGHT_BM + zbm0(:)*(1.-WEIGHT_BM)
    endif
    gam0 = 0.

    ! ** U-Component Wind Diffusion **

    ! Compute tendency for u-component wind diffusion
    if (compute_tend) then
       call difuvdfj(tu,uu,kmsg,zero,uw_ng,zero,zero,bmsg,sg,zsige,tau, &
            typem,1.,ni,ni,ni,nkm1)
       if (phy_error_L) return
    endif

    ! Diagnose atmospheric fluxes for u-component wind
    if (ISREQSTEPL((/'TFUU', 'TFUV'/))) then
       call atmflux4(zturbuf, tu, zsigm, ps, ni, nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! ** V-Component Wind Diffusion **

    ! Compute tendency for v-component wind diffusion
    if (compute_tend) then
       call difuvdfj(tv,vv,kmsg,zero,vw_ng,zero,zero,bmsg,sg,zsige,tau, &
            typem,1.,ni,ni,ni,nkm1)
       if (phy_error_L) return
    endif

    ! Diagnose atmospheric fluxes for v-component wind
    if (ISREQSTEPL((/'TFVV', 'TFUV'/))) then
       call atmflux4(zturbvf, tv, zsigm, ps, ni, nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Diagnose total turbulent momentum flux
    if (ISREQSTEP('TFUV')) &
         zturbuvf(:,:) = sqrt(zturbuf(:,:)**2+zturbvf(:,:)**2)

    ! ** Vertical Motion Diffusion **

    ! Compute tendency for vertical motion diffusion on user request
    if (diffuw .and. compute_tend) then
       call difuvdfj(tw,w,kmsg,zero,zero,zero,zero,zero,sg,sigef,tau, &
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
       call difuvdfj(tconserv_q,conserv_q,ktsg,zgq,wqw_ng,zero,aq,bq,sg,sigef,tau, &
            typet,1.,ni,ni,ni,nkm1)
       if (phy_error_L) return
    endif

    ! Diagnose atmospheric fluxes for moisture
    if (ISREQSTEP('TFHU')) then
       call atmflux4(zturbqf, tconserv_q, zsigt, ps, ni, nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! ** Temperature Diffusion **

    ! Set surface flux boundary terms
    aq(:) = fh_mult(:) * zalfat(:)
    bq(:) = fh_mult(:) * btsg(:)

    ! Compute conserved temperature variable for vertical diffusion
    if (fluvert == 'MOISTKE') then
       ! Liquid water potential temperature (thl)
       call ficemxp(ficelocal,unused,unused,tm,ni,nkm1)
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
       call difuvdfj(tconserv_t,conserv_t,ktsg,zgte,wthl_ng,zero,aq,bq,sg,sigef,tau, &
            typet,1.,ni,ni,ni,nkm1)
       if (phy_error_L) return
    endif

    ! Convert conserved variable tendencies to state (T,Q) variable tendencies
    if (fluvert == 'MOISTKE') then
       call baktotq(tt,tq,tl,thl,qw,tthl,tqw,qclocal,sg,zsigw, &
            ps,zgztherm,tm,ficelocal,ztve,zh,zfc_ag,zqtbl,zfnn,zfbl,zfblgauss, &
            zfblnonloc,zc1pbl,zzn,zzd,zmg,zmrk2,zvcoef,zsigmas,zpblq1,tau,ni,nkm1)
    endif

    ! Compute tendency for condensate diffusion
    if (pbl_diff_condens .and. compute_tend) then
       aq = 0.
       bq = 0.
       call difuvdfj(tqc,zqcplus,ktsg,zgq,wqw_ng,zero,aq,bq,sg,sigef,tau, &
            typet,1.,ni,ni,ni,nkm1)
       if (phy_error_L) return
    else
       tqc = 0.
    endif

    ! Dissipative heating
    call dissheat(dket, uu, vv, tu, tv, kmsg, zsigm, zsigt, zvcoef, tau, ni, nkm1)
    tt(:,1:nkm1) = tt(:,1:nkm1) - (1./CPD) * dket(:,1:nkm1)

    ! Diagnose atmospheric fluxes for temperature
    if (ISREQSTEP('TFTT')) then
       call atmflux4(zturbtf, tt, zsigt, ps, ni, nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Diagnose (implicit) surface fluxes based on updated state
    call sfcflux(zustress, zvstress, zfq, zue, zfsh, zflw, &
         uu, vv, t, q, tu, tv, tt, tq, bmsg, btsg, zalfaq, &
         sfc_density, ps, tau, ni, nkm1)
    if (phy_error_L) return

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
  end subroutine difver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DIFUVDFj(TU, U, KU, GU, JNG, R, ALFA, BETA, S, SK, &
       TAU, itype, F, NU, NR, N, NK)
    use, intrinsic :: iso_fortran_env, only: REAL64
    implicit none
!!!#include <arch_specific.hf>

    integer NU, NR, N, NK
    real TU(NU, NK), U(NU, NK), KU(NR, NK), GU(NR, NK), R(NR,NK)
    real JNG(NR, NK)
    real ALFA(N), BETA(N), S(n,NK), SK(n,NK), TAU, F
    integer itype

    !@Author R. Benoit (Mar 89)
    !@Revisions
    ! 001      R. Benoit (Aug 93) -Local sigma: s and sk become 2D
    ! 002      B. Bilodeau (Dec 94) - "IF" tests on integer
    !          instead of character.
    ! 003      J. Mailhot (Sept 00) - Add itype 4='EB'
    ! 004      A. PLante (June 2003) - IBM conversion
    !             - calls to vrec routine (from massvp4 library)
    ! 005      J. Mailhot/L. Spacek (Dec 07) - Add itype 5='ET' and cleanup
    ! 006      J. Mailhot (Aug 12) - Add argument and modifications (non-gradient flux,
    !                                geometric averaging) and change name to difuvdfj1
    !@Object
    !          to solve a vertical diffusion equation by finite
    !          differences
    !@Arguments
    !          - Output -
    ! TU       U tendency (D/DT U) due to the vertical diffusion and to
    !          term R
    !          - Input -
    ! U        variable to diffuse (U,V,T,Q,E)
    ! KU       diffusion coefficient
    ! GU       optional countergradient term
    ! JNG      optional non-gradient flux term
    ! R        optional inhomogeneous term
    ! ALFA     inhomogeneous term for the surface flux (for itype 1='U', 2='UT' or 5='ET'or 6='ST')
    !          surface boundary condition (for itype 4='EB')
    ! BETA     homogeneous term for the surface flux
    ! S        sigma coordinates of full levels
    ! SK       sigma coordinates of diffusion coefficient
    !          (or staggered variables) levels
    ! TAU      length of timestep
    ! ITYPE    itype of variable to diffuse (1='U',2='UT',3='E',4='EB' or 5='ET' or 6='ST')
    ! F        weighting factor for time 'N+1'
    ! NU       1st dimension of TU and U
    ! NR       1st dimension of KU, GU, JNG and R
    ! N        number of columns to process
    ! NK       vertical dimension
    !@Notes
    !          D/DT U = D(U) + R
    !          D(U) = D/DS ( J(U) + JNG )
    !          J(U) = KU*(D/DS U + GU)
    !          Limiting Conditions where S=ST: J=0(for 'U'/'ET'), D=0(for 'UT'
    !          and ST=1)
    !**        U=0(for 'E' and 'EB')
    !          J=0(for 'E' and 'EB')
    !          Limiting Conditions where S=SB: J=ALFA+BETA*U(S(NK))(for
    !          'U'/'UT'/'ET'), J=0(for 'E'), U=ALFA(for 'EB')
    !          ST = S(1)-1/2 (S(2)-S(1)) (except for 'TU')
    !          SB = SK(NK) = 1.
    !*@/

    integer I, K, NKX
    real ST, SB, HM, HP, KUM, KUP, SCK1
    real, dimension(N) :: HD
    real, dimension(N,NK) :: VHM,VHP,A,B,C,D
    real(KIND=8), dimension(N,NK) :: RHD,RHMD,RHPD
    logical :: SFCFLUX
    character(len=16) :: msg_S
    external DIFUVD1, DIFUVD2

    st(i)=s(i,1)-0.5*(s(i,2)-s(i,1))
    sb(i)=1.

    if (itype.le.2) then
       NKX=NK
       SCK1=1
       if (itype.eq.2) then
          SCK1=0
       endif
    else if (itype.eq.3 .or. itype.eq.4) then
       NKX=NK-1
    else if (itype.eq.5 .or. itype.eq.6) then
       NKX=NK
    else
       write(msg_S, '(i0)') itype
       call physeterror('difuvdfj', 'Type inconnu: '//trim(msg_S))
       return
    endif

    ! (1) CONSTRUIRE L'OPERATEUR TRIDIAGONAL DE DIFFUSION N=(A,B,C)
    !                ET LE TERME CONTRE-GRADIENT (DANS D)

    if (itype.le.2) then

       !     K=1

       HM=0
       do I=1,N
          HP=S(i,2)-S(i,1)
          HD(I)=SK(i,1)-ST(i)
          A(I,1)=0
          C(I,1)=SCK1*KU(I,1)/(HP*HD(I))
          B(I,1)=-A(I,1)-C(I,1)
          D(I,1)=SCK1*(KU(I,1)*GU(I,1)+JNG(I,1))/HD(I)
       enddo

       !     K=2...NK-1

       do K=2,NK-1,1
          do I=1,N
             !              THE FOLLOWING LHS ARE IN REAL
             VHM(I,K)=S(I,K)-S(I,K-1)
             VHP(I,K)=S(I,K+1)-S(I,K)
             HD(I)=SK(I,K)-SK(I,K-1)
             !           THE FOLLOWING LHS ARE IN real(REAL64)
             RHD(I,K)=1.d0/HD(I)
             RHMD(I,K)=1.d0/(VHM(I,K)*HD(I))
             RHPD(I,K)=1.d0/(VHP(I,K)*HD(I))

             A(I,K)=KU(I,K-1)*RHMD(I,K)
             C(I,K)=KU(I,K)*RHPD(I,K)
             B(I,K)=-A(I,K)-C(I,K)
             D(I,K)=( KU(I,K)*GU(I,K)-KU(I,K-1)*GU(I,K-1) &
                  +JNG(I,K)-JNG(I,K-1) )*RHD(I,K)
          enddo
       enddo

       !     K=NK

       HP=0
       do I=1,N
          HM=S(i,NK)-S(i,NK-1)
          HD(I)=SB(i)-SK(i,NK-1)
          A(I,NK)=KU(I,NK-1)/(HM*HD(I))
          C(I,NK)=0
          B(I,NK)=-A(I,NK)-C(I,NK)
          D(I,NK)=(0-KU(I,NK-1)*GU(I,NK-1)-JNG(I,NK-1))/HD(I)
       enddo

    else if (itype.eq.3 .or. itype.eq.4 .or. itype.eq.5 .or. itype.eq.6) then

       !     ITYPE='E' or 'EB' or 'ET'

       !     K=1

       do I=1,N
          HM=SK(i,1)-ST(i)
          HP=SK(i,2)-SK(i,1)
          HD(I)=S(i,2)-S(i,1)
          !        Limiting Conditions at S=ST: U=0(for 'E' or 'EB')`
          !**         KUM=0.5*KU(I,1)
          !        Limiting Conditions at S=S(1): J=0(for 'E' or 'EB' or 'ET')
          KUM=0
          KUP=0.5*(KU(I,1)+KU(I,2))
          A(I,1)=KUM/(HM*HD(I))
          C(I,1)=KUP/(HP*HD(I))
          B(I,1)=-A(I,1)-C(I,1)
          D(I,1)=(KUP*(GU(I,1)+GU(I,2))-KUM*GU(I,1) &
               +(JNG(I,1)+JNG(I,2)) )/(2.*HD(I))
       enddo

       !     K=2...NKX-1

       do K=2,NKX-1,1
          do I=1,N
             !              THE FOLLOWING LHS ARE IN REAL
             VHM(I,K)=SK(I,K)-SK(I,K-1)
             VHP(I,K)=SK(I,K+1)-SK(I,K)
             HD(I)=S(I,K+1)-S(I,K)
          enddo
          if (K==NK-1 .and. itype==6) then !VIRTUAL LEVEL FOR ITYPE='ST'
             do I=1,N
                HD(I) = 0.5*(SK(I,K+1)+SK(I,K))-S(I,K)
             enddo
          endif
          do I=1,N
             RHD(I,K)=1.d0/HD(I)
             RHMD(I,K)=1.d0/(VHM(I,K)*HD(I))
             RHPD(I,K)=1.d0/(VHP(I,K)*HD(I))

             KUM=0.5*(KU(I,K-1)+KU(I,K))
             KUP=0.5*(KU(I,K+1)+KU(I,K))
             A(I,K)=KUM*RHMD(I,K)
             C(I,K)=KUP*RHPD(I,K)
             B(I,K)=-A(I,K)-C(I,K)
             D(I,K)=.5*(KUP*(GU(I,K)+GU(I,K+1)) &
                  -KUM*(GU(I,K-1)+GU(I,K)) &
                  +(JNG(I,K)+JNG(I,K+1)) &
                  -(JNG(I,K-1)+JNG(I,K)))*RHD(I,K)
          enddo
       enddo

       !     K=NKX

       if (itype.eq.3 .or. itype.eq.5 .or. itype.eq.6) then

          !       ITYPE='E' or 'ET' or 'ST'

          if (itype.eq.6) then !virtual level for ITYPE='ST'
             do I=1,N
                HD(I)=SB(i)-0.5*(SK(i,NKX)+SK(i,NKX-1))
             enddo
          else
             do I=1,N
                HD(I)=SB(i)-S(i,NKX)
             enddo
          endif
          do I=1,N
             HM=SK(i,NKX)-SK(i,NKX-1)
             KUM=0.5*(KU(I,NKX)+KU(I,NKX-1))
             KUP=0
             A(I,NKX)=KUM/(HM*HD(I))
             C(I,NKX)=0
             B(I,NKX)=-A(I,NKX)-C(I,NKX)
             D(I,NKX)=(0-KUM*(GU(I,NKX)+GU(I,NKX-1)) &
                  -(JNG(I,NKX)+JNG(I,NKX-1)) )/(2.*HD(I))
          enddo

       else if (itype.eq.4) then

          !       ITYPE='EB'

          do I=1,N
             HM=SK(i,NK-1)-SK(i,NK-2)
             HP=SB(i)-SK(i,NK-1)
             HD(I)=S(i,NK)-S(i,NK-1)
             KUM=0.5*(KU(I,NK-1)+KU(I,NK-2))
             KUP=0.5*(KU(I,NK)+KU(I,NK-1))
             A(I,NKX)=KUM/(HM*HD(I))
             B(I,NKX)=-A(I,NKX) -KUP/(HP*HD(I))
             C(I,NKX)=0
             D(I,NKX)=(KUP*(GU(I,NK)+GU(I,NK-1)) &
                  -KUM*(GU(I,NK-1)+GU(I,NK-2)) &
                  +(JNG(I,NK)+JNG(I,NK-1)) &
                  -(JNG(I,NK-1)*JNG(I,NK-2)))/(2.*HD(I)) &
                  +KUP*ALFA(I)/(HD(I)*HP)
          enddo

       endif

    endif


    ! (2) CALCULER LE COTE DROIT D=TAU*(N(U)+R+D/DS(KU*GU+JNG))

    call DIFUVD1(D, 1., A, B, C, U, D, N, NU, NKX)
    do K=1,NKX
       do I=1,N
          D(I,K)=TAU*(D(I,K)+R(I,K))
       enddo
    enddo

    ! (3) CALCULER OPERATEUR DU COTE GAUCHE

    do K=1,NKX
       do I=1,N
          A(I,K)= -F*TAU*A(I,K)
          B(I,K)=1-F*TAU*B(I,K)
          C(I,K)= -F*TAU*C(I,K)
       enddo
    enddo

    ! (4) AJOUTER TERME DE FLUX DE SURFACE

    SFCFLUX = .true.
    select case (itype)
    case (:2) !ITYPE='U'/'UT'
       do I=1,N
          HD(I) = SB(i) - SK(i,NK-1)
       enddo
    case (5) !ITYPE='ET'
       do I=1,N
          HD(I) = SB(i) - S(i,NKX)
       enddo
    case (6) !ITYPE='ST'
       do I=1,N
          HD(I) = SB(i) - 0.5*(SK(i,NKX)+SK(i,NKX-1))
       enddo
    case DEFAULT
       SFCFLUX = .false.
    end select
    if (SFCFLUX) then
       do I=1,N
          B(I,NKX)=B(I,NKX)-TAU*BETA(I)/HD(I)
          D(I,NKX)=D(I,NKX)+(ALFA(I)+BETA(I)*U(I,NKX))*TAU/HD(I)
       enddo
    endif

    ! (5) RESOUDRE SYSTEME TRIDIAGONAL [A,B,C] X = D. METTRE X DANS TU.

    call DIFUVD2(TU, A, B, C, D, D, NU, N, NKX)

    ! (6) OBTENIR TENDANCE

    do K=1,NKX
       do I=1,N
          TU(I,K)=TU(I,K)/TAU
       enddo
    enddo
    !     K=NKX+1..NK
    do K=NKX+1,NK
       do I=1,N
          TU(I,K)=0
       enddo
    enddo

    return
  end subroutine DIFUVDFj
  
end module pbl_diffuse
