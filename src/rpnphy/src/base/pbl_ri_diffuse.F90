module pbl_ri_diffuse
  implicit none
  private

  ! API constants
  integer, parameter, public :: FIELD_STAGGERED_FROM_ELEV=1
  integer, parameter, public :: FIELD_COLLOCATED_WITH_ELEV=2
  integer, parameter, public :: FIELD_ON_ELEV=3
  
  ! API functions
  public :: diffuse
  public :: diffuseall

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine diffuseall(F_pvars, F_tau, F_ni, F_nk, F_nkm1)
    use debug_mod, only: init2nan
    use tdpack_const, only: CAPPA, CHLC, CHLF, CPD, DELTA, GRAV, RGASD
    use phy_options, except=>pbl_dissheat
    use phy_status, only: phy_error_L
    use phybusidx
    use phymem, only: phyvar
    use ens_perturb, only: ens_nc2d, ens_spp_get
    use atmflux, only: atmflux4
    use pbl_ri_utils, only: baktotq, compute_conserved
    use pbl_utils, only: dissheat, sfcflux
    use integrals, only: INT_TYPE_LINEAR

    implicit none

    ! Argument declarations
    integer, intent(in) :: F_ni                         !Horizontal dimension
    integer, intent(in) :: F_nk                         !Vertical dimension
    integer, intent(in) :: F_nkm1                       !Vertical dimension of operator
    type(phyvar), pointer, contiguous :: F_pvars(:)     !Physics buses
    real, intent(in) :: F_tau                           !Time step (s)

    !@Object   Computed turbulent transport (diffusion) tendencies

    ! Common declarations
#include "phymkptr.hf"
    include "surface.cdk"
    
    ! Local parameter definitions
    integer, parameter :: FLUX_INTTYPE=INT_TYPE_LINEAR

    ! Local declarations
    integer :: i, k, ni, nk, nkm1
    real :: rsg, gam0
    real, dimension(F_ni) :: bmsg, btsg, aq, bq, rhosfc, fm_mult, fh_mult
    real, dimension(F_ni,F_nkm1) ::  kmsg, ktsg, zero, dket

    ! Bus pointer declarations
    real, pointer, dimension(:), contiguous :: zbt_ag, zfv_ag, zqsurf_ag, zfvap_ag, &
         ps, zalfat, zalfaq, zbm, zbm0, zfdsi, zfdss, zfq, ztsrad, zustress, &
         zvstress, zue, zflw, zfsh
    real, pointer, dimension(:,:), contiguous :: tu, tv, tw, tt, tq, tl, uu, vv, w, &
         t, q, zsigt, zsigm, tm, zpblsigs,zpblq1, zqcplus, tqc, zze, zkm, zkt, zqtbl, &
         ztve, zfbl, zzd, zzn, zfc, zfv, zpri, zturbqf, zturbtf, zturbuf, zturbvf, &
         zturbuvf, zmrk2, ztcons, zqcons, zdtcons, zdqcons, zturbtlf, zsige
    real, pointer, dimension(:,:,:), contiguous :: zvcoef

    ! Initialize local arrays to detect undefined fields
    call init2nan(bmsg, btsg, aq, bq, rhosfc, fm_mult, fh_mult)
    call init2nan(kmsg, ktsg, zero, dket)

    ! Initialize array bounds for bus pointer macros
    ni = F_ni; nk = F_nk ; nkm1 = F_nkm1

    ! Pointer assignments
    MKPTR1D(ps, pmoins, F_pvars)
    MKPTR1D(zalfaq, alfaq, F_pvars)
    MKPTR1D(zalfat, alfat, F_pvars)
    MKPTR1D(zbm, bm, F_pvars)
    MKPTR1D(zbm0, bm0, F_pvars)
    MKPTR1D(zfdsi, fdsi, F_pvars)
    MKPTR1D(zfdss, fdss, F_pvars)
    MKPTR1D(zflw, flw, F_pvars)
    MKPTR1D(zfq, fq, F_pvars)
    MKPTR1D(zfsh, fsh, F_pvars)
    MKPTR1D(ztsrad, tsrad, F_pvars)
    MKPTR1D(zue, ue, F_pvars)
    MKPTR1D(zustress, ustress, F_pvars)
    MKPTR1D(zvstress, vstress, F_pvars)

    MKPTR1DK(zbt_ag, bt, indx_agrege, F_pvars)
    MKPTR1DK(zfv_ag, fv, indx_agrege, F_pvars)
    MKPTR1DK(zfvap_ag, fvap, indx_agrege, F_pvars)
    MKPTR1DK(zqsurf_ag, qsurf, indx_agrege, F_pvars)

    MKPTR2DN(zfc, fc, ni, nagrege, F_pvars)
    MKPTR2DN(zfv, fv, ni, nagrege, F_pvars)
    MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, F_pvars)

    MKPTR2D(zturbqf, turbqf, F_pvars)
    MKPTR2D(zturbtf, turbtf, F_pvars)
    MKPTR2D(zturbtlf, turbtlf, F_pvars)
    MKPTR2D(zturbuf, turbuf, F_pvars)
    MKPTR2D(zturbuvf, turbuvf, F_pvars)
    MKPTR2D(zturbvf, turbvf, F_pvars)

    MKPTR2Dm1(q, huplus, F_pvars)
    MKPTR2Dm1(zqcplus, qcplus, F_pvars)
    MKPTR2Dm1(zsige, sige, F_pvars)
    MKPTR2Dm1(zsigt, sigt, F_pvars)
    MKPTR2Dm1(zsigm, sigm, F_pvars)
    MKPTR2Dm1(t, tplus, F_pvars)
    MKPTR2Dm1(tm, tmoins, F_pvars)
    MKPTR2Dm1(tu, udifv, F_pvars)
    MKPTR2Dm1(tv, vdifv, F_pvars)
    MKPTR2Dm1(tt, tdifv, F_pvars)
    MKPTR2Dm1(tq, qdifv, F_pvars)
    MKPTR2Dm1(tqc, qcdifv, F_pvars)
    MKPTR2Dm1(uu, uplus, F_pvars)
    MKPTR2Dm1(vv, vplus, F_pvars)
    MKPTR2Dm1(w, wplus, F_pvars)
    MKPTR2Dm1(zdqcons, dqcons, F_pvars)
    MKPTR2Dm1(zdtcons, dtcons, F_pvars)
    MKPTR2Dm1(zfbl, fbl, F_pvars)
    MKPTR2Dm1(zkm, km, F_pvars)
    MKPTR2Dm1(zkt, kt, F_pvars)
    MKPTR2Dm1(zpri, pri, F_pvars)
    MKPTR2Dm1(zqcons, qcons, F_pvars)
    MKPTR2Dm1(ztcons, tcons, F_pvars)
    MKPTR2Dm1(ztve, tve, F_pvars)
    MKPTR2Dm1(zzd, zd, F_pvars)
    MKPTR2Dm1(zze, ze, F_pvars)
    MKPTR2Dm1(zzn, zn, F_pvars)
    MKPTR2Dm1(tw, wdifv, F_pvars)
    MKPTR2Dm1(tl, ldifv, F_pvars)
    MKPTR2Dm1(zpblsigs, pblsigs, F_pvars)
    MKPTR2Dm1(zpblq1, pblq1, F_pvars)
    
    MKPTR3D(zvcoef, vcoef, F_pvars)

    if (advecqtbl) then
       MKPTR2Dm1(zqtbl, qtblplus, F_pvars)
    else
       MKPTR2Dm1(zqtbl, qtbl, F_pvars)
    endif

    ! Retrieve stochastic parameter information on request
    fm_mult(:) = ens_spp_get('fm_mult', zmrk2, default=1.)
    fh_mult(:) = ens_spp_get('fh_mult', zmrk2, default=1.)

    ! Conversion precalculations for sigma coordinates
    zero = 0.
    rsg = (GRAV/RGASD)
    do k=1,F_nkm1
       do i=1,F_ni
          gam0 = (rsg*zsige(i,k)/ztve(i,k))**2
          kmsg(i,k) = zkm(i,k)*gam0
          ktsg(i,k) = zkt(i,k)*gam0
       end do
    end do

    ! Compute turbulent flux surface boundary conditions
    do i=1,F_ni
       aq(i)=-rsg/(t(i,F_nkm1)*(1. + DELTA * zqsurf_ag(i)))
       bmsg(i)   = fm_mult(i) * zbm(i) * aq(i)
       btsg(i)   = zbt_ag(i) * aq(i)
       zalfat(i) = zalfat(i)  * aq(i)
       zalfaq(i) = zalfaq(i)  * aq(i)
       rhosfc(i) = -aq(i) * ps(i)/GRAV
    end do

    ! Compute conserved thermodynamic variables for vertical diffusion
    call compute_conserved(ztcons, zqcons, t, q, zqtbl, zsigt, F_ni, F_nkm1)
    
    ! Compute tendency for u-component wind and diagnose flux on request
    call diffuse(tu, uu, kmsg, zero, bmsg, zsige, zsigt, F_tau, &
         FIELD_COLLOCATED_WITH_ELEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any([(any((/'tfuu', 'tfuv'/) == phyoutlist_S(i)), i=1,nphyoutlist)])) then
       call atmflux4(zturbuf, tu, zsigm, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Compute tendency for v-component wind and diagnose flux on request
    call diffuse(tv, vv, kmsg, zero, bmsg, zsige, zsigt, F_tau, &
         FIELD_COLLOCATED_WITH_ELEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any([(any((/'tfvv', 'tfuv'/) == phyoutlist_S(i)), i=1,nphyoutlist)])) then
       call atmflux4(zturbvf, tv, zsigm, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Diagnose total turbulent momentum flux
    if (any('tfuv' == phyoutlist_S)) &
         zturbuvf(:,:) = sqrt(zturbuf(:,:)**2+zturbvf(:,:)**2)

    ! Set surface flux boundary terms for moisture diffusion
    aq(:) = fh_mult(:) * zalfaq(:)
    bq(:) = fh_mult(:) * btsg(:)

    ! Compute tendency for moisture and diagnose flux on request
    call diffuse(zdqcons, zqcons, ktsg, aq, bq, zsige, zsigt, F_tau, &
         FIELD_STAGGERED_FROM_ELEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any(phyoutlist_S == 'tfhu')) then
       call atmflux4(zturbqf, zdqcons, zsigt, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif
    
    ! Set surface flux boundary terms for thermal diffusion
    aq(:) = fh_mult(:) * zalfat(:)
    bq(:) = fh_mult(:) * btsg(:)

    ! Compute tendency for theta_l and diagnose flux on request 
    call diffuse(zdtcons, ztcons, ktsg, aq, bq, zsige, zsigt, F_tau, &
         FIELD_STAGGERED_FROM_ELEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any(phyoutlist_S == 'tftl')) then
       call atmflux4(zturbtlf, zdtcons, zsigt, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Convert conserved variable tendencies to state (T, Q) variable tendencies
    call baktotq(tt, tq, tl, ztcons, zqcons, zdtcons, zdqcons, zqtbl, zsigt, &
         ps, tm, ztve, zfbl, zpri, zzn, zzd, zvcoef, zpblsigs, zpblq1, &
         F_tau, F_ni, F_nkm1)

    ! Update PBL cloud condensate
    zqtbl = max(zqtbl, 0.) + tl*F_tau
    
    ! Compute tendency for condensate diffusion  FIXME (replace with zqcdifv diagnosis)
    tqc = 0.

    ! Dissipative heating
    call dissheat(dket, uu, vv, tu, tv, kmsg, zsigm, zsigt, zvcoef, F_tau, F_ni, F_nkm1)
    tt(:,1:F_nkm1) = tt(:,1:F_nkm1) - (1./CPD) * dket(:,1:F_nkm1)

    ! Diagnose atmospheric fluxes for dry-air temperature
    if (any(phyoutlist_S == 'tftt')) then
       call atmflux4(zturbtf, tt, zsigt, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif
    
    ! Diagnose (implicit) surface fluxes based on updated atmospheric state
    if (any([(any((/'f3 ','f2 ','fq ','ue ','fsh','flw'/) == &
         phyoutlist_S(i)), i=1,nphyoutlist)])) then
       call sfcflux(zustress, zvstress, zfq, zue, zfsh, zflw, &
            uu, vv, t, q, tu, tv, tt, tq, bmsg, btsg, zalfaq, &
            rhosfc, ps, F_tau, F_ni, F_nkm1)
       if (phy_error_L) return
    endif

    ! End of subprogram
    return
  end subroutine diffuseall
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine diffuse(F_tend, F_fld, F_kcoef, F_alpha, F_beta, F_sige, F_sigt, &
       F_tau, F_type, F_ni, F_nk)
    use, intrinsic :: iso_fortran_env, only: REAL64
    use debug_mod, only: init2nan
    use tridiag, only: td_solve
    implicit none

    ! Argument description
    integer, intent(in) :: F_ni                 !horizontal dimension
    integer, intent(in) :: F_nk                 !vertical dimension
    real, dimension(F_ni,F_nk) :: F_fld         !field to diffuse
    real, dimension(F_ni,F_nk) :: F_kcoef       !diffusion coefficients on e-levs
    real, dimension(F_ni) :: F_alpha            !field-independent sfc flux term (inhomogeneous)
    real, dimension(F_ni) :: F_beta             !field-dependent sfc flux term (homogeneous)
    real, dimension(F_ni,F_nk) :: F_sige        !sigma values of energy levels
    real, dimension(F_ni,F_nk) :: F_sigt        !sigma values of thermo levels (e-layer interfaces)
    real, intent(in) :: F_tau                   !time step (s)
    integer, intent(in) :: F_type               !type of field to diffuse (use API constants)
    real, dimension(F_ni,F_nk) :: F_tend        !field tendency due to vertical diffusion

    !@Object:  Solve a vertical diffusion equation

    ! Internal variable declarations
    integer :: i, k
    real :: hm, hp, kum, kup, hk
    real, dimension(F_ni,F_nk) :: a, b, c, d, fldupd

    ! Initialize local arrays to detect undefined fields
    call init2nan(a, b, c, d, fldupd)
    
    ! Construct the tridiagonal matrix (N=[A,B,C]) and forcing (D)
    STAGGERING: select case(F_type)

    case (FIELD_STAGGERED_FROM_ELEV)
       ! Field is centered between energy levels

       ! Top level
       do i=1,F_ni
          hp = F_sigt(i,1) - F_sigt(i,2)
          hk = F_sige(i,1) - F_sige(i,2)
          a(i,1) = 0.
          c(i,1) = -F_tau * F_kcoef(i,2) / (hp*hk)
          b(i,1) = 1. - (a(i,1) + c(i,1))
          d(i,1) = F_fld(i,1)
       enddo

       ! Interior of the column
       do k=2,F_nk-1
          do i=1,F_ni
             hm = F_sigt(i,k-1) - F_sigt(i,k)
             hp = F_sigt(i,k) - F_sigt(i,k+1)
             hk = F_sige(i,k) - F_sige(i,k+1)
             a(i,k) = -F_tau * F_kcoef(i,k) / (hm*hk)
             c(i,k) = -F_tau * F_kcoef(i,k+1) / (hp*hk)
             b(i,k) = 1. - (a(i,k) + c(i,k))
             d(i,k) = F_fld(i,k)    
          enddo
       enddo

       ! Bottom level
       do i=1,F_ni
          hm = F_sigt(i,F_nk-1) - F_sigt(i,F_nk)
          hk = F_sige(i,F_nk) - 1.
          a(i,F_nk) = -F_tau * F_kcoef(i,F_nk)/(hm*hk)
          c(i,F_nk) = 0.
          b(i,F_nk) = 1. - (a(i,F_nk) - F_tau * F_beta(i)/hk)
          d(i,F_nk) = F_fld(i,F_nk) - F_tau * F_alpha(i)/hk
       enddo
       
    case (FIELD_COLLOCATED_WITH_ELEV, FIELD_ON_ELEV)
       ! Field is colocated with energy levels

       ! Top level
       do i=1,F_ni
          hp = F_sige(i,1) - F_sige(i,2)
          hk = 2. * (F_sige(i,1) - F_sigt(i,1))
          kup = 0.5*(F_kcoef(i,1)+F_kcoef(i,2))
          a(i,1) = 0.          
          c(i,1) = -F_tau * kup / (hp*hk)
          b(i,1) = 1. - (a(i,1) + c(i,1))
          d(i,1) = F_fld(i,1)
       enddo

       ! Interior of the column 
       do k=2,F_nk-1
          do i=1,F_ni
             hm = F_sige(i,k-1) - F_sige(i,k)
             hp = F_sige(i,k) - F_sige(i,k+1)
             hk = F_sigt(i,k-1) - F_sigt(i,k)
             kum = 0.5 * (F_kcoef(i,k-1) + F_kcoef(i,k))
             kup = 0.5 * (F_kcoef(i,k+1) + F_kcoef(i,k))
             a(i,k) = -F_tau * kum / (hm*hk)
             c(i,k) = -F_tau * kup / (hp*hk)
             b(i,k) = 1. - (a(i,k) + c(i,k))
             d(i,k) = F_fld(i,k)
          enddo
       enddo

       ! Bottom level
       BOTTOM_ELEV: if (F_type == FIELD_ON_ELEV) then
          ! Bottom level is an energy level (surface)
          k = F_nk
          do i=1,F_ni
             hm = F_sige(i,k-1) - F_sige(i,k)
             hp = F_sige(i,k) - 1.
             hk = F_sigt(i,k-1) - F_sigt(i,k)
             kum = 2.*F_kcoef(i,k-1)*F_kcoef(i,k) / (F_kcoef(i,k-1) + F_kcoef(i,k))             
             kup = F_kcoef(i,k) / (2.*sqrt(2.)) * sqrt(1. + F_alpha(i)/(F_fld(i,k)))
             a(i,k) = -F_tau * kum / (hm*hk)
             c(i,k) = 0.
             b(i,k) = 1. - (a(i,k) - F_tau * kup / (hp*hk))
             d(i,k) = F_fld(i,k) + F_tau * kup * F_alpha(i) / (hp*hk)
          enddo   
       else
          ! Bottom level is above-ground
          do i=1,F_ni
             hk = F_sigt(i,F_nk-1) - 1.
             hm = F_sige(i,F_nk-1) - F_sige(i,F_nk)
             kum = 0.5*(F_kcoef(i,F_nk) + F_kcoef(i,F_nk-1))
             a(i,F_nk) = -F_tau * kum/(hm*hk)
             c(i,F_nk) = 0.
             b(i,F_nk) = 1. - (a(i,F_nk) - F_tau * F_beta(i)/hk)
             d(i,F_nk) = F_fld(i,F_nk) - F_tau * F_alpha(i)/hk
          enddo
       endif BOTTOM_ELEV

    case DEFAULT
       call physeterror('pbl_ri_diffuse::diffuse', 'Invalid staggering (F_type) selected')
       return

    end select STAGGERING
    
    ! Solve the tridiagonal system for the updated field value
    call td_solve(fldupd, a, b, c, d, F_ni, F_nk)

    ! Retrieve tendency
    F_tend = (fldupd - F_fld) / F_tau

    return
  end subroutine diffuse

end module pbl_ri_diffuse
