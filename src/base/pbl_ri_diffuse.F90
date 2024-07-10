module pbl_ri_diffuse
  implicit none
  private

  ! API constants
  integer, parameter, public :: FIELD_ON_MLEV=1
  integer, parameter, public :: FIELD_ON_TLEV=2
  integer, parameter, public :: FIELD_IS_TKE=3
  
  ! API functions
  public :: diffuse
  public :: diffuseall

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine diffuseall(F_pvars, F_tau, F_ni, F_nk, F_nkm1)
    use debug_mod, only: init2nan
    use tdpack_const, only: CAPPA, CHLC, CHLF, CPD, DELTA, GRAV, RGASD
    use vintphy, only: vint_thermo2mom
    use phy_options, except=>pbl_dissheat
    use phy_status, only: phy_error_L, PHY_OK
    use phybusidx, except1=>lwc, except2=>iwc, except3=>ti
    use phymem, only: phyvar
    use ens_perturb, only: ens_nc2d, ens_spp_get
    use microphy_utils, only: mp_lwc, mp_iwc
    use atmflux, only: atmflux4
    use pbl_ri_utils, only: baktot, compute_conserved, ktosig, sgsvar, K_M, K_T
    use pbl_utils, only: dissheat, sfcflux
    use integrals, only: INT_TYPE_LINEAR
    use tendency, only: apply_tendencies

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
    real, dimension(F_ni) :: bmsg, btsg, aq, bq, rhosfc, fm_mult, fh_mult, zero
    real, dimension(F_ni,F_nkm1) ::  kmsg, ktsg, dket, tvt, tvm, thl_star, qw_star, &
         lwc, iwc, tl, ti, tzn

    ! Bus pointer declarations
    real, pointer, dimension(:), contiguous :: zbt_ag, zfv_ag, zqsurf_ag, zfvap_ag, &
         ps, zalfat, zalfaq, zbm, zbm0, zfdsi, zfdss, zfq, ztsrad, zustress, &
         zvstress, zue, zflw, zfsh, zdxdy, ztdmaskxdt
    real, pointer, dimension(:,:), contiguous :: tu, tv, tt, tq, uu, vv, w, &
         t, q, zsigt, zsigm, tm, zpblsigs, zqcplus, tqc, zze, zkm, zkt, &
         zzd, zzn, zfc, zfv, zpri, zturbqf, zturbtf, zturbuf, zturbvf, zfxp, &
         zturbuvf, zmrk2, ztcons, zqcons, zdtcons, zdqcons, zturbtlf, zgztherm
    real, pointer, dimension(:,:,:), contiguous :: zvcoef

    ! Initialize local arrays to detect undefined fields
    call init2nan(bmsg, btsg, aq, bq, rhosfc, fm_mult, fh_mult)
    call init2nan(kmsg, ktsg, dket, tl, ti, lwc, iwc, tzn)
    call init2nan(zero)

    ! Initialize array bounds for bus pointer macros
    ni = F_ni; nk = F_nk ; nkm1 = F_nkm1

    ! Pointer assignments
    MKPTR1D(ps, pmoins, F_pvars)
    MKPTR1D(zalfaq, alfaq, F_pvars)
    MKPTR1D(zalfat, alfat, F_pvars)
    MKPTR1D(zbm, bm, F_pvars)
    MKPTR1D(zbm0, bm0, F_pvars)
    MKPTR1D(zdxdy, dxdy, F_pvars)
    MKPTR1D(zfdsi, fdsi, F_pvars)
    MKPTR1D(zfdss, fdss, F_pvars)
    MKPTR1D(zflw, flw, F_pvars)
    MKPTR1D(zfq, fq, F_pvars)
    MKPTR1D(zfsh, fsh, F_pvars)
    MKPTR1D(ztsrad, tsrad, F_pvars)
    MKPTR1D(ztdmaskxdt, tdmaskxdt, F_pvars)
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
    MKPTR2Dm1(zfxp, fxp, F_pvars)
    MKPTR2Dm1(zgztherm, gztherm, F_pvars)
    MKPTR2Dm1(zkm, km, F_pvars)
    MKPTR2Dm1(zkt, kt, F_pvars)
    MKPTR2Dm1(zpri, pri, F_pvars)
    MKPTR2Dm1(zqcons, qcons, F_pvars)
    MKPTR2Dm1(ztcons, tcons, F_pvars)
    MKPTR2Dm1(zzd, zd, F_pvars)
    MKPTR2Dm1(zze, ze, F_pvars)
    MKPTR2Dm1(zpblsigs, pblsigs, F_pvars)
    
    MKPTR3D(zvcoef, vcoef, F_pvars)

    if (znplus > 0) then
       MKPTR2Dm1(zzn, znplus, F_pvars)
    else
       MKPTR2Dm1(zzn, zn, F_pvars)
    endif

    ! Initializations
    zero(:) = 0.
    
    ! Retrieve stochastic parameter information on request
    fm_mult(:) = ens_spp_get('fm_mult', zmrk2, default=1.)
    fh_mult(:) = ens_spp_get('fh_mult', zmrk2, default=1.)

    ! Convert diffusion coefficients to sigma coordinates
    call ktosig(kmsg, zkm, t, q, zsigm, zsigt, zvcoef, K_M, F_ni, F_nkm1)
    call ktosig(ktsg, zkt, t, q, zsigm, zsigt, zvcoef, K_T, F_ni, F_nkm1)

    ! Compute turbulent flux surface boundary conditions
    do i=1,F_ni
       aq(i)     = -(GRAV/RGASD) / (t(i,F_nkm1) * (1. + DELTA * zqsurf_ag(i)))
       bmsg(i)   = fm_mult(i) * zbm(i) * aq(i)
       btsg(i)   = zbt_ag(i) * aq(i)
       zalfat(i) = zalfat(i)  * aq(i)
       zalfaq(i) = zalfaq(i)  * aq(i)
       rhosfc(i) = -aq(i) * ps(i)/GRAV
    end do

    ! Compute condensate masses
    if (mp_lwc(lwc, F_pvars) /= PHY_OK .or. &
         mp_iwc(iwc, F_pvars) /= PHY_OK) then
       call physeterror('pbl_ri_diffuse::diffuseall', 'Cannot compute condensate contents')
       return
    endif
    
    ! Compute conserved thermodynamic variables for vertical diffusion
    call compute_conserved(ztcons, zqcons, t, q, lwc, iwc, zsigt, F_ni, F_nkm1)
    
    ! Compute tendency for u-component wind and diagnose flux on request
    call diffuse(tu, uu, kmsg, zero, bmsg, zsigm, zsigt, F_tau, &
         FIELD_ON_MLEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any([(any((/'TFUU', 'TFUV'/) == phyoutlist_S(i)), i=1,nphyoutlist)])) then
       call atmflux4(zturbuf, tu, zsigm, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Compute tendency for v-component wind and diagnose flux on request
    call diffuse(tv, vv, kmsg, zero, bmsg, zsigm, zsigt, F_tau, &
         FIELD_ON_MLEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any([(any((/'TFVV', 'TFUV'/) == phyoutlist_S(i)), i=1,nphyoutlist)])) then
       call atmflux4(zturbvf, tv, zsigm, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Diagnose total turbulent momentum flux
    if (any('TFUV' == phyoutlist_S)) &
         zturbuvf(:,:) = sqrt(zturbuf(:,:)**2+zturbvf(:,:)**2)
    
    ! Set surface flux boundary terms for thermal diffusion
    aq(:) = fh_mult(:) * zalfat(:)
    bq(:) = fh_mult(:) * btsg(:)

    ! Compute tendency for theta_l and diagnose flux on request 
    call diffuse(zdtcons, ztcons, ktsg, aq, bq, zsigm, zsigt, F_tau, &
         FIELD_ON_TLEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any(phyoutlist_S == 'TFTL')) then
       call atmflux4(zturbtlf, zdtcons, zsigt, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif    

    ! Set surface flux boundary terms for specific humidity diffusion
    aq(:) = fh_mult(:) * zalfaq(:)
    bq(:) = fh_mult(:) * btsg(:)

    ! Compute tendency for moisture and diagnose flux on request
    call diffuse(tq, q, ktsg, aq, bq, zsigm, zsigt, F_tau, &
         FIELD_ON_TLEV, F_ni, F_nkm1)
    if (phy_error_L) return
    if (any(phyoutlist_S == 'TFHU')) then
       call atmflux4(zturbqf, tq, zsigt, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif
    
    ! Compute diffusion tendency for condensate
    call diffuse_condensate(tl, ti, lwc, iwc, ktsg, zsigm, zsigt, ztdmaskxdt, &
         F_pvars, F_tau, F_ni, F_nkm1)

    ! Update conserved variables
    do k=1,F_nkm1
       do i=1,F_ni
          tqc(i,k) = tl(i,k) + ti(i,k)
          thl_star(i,k) = ztcons(i,k) + F_tau*zdtcons(i,k)
          qw_star(i,k) = zqcons(i,k) + F_tau*(tq(i,k) + tqc(i,k))
       enddo
    enddo
    
    ! Convert conserved variable tendencies to temperature tendency
    call baktot(tt, zdtcons, tl, ti, zsigt, F_ni, F_nkm1)
    
    ! Estimate updated subgrid-scale variance
    call sgsvar(zpblsigs, thl_star, qw_star, lwc, iwc, tm, zfxp, zpri, &
         zzn, zzd, zgztherm, zsigt, ps, zdxdy, zvcoef, F_ni, F_nkm1)
    
    ! Dissipative heating
    call dissheat(dket, uu, vv, tu, tv, kmsg, zsigm, zsigt, zvcoef, F_tau, F_ni, F_nkm1)
    tt(:,1:F_nkm1) = tt(:,1:F_nkm1) - (1./CPD) * dket(:,1:F_nkm1)

    ! Diagnose atmospheric fluxes for dry-air temperature
    if (any(phyoutlist_S == 'TFTT')) then
       call atmflux4(zturbtf, tt, zsigt, ps, F_ni, F_nkm1, F_type=FLUX_INTTYPE)
       if (phy_error_L) return
    endif

    ! Diagnose (implicit) surface fluxes based on updated atmospheric state
    call sfcflux(zustress, zvstress, zfq, zue, zfsh, zflw, &
         uu, vv, t, q, tu, tv, tt, tq, bmsg, btsg, zalfaq, &
         rhosfc, ps, F_tau, F_ni, F_nkm1)
    if (phy_error_L) return

    ! End of subprogram
    return
  end subroutine diffuseall
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine diffuse(F_tend, F_fld, F_kcoef, F_alpha, F_beta, F_sigm, F_sigt, &
       F_tau, F_type, F_ni, F_nkm1, F_ngflux)
    use, intrinsic :: iso_fortran_env, only: REAL64
    use debug_mod, only: init2nan
    use tridiag, only: td_solve
    implicit none

    ! Argument description
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fld           !field to diffuse
    real, dimension(F_ni,F_nkm1), intent(in) :: F_kcoef         !diffusion coefficients at layer interface
    real, dimension(F_ni), intent(in) :: F_alpha                !field-independent sfc flux term (inhomogeneous)
    real, dimension(F_ni), intent(in) :: F_beta                 !field-dependent sfc flux term (homogeneous)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigm          !sigma values of momentum levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma values of thermo levels (m-layer interfaces)
    real, intent(in) :: F_tau                                   !time step (s)
    integer, intent(in) :: F_type                               !type of field to diffuse (use API constants)
    real, dimension(F_ni,F_nkm1), intent(in), optional :: F_ngflux !nongradient fluxes at layer interface
    real, dimension(F_ni,F_nkm1), intent(out) :: F_tend         !field tendency due to vertical diffusion

    !@Object:  Solve a vertical diffusion equation

    ! Internal variable declarations
    integer :: i, k
    real :: hm, hp, kp, hk
    real, dimension(F_ni,F_nkm1) :: a, b, c, d, fldupd, ngflux

    ! Initialize local arrays to detect undefined fields
    call init2nan(a, b, c, d, fldupd, ngflux)

    ! Optional argument handling
    if (present(F_ngflux)) then
       ngflux(:,:) = F_ngflux(:,:)
    else
       ngflux(:,:) = 0.
    endif
    
    ! Construct the tridiagonal matrix (N=[A,B,C]) and forcing (D)
    STAGGERING: select case(F_type)

    case (FIELD_ON_TLEV)
       ! Field is on thermo levels (non-gradient fluxes on momentum levels)

       ! Top level
       do i=1,F_ni
          hp = F_sigt(i,1) - F_sigt(i,2)
          hk = F_sigm(i,1) - F_sigm(i,2)
          a(i,1) = 0.
          c(i,1) = -F_tau * F_kcoef(i,2) / (hp*hk)
          b(i,1) = 1. - (a(i,1) + c(i,1))
          d(i,1) = F_fld(i,1) + &
               F_tau * (ngflux(i,1) - ngflux(i,2)) / hk
       enddo

       ! Interior of the column
       do k=2,F_nkm1-1
          do i=1,F_ni
             hm = F_sigt(i,k-1) - F_sigt(i,k)
             hp = F_sigt(i,k) - F_sigt(i,k+1)
             hk = F_sigm(i,k) - F_sigm(i,k+1)
             a(i,k) = -F_tau * F_kcoef(i,k) / (hm*hk)
             c(i,k) = -F_tau * F_kcoef(i,k+1) / (hp*hk)
             b(i,k) = 1. - (a(i,k) + c(i,k))
             d(i,k) = F_fld(i,k) + &
                  F_tau * (ngflux(i,k) - ngflux(i,k+1)) / hk
          enddo
       enddo

       ! Bottom level
       do i=1,F_ni
          hm = F_sigt(i,F_nkm1-1) - F_sigt(i,F_nkm1)
          hk = F_sigm(i,F_nkm1) - 1.
          a(i,F_nkm1) = -F_tau * F_kcoef(i,F_nkm1) / (hm*hk)
          c(i,F_nkm1) = 0.
          b(i,F_nkm1) = 1. - (a(i,F_nkm1) - F_tau * F_beta(i) / hk)
          d(i,F_nkm1) = F_fld(i,F_nkm1) + &
               F_tau * (ngflux(i,F_nkm1) - F_alpha(i)) / hk
       enddo
       
    case (FIELD_ON_MLEV, FIELD_IS_TKE)
       ! Field is on momentum levels (non-gradient fluxes on thermo levels)

       ! Top level
       do i=1,F_ni
          hp = F_sigm(i,1) - F_sigm(i,2)
          hk = F_sigm(i,1) - F_sigt(i,1)
          a(i,1) = 0.          
          c(i,1) = -F_tau * F_kcoef(i,1) / (hp*hk)
          b(i,1) = 1. - (a(i,1) + c(i,1))
          d(i,1) = F_fld(i,1) + &
               F_tau * ngflux(i,1) / hk
       enddo

       ! Interior of the column 
       do k=2,F_nkm1-1
          do i=1,F_ni
             hm = F_sigm(i,k-1) - F_sigm(i,k)
             hp = F_sigm(i,k) - F_sigm(i,k+1)
             hk = F_sigt(i,k-1) - F_sigt(i,k)
             a(i,k) = -F_tau * F_kcoef(i,k-1) / (hm*hk)
             c(i,k) = -F_tau * F_kcoef(i,k) / (hp*hk)
             b(i,k) = 1. - (a(i,k) + c(i,k))
             d(i,k) = F_fld(i,k) + &
                  F_tau * (ngflux(i,k-1) - ngflux(i,k)) / hk
          enddo
       enddo

       ! Bottom level
       BOTTOM_ELEV: if (F_type == FIELD_IS_TKE) then
          ! Bottom level is an energy level (surface)
          k = F_nkm1
          do i=1,F_ni
             hm = F_sigm(i,k-1) - F_sigm(i,k)
             hp = F_sigm(i,k) - 1.
             hk = F_sigt(i,k-1) - F_sigt(i,k)
             kp = F_kcoef(i,k) / (2.*sqrt(2.)) * sqrt(1. + F_alpha(i)/(F_fld(i,k)))
             a(i,k) = -F_tau * F_kcoef(i,F_nkm1-1) / (hm*hk)
             c(i,k) = 0.
             b(i,k) = 1. - (a(i,k) - F_tau * kp / (hp*hk))
             d(i,k) = F_fld(i,k) + F_tau * kp * F_alpha(i) / (hp*hk)
          enddo   
       else
          ! Bottom level is above-ground
          do i=1,F_ni
             hk = F_sigt(i,F_nkm1-1) - 1.
             hm = F_sigm(i,F_nkm1-1) - F_sigm(i,F_nkm1)
             a(i,F_nkm1) = -F_tau * F_kcoef(i,F_nkm1-1) / (hm*hk)
             c(i,F_nkm1) = 0.
             b(i,F_nkm1) = 1. - (a(i,F_nkm1) - F_tau * F_beta(i) / hk)
             d(i,F_nkm1) = F_fld(i,F_nkm1) - &
                  F_tau * (ngflux(i,F_nkm1) + F_alpha(i)) / hk
          enddo
       endif BOTTOM_ELEV

    case DEFAULT
       call physeterror('pbl_ri_diffuse::diffuse', 'Invalid staggering (F_type) selected')
       return

    end select STAGGERING
    
    ! Solve the tridiagonal system for the updated field value
    call td_solve(fldupd, a, b, c, d, F_ni, F_nkm1)

    ! Retrieve tendency
    F_tend = (fldupd - F_fld) / F_tau

    return
  end subroutine diffuse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine diffuse_condensate(F_dliq, F_dice, F_lwc, F_iwc, F_ktsg, F_sigm, F_sigt, &
       F_tdmaskxdt, F_pvars, F_tau, F_ni, F_nkm1)
    use phymem, only: phyvar, phymem_find
    use microphy_utils, only: mp_todiffuse, mp_lwc, mp_iwc
    use phy_status, only: PHY_OK
    use tendency, only: apply_tendencies

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_ktsg          !diffusion coefficient for scalars
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigm          !sigma values of momentum levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma values of thermo levels
    real, dimension(F_ni), intent(in) :: F_tdmaskxdt            !tendency mask times F_tau
    type(phyvar), pointer, contiguous :: F_pvars(:)             !physics buses
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_lwc        !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_iwc        !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dliq         !liquid condensate tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dice         !ice condensate tendency (kg/kg/s)

    ! Local parameters
    integer, parameter :: NMAX_DIFF=256
    
    ! Local variables
    integer :: i, k
    integer, dimension(:), pointer :: dindxl=>null()
    real, dimension(F_ni) :: zero
    real, dimension(F_ni,F_nkm1) :: lwc, iwc, tnd
    real, dimension(:,:), pointer :: busptr3d

    ! Initializations
    zero(:) = 0.
    
    ! Retrieve list of fields to diffuse
    if (mp_todiffuse(dindxl) /= PHY_OK) then
       call physeterror('pbl_ri_diffuse::diffuse_condensate', &
            'Cannot retrieve list of microphysics fields to diffuse')
       return
    endif

    ! Diffuse requested fields
    do i=1,size(dindxl)
       busptr3d(1:F_ni,1:F_nkm1) => F_pvars(dindxl(i))%data(:)
       call diffuse(tnd, busptr3d, F_ktsg, zero, zero, F_sigm, F_sigt, F_tau, &
            FIELD_ON_TLEV, F_ni, F_nkm1)
       call apply_tendencies(busptr3d, tnd, F_tdmaskxdt, F_ni, F_nkm1)
    enddo

    ! Recompute condensate masses
    if (mp_lwc(lwc, F_pvars) /= PHY_OK .or. &
         mp_iwc(iwc, F_pvars) /= PHY_OK) then
       call physeterror('pbl_ri_diffuse::diffuse_condensate', &
            'Cannot compute condensate contents')
       return
    endif

    ! Compute condensate tendency and contents
    do k=1,F_nkm1
       do i=1,F_ni
          F_dliq(i,k) = (lwc(i,k) - F_lwc(i,k)) / F_tau
          F_dice(i,k) = (iwc(i,k) - F_iwc(i,k)) / F_tau
          F_lwc(i,k) = lwc(i,k)
          F_iwc(i,k) = iwc(i,k)
       enddo
    enddo

  end subroutine diffuse_condensate
  
end module pbl_ri_diffuse
