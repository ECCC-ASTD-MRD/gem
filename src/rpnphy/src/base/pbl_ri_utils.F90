module pbl_ri_utils
  use vintphy, only: vint_mom2thermo, vint_thermo2mom
  implicit none
  private

  ! API functions
  public :: baktot                                              !Convert conserved variable tendencies to temperature
  public :: buoyflux                                            !Buoyancy flux calculation
  public :: covareq                                             !Equilibrium solution for SGS variance
  public :: covarstep                                           !Algebraic solution of the SGS variance equation
  public :: kcalc                                               !Compute diffusion coefficient
  public :: ktosig                                              !Convert coefficients from height to sigma coord
  public :: sgsvar                                              !Diagnose subgrid-scale variance (sigma_s)
  public :: tkestep                                             !Algebraic solution of the TKE equation
  public :: tkebudget                                           !Diagnose terms of the TKE budget

  ! API parameters
  integer, parameter, public :: K_M = 0
  integer, parameter, public :: K_T = 1
  real, parameter, public :: PBL_RI_CU = 3.75                   !Friction velocity coefficient
  real, parameter, public :: PBL_RI_CW = 0.2                    !Convective scale coefficient
  real, parameter, public :: PBL_RI_CK = 1./sqrt(PBL_RI_CU)     !Coefficient for diffusivity
  real, parameter, public :: PBL_RI_CE = 1./PBL_RI_CU**1.5      !Dissipation coefficient ce*ck*cu^2=l_e/l=1 for Ri=0
  real, parameter, public :: PBL_RI_CAB = 2.5                   !Coefficient for second-order moments (BFP92 Eq. 19)
  real, parameter, public :: PBL_SDMIN=1e-10                    !Minimum value for unscaled SGS stddev

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine baktot(F_dt, F_dthl, F_dliq, F_dice, F_sigt, F_ni, F_nkm1)
    use tdpack_const, only: CPD, CAPPA, CHLC, CHLF
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dthl          !tendency of theta_l (K/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dliq          !liquid condensate tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dice          !ice condensate tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dt           !tendency of dry air temperature (K/s)

    !@Object
    !          Compute temperature tendency from conserved variables

    ! Local variables
    integer :: i, k
 
    ! Convert back to state variables and tendencies
    do k=1,F_nkm1
       do i=1,F_ni
          F_dt(i,k) = F_dthl(i,k) * F_sigt(i,k)**CAPPA + &
               (CHLC * F_dliq(i,k) + (CHLC+CHLF) * F_dice(i,k)) / CPD
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine baktot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine buoyflux(F_dthv, F_dudz2, F_ri, F_uu, F_vv, F_tt, F_hu, &
       F_lwc, F_iwc, F_fbl, F_tsurf, F_qsurf, F_z0m, F_z0t, &
       F_lat, F_fcor, F_sigt, F_gzm, F_gzt, F_ps, F_vcoef, F_ni, F_nkm1)
    use tdpack_const, only: DELTA, GRAV, CPD, CAPPA
    use phy_options
    use cons_thlqw, only: thlqw_compute, thlqw_thermco
    use sfclayer, only: sl_prelim, sl_sfclayer, SL_OK

    implicit none

    ! Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_uu            !u-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_vv            !v-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_hu            !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc           !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc           !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fbl           !PBL cloud fraction
    real, dimension(F_ni), intent(in) :: F_tsurf                !surface air temperature (K)
    real, dimension(F_ni), intent(in) :: F_qsurf                !surface specific humidity (kg/kg)
    real, dimension(F_ni), intent(in) :: F_z0m                  !momentum roughness length (m)
    real, dimension(F_ni), intent(in) :: F_z0t                  !thermal roughness length (m)
    real, dimension(F_ni), intent(in) :: F_lat                  !latitude (rad)
    real, dimension(F_ni), intent(in) :: F_fcor                 !Coriolis factor (/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzm           !heights of momentum levels (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzt           !heights of thermo levels (m)
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dudz2        !square of vertical wind shear (s^-2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_ri           !gradient Richardson number
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dthv         !square of buoyancy frequency (s^-2)

    !Object
    !          Compute buoyancy and shear (squared) gradients for TKE sources

    ! External common blocks
    include "surface.cdk"
    
    ! Local variable declarations
    integer :: i, k
    real :: idz, dqwdz, dthldz, exner, alpha, bet, mu, dudz, dvdz, twc
    real, dimension(F_ni) :: wspd, wdir, theta
    real, dimension(F_ni,F_nkm1) :: thl, qw, acoef, bcoef, coefthl, coefqw, &
         thv, ut, vt, leff, thvm

    ! Compute thermodynamic coefficients on thermo levels following Bechtold and Siebsma (JAS 1998)
    call thlqw_compute(thl, qw, F_tt, F_hu, F_lwc, F_iwc, F_sigt, F_ni, F_nkm1)
    call thlqw_thermco(thl, qw, F_tt, F_lwc, F_iwc, F_sigt, F_ps, &
         F_ni, F_nkm1, F_acoef=acoef, F_bcoef=bcoef, F_leff=leff)

    ! Compute terms of buoyancy flux equation (Eq. 4 of Bechtold and Siebsma (JAS 1998))
    do k=1,F_nkm1
       do i=1,F_ni
          exner = F_sigt(i,k)**CAPPA
          theta(i) = F_tt(i,k) / exner
          alpha = DELTA * theta(i)
          twc = F_lwc(i,k) + F_iwc(i,k)
          mu = DELTA*qw(i,k) - (1. + DELTA)*twc
          bet = leff(i,k)/CPD / exner * (1. + mu) - (1. + DELTA)*theta(i)
          thv(i,k) = thl(i,k) + alpha*qw(i,k) + bet*twc
          coefthl(i,k) = 1.0 + mu - bet*bcoef(i,k)*F_fbl(i,k)          
          coefqw(i,k) = alpha + bet*acoef(i,k)*F_fbl(i,k)
       enddo
    enddo

    ! Interpolate inputs to required level sets
    call vint_mom2thermo(ut,   vt, &
         &               F_uu, F_vv, F_vcoef, F_ni, F_nkm1)
    call vint_thermo2mom(thvm, coefthl, coefqw, &
         &               thv,  coefthl, coefqw, F_vcoef, F_ni, F_nkm1)
    if (sl_prelim(F_tt(:,F_nkm1), F_hu(:,F_nkm1), F_uu(:,F_nkm1), F_vv(:,F_nkm1), &
         F_ps, F_gzm(:,F_nkm1), spd_air=wspd, dir_air=wdir, min_wind_speed=VAMIN) /= SL_OK) &
         call physeterror('pbl_ri_utils::buoyflux', 'Problem preparing surface layer')
    do i=1,F_ni
       theta(i) = F_tt(i,F_nkm1) / F_sigt(i,F_nkm1)**CAPPA
    enddo
    if (sl_sfclayer(theta, F_hu(:,F_nkm1), wspd, wdir, F_gzm(:,F_nkm1), &
         F_gzt(:,F_nkm1), F_tsurf, F_qsurf, F_z0m, F_z0t, F_lat, F_fcor, &
         hghtm_diag_row=F_gzt(:,F_nkm1), u_diag=ut(:,F_nkm1), v_diag=vt(:,F_nkm1)) /= SL_OK) &
         call physeterror('pbl_ri_utils::buoyflux', 'Problem computing surface layer')
    
    ! Compute components of the Richardson number on momentum levels
    do k=2,F_nkm1
       do i=1,F_ni
          ! Vertical derivatives of conserved variables
          idz = 1. / (F_gzt(i,k) - F_gzt(i,k-1))
          dthldz = (thl(i,k) - thl(i,k-1)) * idz
          dqwdz = (qw(i,k) - qw(i,k-1)) * idz
          dudz = (ut(i,k) - ut(i,k-1)) * idz
          dvdz = (vt(i,k) - vt(i,k-1)) * idz
          ! Buoyancy frequency squared
          F_dthv(i,k) = (coefthl(i,k)*dthldz + coefqw(i,k)*dqwdz) * (GRAV/thvm(i,k))
          ! Shear-squared
          F_dudz2(i,k) = dudz**2 + dvdz**2
          ! Richardson number
          F_ri(i,k) = F_dthv(i,k) / max(F_dudz2(i,k), 1e-6)
       enddo
    enddo
    do i=1,F_ni
       F_dthv(i,1) = 0.
       F_dudz2(i,1) = 0.
       F_ri(i,1) = 0.
    enddo

    return
  end subroutine buoyflux
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine covareq(F_coveq, F_flda, F_fldb, F_en, F_kt, F_zd, &
       F_gzt, F_vcoef, F_ni, F_nkm1, F_taud)
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_flda          !gridscale field "a" value
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fldb          !gridscale field "b" value
    real, dimension(F_ni,F_nkm1), intent(in) :: F_en            !turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_kt            !diffusion coefficient (m2/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zd            !dissipation length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzt           !heights of thermo levels (m)
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(out) :: F_coveq        !equilibrium SGS (a'b') covariance
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_taud !dissipation time scale (s)

    !@Objective Compute equilibrium covariance of fields "a" and "b"

    ! Local variables
    integer :: i, k
    real, dimension(F_ni,F_nkm1) :: taud, dabdz

    ! Compute equilibrium state based on vertical gradients
    do k=2,F_nkm1
       do i=1,F_ni
          dabdz(i,k) = (F_flda(i,k-1) - F_flda(i,k)) * (F_fldb(i,k-1) - F_fldb(i,k)) / &
               (F_gzt(i,k-1) - F_gzt(i,k))**2
       enddo
    enddo
    dabdz(:,1) = dabdz(:,2)
    do k=1,F_nkm1
       do i=1,F_ni
          taud(i,k) = F_zd(i,k) / (PBL_RI_CE * sqrt(F_en(i,k)))
          F_coveq(i,k) = 2. * F_kt(i,k) * taud(i,k) * dabdz(i,k)          
       enddo
    enddo
    call vint_mom2thermo(F_coveq, taud, &
         &               F_coveq, taud, F_vcoef, F_ni, F_nkm1)
    
    ! Return dissipation timescale on request
    if (present(F_taud)) F_taud(:,:) = taud(:,:)

    ! End of subprogram
    return
  end subroutine covareq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine covarstep(F_covstar, F_cov, F_flda, F_fldb, &
       F_en, F_kt, F_zd, F_gzt, F_vcoef, F_tau, F_ni, F_nkm1)
    implicit none
    
    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_cov           !current field SGS (a'b') covariance
    real, dimension(F_ni,F_nkm1), intent(in) :: F_flda          !gridscale field "a" value
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fldb          !gridscale field "b" value
    real, dimension(F_ni,F_nkm1), intent(in) :: F_en            !turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_kt            !diffusion coefficient (m2/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zd            !dissipation length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzt           !heights of thermo levels (m)
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_covstar      !updated (a'b') covariance

    !@Object Update covariance of "a" and "b"

    ! Local variables
    integer :: i, k
    real :: expt
    real, dimension(F_ni,F_nkm1) :: taud, xe

    ! Compute equilibrium state
    call covareq(xe, F_flda, F_fldb, F_en, F_kt, F_zd, F_gzt, &
         F_vcoef, F_ni, F_nkm1, F_taud=taud)

    ! Compute updated covariance
    do k=1,F_nkm1
       do i=1,F_ni
          expt = exp(-F_tau / taud(i,k))
          F_covstar(i,k) = F_cov(i,k) * expt + xe(i,k)*(1. - expt)
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine covarstep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kcalc(F_kcoef, F_zn, F_en, F_pri, F_vcoef, F_type, F_ni, F_nkm1)
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zn            !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_en            !turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_pri           !inverse Prandtl number
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    integer, intent(in) :: F_type                               !type of coefficient (K_M or K_T)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_kcoef        !diffusion coefficient (m2/s)

    !@Object Compute diffusion coefficients

    ! Local variables
    integer :: i, k
    real, dimension(F_ni,F_nkm1) :: kt, znt, ent
    
    ! Compute coefficients on momentum levels
    do k=1,F_nkm1
       do i=1,F_ni
          F_kcoef(i,k) = PBL_RI_CK * F_zn(i,k) * sqrt(F_en(i,k))          
       enddo
    enddo

    ! Interpolate to thermo levels if required
    if (F_type == K_M) then
!!$       do k=1,F_nkm1-1
!!$          do i=1,F_ni
!!$             kt(i,k) =  0.5*(F_kcoef(i,k) + F_kcoef(i,k+1))  !arithmetic
!!$!             kt(i,k) = sqrt(F_kcoef(i,k) * F_kcoef(i,k+1))  !geometric
!!$!             kt(i,k) = 2. / (1./F_kcoef(i,k) + 1./F_kcoef(i,k+1))  !harmonic
!!$          enddo
!!$       enddo
       call vint_mom2thermo(znt,  ent, &
            &               F_zn, F_en, F_vcoef, F_ni, F_nkm1)
       do k=1,F_nkm1-1
          do i=1,F_ni
             kt(i,k) = PBL_RI_CK * znt(i,k) * sqrt(ent(i,k))
          enddo
       enddo
       ! Lowest thermo-level K is interpolated from 0 at the surface
       kt(:,F_nkm1) = 0.5 * F_kcoef(:,F_nkm1)
       F_kcoef(:,:) = kt(:,:)
    endif

    ! Apply Prandtl-number adjustment if required
    if (F_type == K_T) F_kcoef(:,:) = F_kcoef(:,:) * F_pri(:,:)

  end subroutine kcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ktosig(F_ksig, F_kcoef, F_tt, F_hu, F_sigm, F_sigt, &
       F_vcoef, F_type, F_ni, F_nkm1)
    use tdpack, only: GRAV, RGASD
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_kcoef         !diffusion coefficient in height coord (m2/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_hu            !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), target, intent(in) :: F_sigm  !sigma of momentum levels
    real, dimension(F_ni,F_nkm1), target, intent(in) :: F_sigt  !sigma of thermodynamic levels
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    integer, intent(in) :: F_type                               !type of coefficient (K_M or K_T)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_ksig         !diffusion coefficient in sigma coord

    !@Object Convert height-coordinate diffusivity to sigma coordinates
    
    ! Local variables
    integer :: i, k
    real :: goverr
    real, dimension(F_ni,F_nkm1) :: tv
    real, dimension(:,:), pointer :: sig

    ! Initializations
    goverr = GRAV/RGASD
    if (F_type == K_T) then
       sig => F_sigm
    else
       sig => F_sigt
    endif
    
    ! Compute virtual temperature on appropriate levels    
    call mfotvt(tv, F_tt, F_hu, F_ni, F_nkm1, F_ni)
    if (F_type == K_T) call vint_thermo2mom(tv, tv, F_vcoef, F_ni, F_nkm1)

    ! Convert diffusion coefficients
    do k=1,F_nkm1
       do i=1,F_ni
          F_ksig(i,k) = F_kcoef(i,k) * (goverr * sig(i,k)/tv(i,k))**2          
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine ktosig
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sgsvar(F_sigmas, F_thl, F_qw, F_lwc, F_iwc, F_tt, F_pri, &
       F_zn, F_ze, F_gzt, F_sigt, F_ps, F_dxdy, F_vcoef, F_ni, F_nkm1)
    use tdpack_const
    use pbl_utils, only: dvrtdf
    use cons_thlqw, only: thlqw_thermco
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_thl           !liquid water potential temperature (K; theta_l)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw            !total water mixing ratio (kg/kg; q_tot)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc           !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc           !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_pri           !inverse Prandtl number
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zn            !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_ze            !dissipation length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzt           !heights of thermo levels (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(F_ni), intent(in) :: F_dxdy                 !grid cell area (m2)
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(out) :: F_sigmas       !subgrid moisture variance (kg/kg)

    !@Object Calculate the boundary layer subgrid-scale properties

    ! Local parameters
    real, parameter :: REF_DX=500.                              !grid mesh for variance scaling

    ! Local variables
    integer :: i,k
    real :: c1coef, idz
    real, dimension(F_ni,F_nkm1) :: dqwdz, dthldz, acoef, bcoef, ccoef
    
    ! Pre-compute vertical derivatives of conserved variables (valid on m-levels)
    do k=2,F_nkm1
       do i=1,F_ni
          idz = 1. / (F_gzt(i,k) - F_gzt(i,k-1))
          dthldz(i,k) = (F_thl(i,k) - F_thl(i,k-1)) * idz
          dqwdz(i,k) = (F_qw(i,k) - F_qw(i,k-1)) * idz
       enddo
    enddo

    ! Compute thermodynamic coefficients from conserved variables
    call thlqw_thermco(F_thl, F_qw, F_tt, F_lwc, F_iwc, F_sigt, F_ps, &
         F_ni, F_nkm1, F_acoef=acoef, F_bcoef=bcoef, F_ccoef=ccoef)
    
    ! Compute subgrid standard deviation of supersaturation (s) following Eq. 10 of BCMT95
    call vint_thermo2mom(acoef, bcoef, &
         &               acoef, bcoef, F_vcoef, F_ni, F_nkm1)
    do k=2,F_nkm1
       do i=1,F_ni
          c1coef = sqrt(2. * PBL_RI_CK * F_pri(i,k) / PBL_RI_CAB * F_zn(i,k) * F_ze(i,k))
          F_sigmas(i,k) = sqrt(sqrt(F_dxdy(i))/REF_DX) * max(c1coef * abs(acoef(i,k)*dqwdz(i,k) &
               - bcoef(i,k)*dthldz(i,k)), PBL_SDMIN)
       end do
    end do
    do i=1,F_ni
       F_sigmas(i,1) = sqrt(sqrt(F_dxdy(i))/REF_DX) * PBL_SDMIN
    enddo
    call vint_mom2thermo(F_sigmas, F_sigmas, F_vcoef, F_ni, F_nkm1)
    
    return
  end subroutine sgsvar
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tkestep(F_estar, F_en, F_b, F_c, F_tau, F_ni, F_nkm1)
    use phy_options, only: etrmin2
    implicit none

    ! Argument declarations
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_en            !current turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_b             !B (generation) term of TKE equation (m/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_c             !C (dissipation) term of TKE eqaution (1/m)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_estar        !updated turbulent kinetic energy (m2/s2)

    !Object
    !          Update TKE using the analytical part of the TKE equation
    
    ! Local variable declarations
    integer :: i, k
    complex :: ye, itbc, xi0, gamma, yz
    real :: sgn, tmax, bsafe

    ! Update TKE using an analytical solution of the non-diffusive TKE equation
    do k=1,F_nkm1
       do i=1,F_ni
          bsafe = sign(max(abs(F_b(i,k)), etrmin2/F_tau), F_b(i,k))
          ye = sqrt(cmplx(bsafe / F_c(i,k)))
          itbc = 0.5 * sqrt(cmplx(bsafe * F_c(i,k)))  
          xi0 = sqrt(F_en(i,k)) / ye
          sgn = sign(1., 1. - xi0%RE)
          if (abs(xi0%RE) == 1.) then
             F_estar(i,k) = F_en(i,k)
          else             
             gamma = atanh(xi0**sgn)
             tmax = -min(gamma%IM, -tiny(gamma%IM)*F_tau) / max(itbc%IM, tiny(itbc%IM))
             yz = ye * tanh(gamma + itbc * min(F_tau, tmax))**sgn
             F_estar(i,k) = max(yz%RE**2, etrmin2)
          endif
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine tkestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tkebudget(F_en, F_estar, F_b, F_c, F_zn, F_buoy_flux, F_dvdz2, &
       F_buoy, F_shr, F_diss, F_tau, F_ni, F_nkm1)
    implicit none

    ! Argument declarations
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_en            !current turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_estar         !updated turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_b             !B (generation) coefficient of TKE equation (m/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_c             !C (dissipation) termcoefficient of TKE equation (1/m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zn            !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_buoy_flux     !square of buoyancy frequency (1/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dvdz2         !square of vertical shear (1/s2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_buoy         !buoyancy generation/suppression term (m2/s3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_shr          !shear generation term (m2/s3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_diss         !viscous dissipation term (m2/s3)

    !Object
    !          Compute budget terms for TKE equation
    
    ! Local variable declarations
    integer :: i, k
    real :: clambda, bterm, b_over_c
          
    ! Compute TKE budget diagnostics
    ! AZ: if we use hysteresis, and depending on how it is implemented, this budget calculation
    !     may need to be revised - tbc
    do k=1,F_nkm1
       do i=1,F_ni
          b_over_c = F_b(i,k) / F_c(i,k)
          if( F_en(i,k) /= 0. .and. F_en(i,k) /= b_over_c .and. F_estar(i,k) /= b_over_c ) then
             bterm = log( abs( (F_en(i,k) - b_over_c) / (F_estar(i,k) - b_over_c) ) ) &
                    / (sqrt(F_en(i,k)) * F_c(i,k) * F_tau)
          else
             bterm = 1.
          endif
          clambda = PBL_RI_CK * F_zn(i,k) * sqrt(F_en(i,k)) * bterm
          F_shr(i,k) = clambda * F_dvdz2(i,k)
          F_buoy(i,k) = -clambda * F_buoy_flux(i,k)
          F_diss(i,k) = (F_estar(i,k) - F_en(i,k)) / F_tau - (F_shr(i,k) + F_buoy(i,k))
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine tkebudget  
  
end module pbl_ri_utils
