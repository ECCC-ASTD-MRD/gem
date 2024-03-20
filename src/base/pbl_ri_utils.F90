module pbl_ri_utils
  implicit none
  private

  ! API functions
  public :: baktotq                                             !Convert conserved variable tendencies to state tendencies
  public :: buoyflux                                            !Buoyancy flux calculation
  public :: sgsvar                                              !Diagnose subgrid-scale variance (sigma_s)
  public :: thermco                                             !Compute thermodynamic coefficients
  public :: tkestep                                             !Algebraic solution of the TKE equation
  public :: tkebudget                                           !Diagnose terms of the TKE budget
  public :: compute_conserved                                   !Compute conserved variables
  
  ! API parameters
  integer, parameter, public :: IMPLICIT_CLOUD = 0
  integer, parameter, public :: EXPLICIT_CLOUD = 1
  integer, parameter, public :: COMBINED_CLOUD = 2
  real, parameter, public :: PBL_RI_CU = 3.75                   !Friction velocity coefficient
  real, parameter, public :: PBL_RI_CW = 0.2                    !Convective scale coefficient
  real, parameter, public :: PBL_RI_CK = 1./sqrt(PBL_RI_CU)     !Coefficient for diffusivity
  real, parameter, public :: PBL_RI_CE = 0.714                  !Dissipation coefficient (BFP92 Eq. 8)
  real, parameter, public :: PBL_RI_CAB = 2.5                   !Coefficient for second-order moments (BFP92 Eq. 19)
  real, parameter, public :: PBL_RI_AE = 2.                     !Multiplicative factor on K for TKE diffusion

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine baktotq(F_dt, F_dqv, F_dqc, F_thl, F_qw, F_dthl, F_dqw, F_qc, F_sigt, &
       F_ps, F_tif, F_tve, F_fn, F_pri, F_zn, F_ze, &
       F_vcoef, F_pblsigs, F_pblq1, F_tau, F_ni, F_nkm1)
    use tdpack_const, only: CPD, CAPPA, CHLC, CHLF
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step length (s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_thl           !liquid water potential temperature (K; theta_l)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw            !total water mixing ratio (kg/kg; q_tot)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qc            !PBL cloud water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni),    intent(in) :: F_ps                !surface pressure (Pa)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tif           !temperature used for ice fraction (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dthl          !tendency of theta_l (K/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dqw           !tendency of q_tot (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tve           !virtual temperature on e-levs (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_pri           !inverse Prandtl number
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zn            !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_ze            !dissipation length (m)
    real, dimension(*), intent(in) :: F_vcoef                   !coefficients for vertical interpolation   
    real, dimension(F_ni,F_nkm1), intent(out) :: F_fn           !cloud fraction
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dt           !tendency of dry air temperature (K/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dqv          !tendency of specific humidity (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dqc          !tendency of cloud water content (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_pblsigs      !Subgrid moisture variance
    real, dimension(F_ni,F_nkm1), intent(out) :: F_pblq1        !Normalized saturation deficit

    !@Object
    !          Compute state variable tendencies from conserved variables

    ! Local variables
    integer :: i, k
    real :: tauinv
    real, dimension(F_ni,F_nkm1) :: exner, thl_star, qw_star, acoef, bcoef, ccoef, qcp, &
         qv, leff

    ! Update conserved variables after diffusion
    do k=1,F_nkm1
       do i=1,F_ni
          exner(i,k) = F_sigt(i,k)**CAPPA
          thl_star(i,k) = F_thl(i,k) + F_tau*F_dthl(i,k)
          qw_star(i,k) = F_qw(i,k) + F_tau*F_dqw(i,k)
       enddo
    enddo

    ! Compute thermodynamic coefficients from conserved variables
    call thermco(thl_star, qw_star, F_tif, F_fn, F_sigt, F_ps, IMPLICIT_CLOUD, F_ni, F_nkm1, &
         F_acoef=acoef, F_bcoef=bcoef, F_ccoef=ccoef, F_leff=leff)

    ! Compute subgrid-scale variances
    call sgsvar(thl_star, F_tve, qw_star, qcp, F_fn, F_pri, F_zn, F_ze, F_sigt, &
         acoef, bcoef, ccoef, F_vcoef, F_pblsigs, F_pblq1, F_ni, F_nkm1)

    ! Convert back to state variables and tendencies
    tauinv = 1./F_tau
    do k=1,F_nkm1
       do i=1,F_ni
          qv(i,k) = F_qw(i,k) - max(0.,F_qc(i,k))
          F_dqc(i,k) = (max(0.,qcp(i,k)) - max(F_qc(i,k), 0.)) * tauinv
          F_dqc(i,k) = max(F_dqc(i,k), -max(F_qc(i,k), 0.) * tauinv) !prevent negative values of qc
          F_dqc(i,k) = min(F_dqc(i,k), F_dqw(i,k) + max(qv(i,k), 0.) * tauinv) !prevent over-depletion of water vapour by PBL cloud
          F_dt(i,k) = exner(i,k)*F_dthl(i,k) + (leff(i,k)/CPD)*F_dqc(i,k)
          F_dqv(i,k) = max((F_dqw(i,k) - F_dqc(i,k)), -max(qv(i,k), 0.) * tauinv)
       enddo
    enddo
    
    return
  end subroutine baktotq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine buoyflux(F_u, F_v, F_t, F_tve, F_q, F_qc, F_fbl, F_dudz, F_dvdz, &
       F_sigt, F_ps, F_dudz2, F_ri, F_dthv, F_vcoef, F_ni, F_nkm1)
    use tdpack_const, only: DELTA, GRAV, RGASD, CPD, CAPPA
    use phy_options
    use vintphy, only: vint_mom2thermo2

    implicit none

    ! Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_u             !u-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_v             !v-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_t             !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tve           !virt. temperature on energy levels (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_q             !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qc            !PBL cloud water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(*), intent(in) :: F_vcoef                   !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fbl           !PBL cloud fraction
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dudz         !x-compontent vertical wind shear (s^-1)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dvdz         !y-compontent vertical wind shear (s^-1)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dudz2        !square of vertical wind shear (s^-2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_ri           !gradient Richardson number
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dthv         !buoyancy flux (m2/s2)

    !Object
    !          Compute buoyancy and shear (squared) gradients for TKE sources

    ! Local variable declarations
    integer :: i, k
    real :: idz, dqwdz, dthldz, exner, theta, alpha, beta
    real, dimension(F_ni,F_nkm1) :: thl, qw, acoef, bcoef, coefthl, coefqw, &
         thv, ut, vt, leff

    ! Compute thermodynamic coefficients on thermo levels following Bechtold and Siebsma (JAS 1998)
    call compute_conserved(thl, qw, F_t, F_q, F_qc, F_sigt, F_ni, F_nkm1)
    call thermco(thl, qw, F_t, F_fbl, F_sigt, F_ps, IMPLICIT_CLOUD, F_ni, F_nkm1, &
         F_acoef=acoef, F_bcoef=bcoef, F_leff=leff)

    ! Compute terms of buoyancy flux equation (Eq. 4 of Bechtold and Siebsma (JAS 1998))
    do k=1,F_nkm1
       do i=1,F_ni
          exner = F_sigt(i,k)**CAPPA
          theta = F_t(i,k) / exner
          alpha = DELTA * theta
          beta = leff(i,k)/CPD / exner - (1. + DELTA)*theta
          thv(i,k) = thl(i,k) + alpha*qw(i,k) + beta*F_qc(i,k)
          coefthl(i,k) = 1.0 + DELTA*qw(i,k) - beta*bcoef(i,k)*F_fbl(i,k)
          coefqw(i,k) = alpha + beta*acoef(i,k)*F_fbl(i,k)
       enddo
    enddo

    ! Compute components of the Richardson number on energy levels
    call vint_mom2thermo2(ut,  vt, &
         &                F_u, F_v, F_vcoef, F_ni, F_nkm1)
    do k=2,F_nkm1
       do i=1,F_ni
          ! Vertical derivatives of conserved variables
          idz = -GRAV / (RGASD * F_tve(i,k) * log(F_sigt(i,k) / F_sigt(i,k-1)))
          dthldz = (thl(i,k) - thl(i,k-1)) * idz
          dqwdz = (qw(i,k) - qw(i,k-1)) * idz
          F_dudz(i,k) = (ut(i,k) - ut(i,k-1)) * idz
          F_dvdz(i,k) = (vt(i,k) - vt(i,k-1)) * idz
          ! Buoyancy flux
          F_dthv(i,k) = (coefthl(i,k)*dthldz + coefqw(i,k)*dqwdz) * (GRAV/thv(i,k))
          ! Shear-squared
          F_dudz2(i,k) = F_dudz(i,k)**2 + F_dvdz(i,k)**2
          ! Richardson number
          F_ri(i,k) = F_dthv(i,k) / max(F_dudz2(i,k), 1e-6)
       enddo
    enddo
    
    return
  end subroutine buoyflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sgsvar(F_thl, F_tve, F_qw, F_qc, F_frac, F_pri, F_zn, F_ze,  &
       F_sigt, F_acoef, F_bcoef, &
       F_ccoef, F_vcoef, F_pblsigs, F_pblq1, F_ni, F_nkm1)
    use tdpack_const
    use phy_options
    use pbl_utils, only: dvrtdf
    use vintphy, only: vint_thermo2mom2
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_thl           !liquid water potential temperature (K; theta_l)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tve           !virtual temperature on e-levs (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw            !total water mixing ratio (kg/kg; q_tot)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_pri           !inverse Prandtl number
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zn            !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_ze            !dissipation length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermodynamic levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_acoef         !thermodynamic coefficient A
    real, dimension(F_ni,F_nkm1), intent(in) :: F_bcoef         !thermodynamic coefficient B
    real, dimension(F_ni,F_nkm1), intent(in) :: F_ccoef         !thermodynamic coefficient C
    real, dimension(*), intent(in) :: F_vcoef                   !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(out) :: F_qc           !PBL cloud water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_frac         !cloud fraction
    real, dimension(F_ni,F_nkm1), intent(out) :: F_pblsigs      !Subgrid moisture variance
    real, dimension(F_ni,F_nkm1), intent(out) :: F_pblq1        !Normalized saturation deficit

    !@Object Calculate the boundary layer subgrid-scale properties

    ! Local parameters
    real, parameter :: EPS=1e-10

    ! Local variables
    integer :: i,k
    real :: lnsig, c1coef
    real, dimension(F_ni,F_nkm1) :: dz, dqwdz, dthldz, sigmas, q1, thlm, qwm

    ! Pre-compute constants and vertical derivatives
    call vint_thermo2mom2(thlm,  qwm,  &
         &                F_thl, F_qw, F_vcoef, F_ni, F_nkm1)
    do k=1,F_nkm1-1
       do i=1,F_ni
          lnsig   = log(F_sigt(i,k+1)/F_sigt(i,k))
          dz(i,k) = -RGASD*F_tve(i,k)*lnsig/GRAV
       end do
    end do
    dz(:,F_nkm1) = 0.0
    call dvrtdf(dthldz, thlm, dz, F_ni, F_ni, F_ni, F_nkm1)
    call dvrtdf(dqwdz , qwm , dz, F_ni, F_ni, F_ni, F_nkm1)

    ! Compute subgrid standard deviation of supersaturation (s) on e-levels following Eq. 10 of BCMT95
    sigmas = 0.
    do k=1,F_nkm1-1
       do i=1,F_ni
          c1coef = sqrt(2. * PBL_RI_CK * F_pri(i,k) / PBL_RI_CAB * F_zn(i,k) * F_ze(i,k))
          sigmas(i,k) = max(c1coef * abs(F_acoef(i,k)*dqwdz(i,k) - F_bcoef(i,k)*dthldz(i,k)), EPS)
       end do
    end do

    ! Compute the normalized saturation deficit (Q1)
    q1(:,1:F_nkm1-1) = F_ccoef(:,1:F_nkm1-1) / sigmas(:,1:F_nkm1-1)
    q1(:,1) = 0. ; q1(:,F_nkm1) = 0. ;
    sigmas(:,1) = 0.; sigmas(:,F_nkm1) = 0.
    F_pblq1 = q1(:,:)
    F_pblsigs = sigmas(:,:)

    ! Compute cloud properties for local scalings
    do k=2,F_nkm1-1
       do i=1,F_ni

          ! Compute cloud fraction
          F_frac(i,k) = 0.5 * (1. + erf(q1(i,k)/sqrt(2.)))

          ! Compute condensate (note factor 2 difference in sigmas in BFP92 vs BCMT95)
          F_qc(i,k) = sigmas(i,k) * (F_frac(i,k)*q1(i,k) + exp((-q1(i,k)**2/2.)/sqrt(2.*PI)))

       end do
    end do

    ! Fill array extrema
    F_frac(1:F_ni,1) = 0. ; F_frac(1:F_ni,F_nkm1) = 0.
    F_qc  (1:F_ni,1) = 0. ; F_qc  (1:F_ni,F_nkm1) = 0.

    return
  end subroutine sgsvar
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thermco(F_thl, F_qw, F_tif, F_fn, F_sigt, F_ps, F_type, F_ni, F_nkm1, &
       F_acoef, F_bcoef, F_ccoef, F_leff)
    use tdpack
    use pbl_utils, only: ficemxp
    implicit none

    ! Arguments
    integer, intent(in) :: F_ni                                         !horizontal dimension
    integer, intent(in) :: F_nkm1                                       !vertical dimension
    integer, intent(in) :: F_type                                       !cloud type switch (0=implicit; 1=explicit; 2=both)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_thl                   !liquid water potential temperature (K; theta_l)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw                    !total water content (kg/kg; q_tot)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tif                   !temperature used for ice fraction (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fn                    !cloud fraction    
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt                  !sigma of thermodynamic levels
    real, dimension(F_ni), intent(in) :: F_ps                           !surface pressure (Pa)
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_acoef      !thermodynamic coefficient A
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_bcoef      !thermodynamic coefficient B
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_ccoef      !thermodynamic coefficient C
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_leff       !effective latent heat

    !Object
    !          Calculate the thermodynamic coefficients used in the presence of clouds
    !          and the conservative variables.

    !Notes
    !          See definitions in:
    !          - Bechtold and Siebesma 1998, JAS 55, 888-895

    ! Local variable declarations
    integer :: i,k
    real, dimension(F_ni,F_nkm1) :: pres, exner, qsat, dqsat, tl, tfice, dfice, work, &
         acoef, bcoef, ccoef, fice, leff

    ! Diagnose implicit ice fraction
    call ficemxp(fice, tfice, dfice, F_tif, F_ni, F_nkm1)
    
    ! Precompute pressure and Exner function
    do k=1,F_nkm1
       pres(:,k) = F_sigt(:,k) * F_ps(:)
       exner(:,k) = F_sigt(:,k)**CAPPA       
    enddo

    ! Compute liquid water temperature
    tl(:,:) = exner(:,:)*F_thl(:,:)

    ! Effective latent heat release
    leff(:,:) = CHLC + fice(:,:)*CHLF
    if (present(F_leff)) F_leff(:,:) = leff(:,:)
    
    ! Compute thermodynamic coefficients following Bechtold and Siebesma (JAS 1998) Appendix A
    COEFS: if (present(F_acoef) .or. present(F_bcoef) .or. present(F_ccoef)) then

       ! Compute saturation specific humidity for selected cloud type
       CLOUD_TYPE: select case (F_type)

       case(IMPLICIT_CLOUD) !Implicit clouds
          call ficemxp(work, tfice, dfice, F_tif, F_ni, F_nkm1)
          do k=1,F_nkm1
             do i=1,F_ni                
                qsat(i,k) = fqsmx(tl(i,k), pres(i,k), tfice(i,k))
             enddo
          enddo
          call mfdlesmx(work, tl, tfice, dfice, F_ni, F_nkm1)
          do k=1,F_nkm1
             do i=1,F_ni
                dqsat(i,k) = fdqsmx(qsat(i,k), work(i,k) )
             enddo
          enddo

       case (EXPLICIT_CLOUD) !Explicit clouds
          do k=1,F_nkm1
             do i=1,F_ni
                qsat(i,k) = foqsa(tl(i,k), pres(i,k))
                dqsat(i,k) = fodqa(qsat(i,k), tl(i,k))
             enddo
          enddo

       case (COMBINED_CLOUD) !Combined implicit and explicit clouds
          call ficemxp(work, tfice, dfice, F_tif, F_ni, F_nkm1)
          call mfdlesmx(work, tl, tfice, dfice, F_ni, F_nkm1)
          do k=1,F_nkm1
             do i=1,F_ni
                if (F_fn(i,k) < 1.0) then
                   qsat(i,k) = fqsmx(tl(i,k), pres(i,k), tfice(i,k))
                   dqsat(i,k) = fdqsmx(qsat(i,k), work(i,k) )
                else
                   qsat(i,k) = foqsa(tl(i,k), pres(i,k))
                   dqsat(i,k) = fodqa(qsat(i,k), tl(i,k))
                endif
             enddo
          enddo

       end select CLOUD_TYPE

       ! Comopute buoyancy flux coefficients
       do k=1,F_nkm1
          do i=1,F_ni             
             acoef(i,k) = 1.0 / ( 1.0 + (leff(i,k)/CPD)*dqsat(i,k))
             bcoef(i,k) = acoef(i,k)*exner(i,k)*dqsat(i,k)
             ccoef(i,k) = acoef(i,k)*(F_qw(i,k) - qsat(i,k))
          enddo
       enddo
       if (present(F_acoef)) F_acoef(:,:) = acoef(:,:)
       if (present(F_bcoef)) F_bcoef(:,:) = bcoef(:,:)
       if (present(F_ccoef)) F_ccoef(:,:) = ccoef(:,:)
    endif COEFS

    return
  end subroutine thermco

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
    !          Update TKE using the analytic part of the TKE equation
    
    ! Local variable declarations
    integer :: i, k
    complex :: ye, itbc, xi0, gamma, yz
    real :: sgn, tmax, bsafe

    ! Update TKE using an anlytic solution of the non-diffusive TKE equation
    do k=1,F_nkm1
       do i=1,F_ni
          bsafe = sign(max(abs(F_b(i,k)), etrmin2/F_tau), F_b(i,k))
          ye = sqrt(cmplx(bsafe / F_c(i,k)))
          itbc = 0.5 * sqrt(cmplx(bsafe * F_c(i,k)))  
          xi0 = sqrt(F_en(i,k)) / ye
          sgn = sign(1., 1. - xi0%RE)
          gamma = atanh(xi0**sgn)
          tmax = -min(gamma%IM, -tiny(gamma%IM)*F_tau) / max(itbc%IM, tiny(itbc%IM))
          yz = ye * tanh(gamma + itbc * min(F_tau, tmax))**sgn
          F_estar(i,k) = max(yz%RE**2, etrmin2)
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine tkestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tkebudget(F_en, F_estar, F_b, F_c, F_zn, F_dvdz2, F_buoy_flux, &
       F_shr, F_buoy, F_diss, F_tau, F_ni, F_nkm1)
    implicit none

    ! Argument declarations
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_en            !current turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_estar         !updated turbulent kinetic energy (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_b             !B (generation) term of TKE equation (m/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_c             !C (dissipation) term of TKE eqaution (1/m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zn            !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dvdz2         !square of vertical shear (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_buoy_flux     !buoyancy flux (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_shr          !shear generation term (m2/s3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_buoy         !buoyancy generation/suppression term (m2/s3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_diss         !viscous dissipation term (m2/s3)

    !Object
    !          Compute budget terms for TKE equation
    
    ! Local variable declarations
    integer :: i, k
    real :: clambda, bterm, b_over_c
          
    ! Compute TKE budget diagnostics
    do k=1,F_nkm1
       do i=1,F_ni
          b_over_c = F_b(i,k) / F_c(i,k)
          if( F_en(i,k) /= b_over_c .and. F_estar(i,k) /= b_over_c ) then
             bterm = log( abs( (F_estar(i,k) - b_over_c) / (F_en(i,k) - b_over_c) ) )
          else
             bterm = 0.
          endif
          clambda = PBL_RI_CK * F_zn(i,k) * bterm / (F_c(i,k) * F_tau)
          F_shr(i,k) = -clambda * F_dvdz2(i,k)
          F_buoy(i,k) = clambda * F_buoy_flux(i,k)
          F_diss(i,k) = (F_estar(i,k) - F_en(i,k)) / F_tau - (F_shr(i,k) + F_buoy(i,k))
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine tkebudget
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_conserved(F_tcon, F_qcon, F_tt, F_hu, F_qc, F_sigt, F_ni, F_nkm1, F_qi)
    use tdpack, only: CHLC, CHLF, CPD, CAPPA
    use pbl_utils, only: ficemxp
    implicit none

    ! Argument declarations
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_hu            !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qc            !cloud water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma of thermodynamic levels
    real, dimension(F_ni,F_nkm1), intent(in), optional :: F_qi  !cloud ice content (kg/kg) [diagnosed]
    real, dimension(F_ni,F_nkm1), intent(out) :: F_tcon         !conserved heat variable (theta_l; K)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_qcon         !conserved moisture variable (q_tot; kg/kg)

    !Object
    !          Compute conserved variables

    ! Local parameter definitions
    integer :: i,k
    real, dimension(F_ni,F_nkm1) :: fice, unused, condens
    
    ! Handle optional arguments
    if (present(F_qi)) then
       do k=1,F_nkm1
          do i=1,F_ni
             condens(i,k) = max(F_qc(i,k) + F_qi(i,k), 0.)
             fice(i,k) = F_qi(i,k) / max(condens(i,k), tiny(condens))
          enddo
       enddo
    else
       call ficemxp(fice, unused, unused, F_tt, F_ni, F_nkm1)
       condens(:,:) = max(F_qc(:,:), 0.)
    endif

    ! Compute conserved variables
    do k=1,F_nkm1
       do i=1,F_ni
          F_qcon(i,k) = F_hu(i,k) + condens(i,k)          
          F_tcon(i,k) = (F_sigt(i,k)**(-CAPPA)) * (F_tt(i,k) - &
               ((CHLC + fice(i,k)*CHLF) / CPD) * condens(i,k))
       enddo
    enddo

  end subroutine compute_conserved

end module pbl_ri_utils
