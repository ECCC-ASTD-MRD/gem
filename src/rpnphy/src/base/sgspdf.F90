module sgspdf
  implicit none
  private

  public :: sgspdf_bkg                  !Compute background component of subgrid-scale variance
  public :: sgspdf_compute              !Compute std dev of subgrid-scale saturation deficit

contains

  subroutine sgspdf_compute(F_sigmas, F_tt, F_hu, F_lwc, F_iwc, F_pri, &
       F_zn, F_zd, F_gzt, F_sigt, F_ps, F_hpbl, &
       F_dxdy, F_mrk, F_vcoef, F_ni, F_nkm1)
    use tdpack_const, only: CAPPA
    use tdpack, only: fqsmx
    use debug_mod, only: init2nan
    use cons_thlqw, only: thlqw_compute
    use pbl_utils, only: dvrtdf
    use pbl_ri_utils, only: pblri_sgspdf
    use microphy_utils, only: mp_icefrac
    use phy_status, only: PHY_OK
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_hu            !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc           !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc           !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_pri           !inverse Prandtl number
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zn            !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_zd            !dissipation length (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzt           !heights of thermo levels (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(F_ni), intent(in) :: F_hpbl                 !boundary layer depth (m)      
    real, dimension(F_ni), intent(in) :: F_dxdy                 !grid cell area (m2)
    real, dimension(F_ni), intent(in) :: F_mrk                  !Markov chains for stochastic parameters
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(out) :: F_sigmas       !subgrid moisture std dev (kg/kg)

    !@Object Calculate the standard deviation of the subgrid-scale saturation deficit

    ! Local parameters
    real, parameter :: REF_DX=500.                              !grid mesh for variance scaling

    ! Local variables
    integer :: i,k
    real :: c1coef, idz, qsat, tl, pres
    real, dimension(F_ni,F_nkm1) :: thl, qw, dqwdz, dthldz, acoef, bcoef, ccoef, &
         sigs_pbl, sigs_bkg, fice

    ! Basic initializations
    call init2nan(sigs_pbl, sigs_bkg)
    
    ! Compute conserved variables
    call thlqw_compute(thl, qw, F_tt, F_hu, F_lwc, F_iwc, F_sigt, F_ni, F_nkm1)
    
    ! Estimate turbulence component
    call pblri_sgspdf(sigs_pbl, thl, qw, F_lwc, F_iwc, F_tt, F_pri, &
         F_zn, F_zd, F_gzt, F_sigt, F_ps, F_dxdy, F_vcoef, F_ni, F_nkm1)

    ! Estimate background component
    call sgspdf_bkg(sigs_bkg, F_tt, qw, F_sigt, F_gzt, F_ps, F_hpbl, F_mrk, F_ni, F_nkm1)
    
    ! Combine for estimate of subgrid-scale variability
    do k=1,F_nkm1
       do i=1,F_ni
          F_sigmas(i,k) = sqrt(sigs_pbl(i,k)**2 + sigs_bkg(i,k)**2)          
       enddo
    enddo

    ! Impose physical limit on subgrid-scale variability assuming T'=p'=0
    if (mp_icefrac(fice, F_tt, F_lwc, F_iwc, F_ni, F_nkm1) /= PHY_OK) then
       call physeterror('sgspdf_compute', 'Cannot compute ice fraction')
       return
    endif
    do k=1,F_nkm1
       do i=1,F_ni
          pres = F_sigt(i,k) * F_ps(i)
          tl = max(thl(i,k) * F_sigt(i,k)**CAPPA, 100.)
          qsat = fqsmx(tl, pres, fice(i,k))
          F_sigmas(i,k) = min(F_sigmas(i,k), qsat / (2.*sqrt(6.)))
       enddo
    enddo
    
    return
  end subroutine sgspdf_compute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sgspdf_bkg(F_sigmas, F_tt, F_qw, F_sigt, F_gzt, F_ps, F_hpbl, F_mrk, F_ni, F_nkm1) 
    use phy_options, only: cond_hu0min, cond_hu0max
    use ens_perturb, only: ens_nc2d, ens_spp_get
    use tdpack, only: foqst

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw            !total water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzt           !height for thermo levels (m)
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(F_ni), intent(in) :: F_hpbl                 !boundary layer depth (m)      
    real, dimension(F_ni,ens_nc2d), intent(in) :: F_mrk         !Markov chains for stochastic parameters
    real, dimension(F_ni,F_nkm1), intent(out) :: F_sigmas       !subgrid moisture std dev (kg/kg)

    !@Object Calculate a background subgrid-scale variance

    ! Local variables
    integer :: i, k
    real :: rhc, tfact, tadj, qsat
    real, dimension(F_ni) :: hu0min, hu0max, dzslope, huslope

    ! Compute profile of RH-based SGS std dev
    hu0min(:) = ens_spp_get('hu0min', F_mrk, default=cond_hu0min)
    hu0max(:) = ens_spp_get('hu0max', F_mrk, default=cond_hu0max)
    dzslope(:) = max(F_hpbl(:), 2000.)  !slope from h to 2*h (min 2000m)
    huslope(:) = -(hu0max(:) - hu0min(:)) / dzslope
    do k=1,F_nkm1
       do i=1,F_ni
          qsat = foqst(F_tt(i,k), F_sigt(i,k)*F_ps(i))
          ! Set basic profile for rhcrit
          rhc = hu0max(i) + huslope(i) * (F_gzt(i,k) - F_hpbl(i))
          rhc = max(hu0min(i), min(hu0max(i), rhc))
          ! Adjust rhcrit for cold temperatures
          tfact = 1. + 0.15 * max(238. - F_tt(i,k), 0.)
          tadj = (hu0max(i) - rhc) * (1. - 1./tfact)
          rhc = min(rhc + tadj, hu0max(i))
          ! Limit rhcrit in very dry conditions
          rhc = max(rhc, 1. - F_qw(i,k)/qsat)
          ! Convert to sigma_s for a triangular PDF
          F_sigmas(i,k) = qsat * (1.-rhc) / sqrt(6.) 
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine sgspdf_bkg
    
end module sgspdf
