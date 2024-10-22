module cons_thlqw
  use phy_status, only: PHY_OK
  implicit none
  private

  ! API functions
  public :: thlqw_compute                                       !Compute conserved variables
  public :: thlqw_thermco                                       !Compute thermodynamic coefficients
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thlqw_compute(F_tcon, F_qcon, F_tt, F_hu, F_lwc, F_iwc, F_sigt, &
       F_ni, F_nkm1)
    use tdpack, only: CHLC, CHLF, CPD, CAPPA
    implicit none

    ! Argument declarations
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_hu            !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc           !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc           !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma of thermodynamic levels
    real, dimension(F_ni,F_nkm1), intent(out) :: F_tcon         !conserved heat variable (theta_l; K)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_qcon         !conserved moisture variable (q_tot; kg/kg)

    !Object Compute theta_l and qw conserved variables

    ! Local variables
    integer :: i,k
    real :: twc
    
    ! Compute conserved variables
    do k=1,F_nkm1
       do i=1,F_ni
          twc = F_lwc(i,k) + F_iwc(i,k)
          F_qcon(i,k) = F_hu(i,k) + twc  
          F_tcon(i,k) = (F_sigt(i,k)**(-CAPPA)) * (F_tt(i,k) - &
               (CHLC * twc + CHLF * F_iwc(i,k)) / CPD)
       enddo
    enddo
    
    ! End of subprogram
    return
  end subroutine thlqw_compute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thlqw_thermco(F_thl, F_qw, F_tt, F_lwc, F_iwc, F_sigt, F_ps, &
       F_ni, F_nkm1, F_acoef, F_bcoef, F_ccoef, F_leff, F_rhc)
    use tdpack
    use microphy_utils, only: mp_icefrac
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                         !horizontal dimension
    integer, intent(in) :: F_nkm1                                       !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_thl                   !liquid water potential temperature (K; theta_l)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw                    !total water content (kg/kg; q_tot)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt                    !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc                   !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc                   !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt                  !sigma of thermodynamic levels
    real, dimension(F_ni), intent(in) :: F_ps                           !surface pressure (Pa)
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_acoef      !thermodynamic coefficient A
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_bcoef      !thermodynamic coefficient B
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_ccoef      !thermodynamic coefficient C
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_leff       !effective latent heat (J/kg)
    real, dimension(F_ni,F_nkm1), intent(in), optional :: F_rhc         !critical relative humidity [1.]
    integer :: F_istat                                                  !return status of function

    !@Object  Calculate the thermodynamic coefficients used in the presence of clouds
    !         following Bechtold and Siebesma (1998; JAS 888-895)

    ! Local variable declarations
    integer :: i,k
    real, dimension(F_ni,F_nkm1) :: pres, exner, qsat, dqsat, tl, dfice, esat, &
         acoef, bcoef, ccoef, fice, leff, rhc
    
    ! Set optional values
    rhc(:,:) = 1.
    if (present(F_rhc)) rhc(:,:) = F_rhc(:,:)
    
    ! Diagnose implicit ice fraction
    if (mp_icefrac(fice, F_tt, F_lwc, F_iwc, F_ni, F_nkm1) /= PHY_OK) then
       call physeterror('cons_thlqw::thlqw_thermco', 'Failed to compute ice fraction')
       return
    endif
    
    ! Precompute pressure and Exner function
    do k=1,F_nkm1
       pres(:,k) = F_sigt(:,k) * F_ps(:)
       exner(:,k) = F_sigt(:,k)**CAPPA       
    enddo

    ! Compute liquid water temperature
    tl(:,:) = max(exner(:,:)*F_thl(:,:), 100.)

    ! Effective latent heat release
    leff(:,:) = CHLC + fice(:,:)*CHLF
    if (present(F_leff)) F_leff(:,:) = leff(:,:)
    
    ! Compute thermodynamic coefficients following Bechtold and Siebesma (JAS 1998) Appendix A
    COEFS: if (present(F_acoef) .or. present(F_bcoef) .or. present(F_ccoef)) then

       ! Compute saturation specific humidity and its derivative wrt temperature
       if (mp_icefrac(fice, F_tt, F_lwc, F_iwc, F_ni, F_nkm1, F_difdt=dfice) /= PHY_OK) then
          call physeterror('cons_thlqw::thlqw_thermco', 'Failed to compute ice fraction delta')
          return
       endif
       do k=1,F_nkm1
          do i=1,F_ni             
             qsat(i,k) = rhc(i,k) * fqsmx(tl(i,k), pres(i,k), fice(i,k))
          enddo
       enddo
       call mfdlesmx(esat, tl, fice, dfice, F_ni, F_nkm1)
       do k=1,F_nkm1
          do i=1,F_ni
             dqsat(i,k) = fdqsmx(qsat(i,k), rhc(i,k)*esat(i,k))
          enddo
       enddo
       
       ! Compute buoyancy flux coefficients
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

    ! End of subprogram
    return
  end subroutine thlqw_thermco

end module cons_thlqw
