module microphy_statcond
  use phy_status, only: PHY_OK, PHY_ERROR
  implicit none
  private

  ! API functions
  public :: sc_adjust                                           !Statistical condensation adjustment
  public :: sc_condense                                         !Compute condensation tendencies

  ! Generic functions
  interface sc_condense
     module procedure sc_condcons
     module procedure sc_condonly
  end interface sc_condense

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sc_adjust(F_dtt, F_dhu, F_dqc, F_pvars, F_tau, F_ni, F_nkm1)
    use phybusidx
    use phymem, only: phyvar
    use tendency, only: apply_tendencies
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step (s)
    type(phyvar), dimension(:), pointer, contiguous :: F_pvars  !physics variables
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dtt          !air temperature tendency (K/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dhu          !specific humidity tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dqc          !condensate tendency (kg/kg/s)

    !@Objective Adjust the profile to account for condensation processes

    ! Common declarations
#include "phymkptr.hf"

    !@Local variables
    real, dimension(:), pointer, contiguous :: ztdmaskxdt
    real, dimension(F_ni,F_nkm1) :: tqc, tqi
    real, dimension(:,:), pointer, contiguous :: ztplus, zhuplus, zqcplus

    ! Assign bus pointers
    MKPTR1D(ztdmaskxdt, tdmaskxdt, F_pvars)
    MKPTR2Dm1(zhuplus, huplus, F_pvars)
    MKPTR2Dm1(zqcplus, qcplus, F_pvars)
    MKPTR2Dm1(ztplus, tplus, F_pvars)
    
    ! Compute condensation
    call sc_condense(F_dtt, F_dhu, tqc, tqi, F_pvars, F_tau, F_ni, F_nkm1)

    ! Apply condensation tendencies
    call apply_tendencies(ztplus,  zhuplus, &
         &                F_dtt, F_dhu, ztdmaskxdt, F_ni, F_nkm1+1, F_nkm1)
    if (associated(zqcplus)) then
       F_dqc(:,1:F_nkm1) = tqc(:,1:F_nkm1) + tqi(:,1:F_nkm1)  !for S2 only       
       call apply_tendencies(zqcplus, F_dqc, ztdmaskxdt, F_ni, F_nkm1+1, F_nkm1)       
    endif    

    ! End of subprogram
    return
  end subroutine sc_adjust
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sc_condonly(F_dtt, F_dhu, F_dqc, F_dqi, F_pvars, F_tau, F_ni, F_nkm1)
    use phymem, only: phyvar
    use phybusidx, except1=>lwc, except2=>iwc
    use cons_thlqw, only: thlqw_compute
    use microphy_utils, only: mp_lwc, mp_iwc
    use phy_status, only: PHY_OK
    implicit none
    
    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step (s)
    type(phyvar), dimension(:), pointer, contiguous :: F_pvars  !physics variables
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dtt          !air temperature tendency (K/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dhu          !specific humidity tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dqc          !liquid condensate tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dqi          !ice condensate tendency (kg/kg/s)
    integer :: F_istat                                          !return status of function

    !@Objective Compute tendencies resulting from cloud condensation/evaporation
    !           from tendencies of model state variables.
    
    ! Common declarations
#include "phymkptr.hf"

    !@Local variables
    integer :: i, k
    real, dimension(F_ni,F_nkm1) :: thl, qw, lwc, iwc, zero
    real, dimension(:,:), pointer, contiguous :: ztplus, zhuplus, zsigt

    ! Assign bus pointers
    MKPTR2Dm1(zhuplus, huplus, F_pvars)
    MKPTR2Dm1(zsigt, sigt, F_pvars)
    MKPTR2Dm1(ztplus, tplus, F_pvars)

    ! Compute condensate masses
    if (mp_lwc(lwc, F_pvars) /= PHY_OK .or. &
         mp_iwc(iwc, F_pvars) /= PHY_OK) then
       call physeterror('microphy_statcond::sc_condonly', &
            'Cannot compute condensate contents')
       return
    endif
    
    ! Compute conserved variables and tendencies
    call thlqw_compute(thl, qw, ztplus, zhuplus, lwc, iwc, zsigt, &
         F_ni, F_nkm1)

    ! Compute condensation
    zero(:,:) = 0.
    call sc_condcons(F_dtt, F_dhu, F_dqc, F_dqi, thl, zero, &
         qw, zero, F_pvars, F_tau, F_ni, F_nkm1, F_lwc=lwc, F_iwc=iwc)    

    ! End of subprogram
    return
  end subroutine sc_condonly
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sc_condcons(F_dtt, F_dhu, F_dqc, F_dqi, F_thl, F_dthl, &
       F_qw, F_dqw, F_pvars, F_tau, F_ni, F_nkm1, F_lwc, F_iwc)
    use phymem, only: phyvar
    use phybusidx, except1=>lwc, except2=>iwc
    use microphy_utils, only: mp_lwc, mp_iwc
    use tdpack_const, only: CHLC, CHLF, CAPPA, CPD
    implicit none
    
    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_thl           !liquid water potential temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dthl          !tendency of theta_l (K/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw            !specific total water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_dqw           !tendency of q_w (kg/kg/s)
    type(phyvar), dimension(:), pointer, contiguous :: F_pvars  !physics variables
    real, dimension(F_ni,F_nkm1), intent(in), optional :: F_lwc !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in), optional :: F_iwc !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dtt          !air temperature tendency (K/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dhu          !specific humidity tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dqc          !liquid condensate tendency (kg/kg/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dqi          !ice condensate tendency (kg/kg/s)
    integer :: F_istat                                          !return status of function

    !@Objective Compute tendencies resulting from cloud condensation/evaporation
    !           in conjunction with changes to conserved variables
    
    ! Common declarations
#include "phymkptr.hf"
    
    !@Local variables
    integer :: i,k
    real :: qv, twc, twcp, qc, qi, dqc, dqi
    real, dimension(:), pointer, contiguous :: zpmoins
    real, dimension(F_ni,F_nkm1) :: qcp, qip, thlp, qwp, &
         lwc, iwc
    real, dimension(:,:), pointer, contiguous :: zsigmas, ztplus, &
         zsigt, zfxp
    
    ! Assign bus pointers
    MKPTR1D(zpmoins, pmoins, F_pvars)
    MKPTR2Dm1(zfxp, fxp, F_pvars)
    MKPTR2Dm1(zsigmas, sigmas, F_pvars)
    MKPTR2Dm1(zsigt, sigt, F_pvars)
    MKPTR2Dm1(ztplus, tplus, F_pvars)

    ! Compute condensate masses
    if (present(F_lwc)) then
       lwc(:,:) = F_lwc(:,:)
    else
       if (mp_lwc(lwc, F_pvars) /= PHY_OK) then
          call physeterror('microphy_statcond::sc_condcons', &
               'Cannot compute liquid condensate contents')
          return
       endif
    endif
    if (present(F_iwc)) then
       iwc(:,:) = F_iwc(:,:)
    else
       if (mp_iwc(iwc, F_pvars) /= PHY_OK) then
          call physeterror('microphy_statcond::sc_condcons', &
               'Cannot compute ice condensate contents')
          return
       endif
    endif
    
    ! Diagnose updated cloud state
    do k=1,F_nkm1
       do i=1,F_ni
          thlp(i,k) = F_thl(i,k) + F_tau * F_dthl(i,k)          
          qwp(i,k) = F_qw(i,k) + F_tau * F_dqw(i,k)
       enddo
    enddo
    if (sc_diag(qcp, qip, zfxp, thlp, qwp, ztplus, zsigmas, lwc, iwc, &
         zsigt, zpmoins, F_ni, F_nkm1) /= PHY_OK) then
       call physeterror('microphy_statcond::sc_condcons', &
            'Failed to diagnose current cloud state')
       return
    endif

    ! Compute state variable tendencies
    do k=1,F_nkm1
       do i=1,F_ni
          qc = max(lwc(i,k), 0.)
          qi = max(iwc(i,k), 0.)
          twc = qc + qi
          twcp = qcp(i,k) + qip(i,k)
          qv = max(F_qw(i,k) - twc, 0.)
          dqc = (qcp(i,k) - qc) / F_tau          
          dqc = max(dqc, -qc/F_tau)                                  !prevent negative qc
          F_dqc(i,k) = min(dqc, F_dqw(i,k) + qv / F_tau)             !prevent over-depletion of vapour
          dqi = (qip(i,k) - qi) / F_tau
          dqi = max(dqi, -qi/F_tau)                                  !prevent negative qi
          F_dqi(i,k) = min(dqi, F_dqw(i,k) + qv / F_tau)             !prevent over-depletion of vapour
          F_dtt(i,k) = F_dthl(i,k) * zsigt(i,k)**CAPPA &
               + CHLC/CPD * F_dqc(i,k) + (CHLC+CHLF)/CPD * F_dqi(i,k)
          F_dhu(i,k) = max(F_dqw(i,k) - (F_dqc(i,k)+F_dqi(i,k)), -qv / F_tau)
       enddo
    enddo
    
    ! End of subprogram
    return
  end subroutine sc_condcons

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sc_diag(F_qc, F_qi, F_fn, F_thl, F_qw, F_tt, F_sigmas, F_lwc, F_iwc, &
       F_sigt, F_ps, F_ni, F_nkm1, F_icefrac) result(F_istat)
    use cons_thlqw, only: thlqw_thermco
    use microphy_utils, only: mp_icefrac
    use pbl_ri_utils, only: PBL_SDMIN
    use tdpack_const, only: CAPPA, PI
    use tdpack, only: fqsmx
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_thl           !liquid water potential temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qw            !specific total water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigmas        !subgrid moisture variance (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc           !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc           !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma of thermodynamic levels
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_fn         !cloud fraction
    real, dimension(F_ni,F_nkm1), intent(out) :: F_qc           !specific liquid condensate (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_qi           !specific ice condensate (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_icefrac !ice fraction
    integer :: F_istat                                          !return status of function

    !@Object Diagnose statistical clouds from conserved variables

    !@Local variables
    integer :: i,k
    real :: pres, tl, qsat, qliq
    real, dimension(F_ni,F_nkm1) :: acoef, q1, fice, sigmas

    ! Initialize return value
    F_istat = PHY_ERROR

    ! Compute conserved-variable constants
    call thlqw_thermco(F_thl, F_qw, F_tt, F_lwc, F_iwc, F_sigt, F_ps, &
         F_ni, F_nkm1, F_acoef=acoef)

    ! Diagnose ice fraction
    if (mp_icefrac(fice, F_tt, F_lwc, F_iwc, F_ni, F_nkm1) /= PHY_OK) then
       call physeterror('microphy_statcond::sc_diag', &
            'Failed to compute ice fraction')
       return
    endif
    
    ! Compute Q1
    do k=1,F_nkm1
       do i=1,F_ni
          pres = F_sigt(i,k) * F_ps(i)
          tl = max(F_thl(i,k) * F_sigt(i,k)**CAPPA, 100.)
          qsat = fqsmx(tl, pres, fice(i,k))
          sigmas(i,k) = max(F_sigmas(i,k), PBL_SDMIN)
          q1(i,k) = acoef(i,k) * (F_qw(i,k) - qsat) / sigmas(i,k)
       enddo
    enddo

    ! Compute cloud properties based on a Gaussian SGS distribution
    do k=1,F_nkm1
       do i=1,F_ni
          F_fn(i,k) = 0.5 * (1. + erf(q1(i,k)/sqrt(2.)))
          qliq = sigmas(i,k) * (F_fn(i,k) * q1(i,k) + &
               exp(-0.5*q1(i,k)**2)/sqrt(2.*PI))
          F_qi(i,k) = qliq * fice(i,k)
          F_qc(i,k) = qliq - F_qi(i,k)
       enddo
    enddo

    ! Fill optional return values
    if (present(F_icefrac)) F_icefrac(:,:) = fice(:,:)
    
    ! Successful completion
    F_istat = PHY_OK
    return
  end function sc_diag

end module microphy_statcond
