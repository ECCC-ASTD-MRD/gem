module pbl_ri
  implicit none
  private
  
  ! Public API functions
  public :: rpnint

contains

  subroutine rpnint(F_en, F_km, F_kt, F_pri, F_rif, F_rig, F_buoy, F_shr2, &
       F_buoyen, F_shren, F_diffen, F_dissen, F_qwvar, F_zn, F_zd, &
       F_enold, F_uu, F_vv, F_tt, F_hu, F_lwc, F_iwc, F_fbl, F_turbreg, &
       F_hpbl, F_lh, F_ps, F_z0m, F_z0t, F_frv, F_wstar, F_qstar, &
       F_hpar, F_tsurf, F_qsurf, F_lat, F_fcor, &
       F_sigm, F_sigt, F_gzm, F_gzt, F_dxdy, F_mrk2, F_vcoef, &
       F_tau, F_kount, F_ni, F_nkm1)

    use, intrinsic :: iso_fortran_env, only: INT64
    use tdpack, only: GRAV, KARMAN, RGASD
    use phy_options, only: ilongmel, etrmin2, phystepoutlist_s, nphystepoutlist, &
         pbl_progvar, pbl_ae, pbl_turbsl_depth, pbl_ricrit, debug_alldiag_L
    use phy_status, only: phy_error_L, PHY_OK
    use mixing_length, only: ml_compute, ML_LMDA
    use ens_perturb, only: ens_nc2d
    use pbl_ri_utils, only: buoyflux, tkestep, tkebudget, kcalc, covareq, covarstep, &
         ktosig, K_M, K_T, PBL_RI_CK, PBL_RI_CU, PBL_RI_CW, &
         PBL_RI_CE, PBL_RI_CAB
    use pbl_ri_diffuse, only: diffuse, FIELD_IS_TKE, FIELD_ON_TLEV
    use pbl_utils, only: LAMINAR, TURBULENT
    use cons_thlqw, only: thlqw_compute
    implicit none
#include <rmnlib_basics.hf>
#include "phymkptr.hf"

    !Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    integer, intent(in) :: F_kount                              !time step number
    real, intent(in) :: F_tau                                   !time step length (s)
    real, dimension(F_ni), intent(in) :: F_hpbl                 !height of the PBL (m)
    real, dimension(F_ni), intent(in) :: F_lh                   !launching height (m)
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(F_ni), intent(in) :: F_z0m                  !momentum roughness length (m)
    real, dimension(F_ni), intent(in) :: F_z0t                  !thermal roughness length (m)
    real, dimension(F_ni), intent(in) :: F_frv                  !friction velocity (m/s)
    real, dimension(F_ni), intent(in) :: F_wstar                !convective velocity scale (m/s)
    real, dimension(F_ni), intent(in) :: F_qstar                !moisture scale (kg/kg)
    real, dimension(F_ni), intent(in) :: F_tsurf                !Surface air temperature (K)
    real, dimension(F_ni), intent(in) :: F_qsurf                !Surface specific humidity (kg/kg)
    real, dimension(F_ni), intent(in) :: F_lat                  !Latitude (rad)
    real, dimension(F_ni), intent(in) :: F_fcor                 !Coriolis factor (/s)
    real, dimension(F_ni), intent(in) :: F_dxdy                 !horizontal grid area (m^2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_enold         !TKE of previous time step (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_hu            !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_uu            !u-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_vv            !v-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigm          !sigma for momenntum levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzm           !height of momentum levels (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_gzt           !height of thermo levels (m)
    real, dimension(F_ni,ens_nc2d), intent(in) :: F_mrk2        !Markov chains for stochastic parameters
    real, pointer, dimension(:,:,:), contiguous :: F_vcoef      !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc           !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc           !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fbl           !PBL cloud fraction
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_turbreg    !turbulence regime
    real, dimension(F_ni), intent(inout) :: F_hpar              !height of parcel ascent (m)
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_zn         !momentum mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_zd           !dissipation length (m)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_en           !updated TKE (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_km           !eddy diffusivity for momentum (m2/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_kt           !eddy diffusivity for heat/moisture (m2/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_pri          !inverse Prandtl number (ratio of KT/KM diffusion coefficients)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_rif          !flux Richardson number
    real, dimension(F_ni,F_nkm1), intent(out) :: F_rig          !gradient Richardson number
    real, dimension(F_ni,F_nkm1), intent(out) :: F_buoy         !buoyancy flux (m^2/s^3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_shr2         !square of vertical wind shear (s^-2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_buoyen       !buoyancy term of TKE budget (m^2/s^3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_shren        !shear term of TKE budget (m^2/s^3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_diffen       !diffusion term of TKE budget (m^2/s^3)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dissen       !dissipation term of TKE budget (m^2/s^3)
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_qwvar      !subgrid-scale moisture variance (kg/kg)^2

    !@Object
    !          Update TKE and diffusion coefficients for the RPN Integrated PBL scheme
    include "phyinput.inc"

    ! External symbols
    integer, external :: neark
    
    ! Local variable declarations
    integer :: stat, i, k
    integer, dimension(F_ni) :: mlen, slk
    real :: dae
    real, dimension(F_ni) :: tke_sfc, zero_sfc
    real, dimension(F_ni,F_nkm1) :: e_star, ke, b_term, c_term, tv, zero, one, qw, thl, &
         qwvarstar, ktsg, qwvardiff
    real, dimension(F_ni,F_nkm1,3) :: zero3d

    ! Initialization
    zero_sfc = 0.
    zero = 0.
    zero3d = 0.
    one = 1.
    if (F_kount == 0) then
       if (.not.ISPHYIN('zn')) then
          do k=1,F_nkm1
             F_zn(:,k) = min(KARMAN*(F_gzm(:,k) + F_z0m(:)), ML_LMDA)
          enddo
       endif
    endif

    ! Estimate buoyancy and shear gradients (not buoyancy flux until later)
    call buoyflux(F_buoy, F_shr2, F_rig, F_uu, F_vv, F_tt, F_hu, &
         F_lwc, F_iwc, F_fbl, F_tsurf, F_qsurf, F_z0m, F_z0t, &
         F_lat, F_fcor, F_sigt, F_gzm, F_gzt, F_ps, F_vcoef, F_ni, F_nkm1)

    ! Determine turbulence regime
    if (pbl_turbsl_depth > 0.) then
       stat = neark(F_sigt, F_ps, pbl_turbsl_depth, F_ni, F_nkm1, slk)
    else
       slk(:) = F_nkm1
    endif
    do k=1,F_nkm1
       do i=1,F_ni
          if (k <= slk(i)) then
             if (F_rig(i,k) > pbl_ricrit(1)) then
                F_turbreg(i,k) = LAMINAR
             else
                F_turbreg(i,k) = TURBULENT
             endif
          else
             F_turbreg(i,k) = TURBULENT
          endif
       enddo
    enddo
    
    ! Compute mixing and dissipation length scales
    mlen = ilongmel
    stat = ml_compute(F_zn, F_zd, F_pri, mlen, F_tt, F_hu, F_lwc, F_iwc, F_fbl, &
         F_gzt, F_gzm, F_gzm, F_sigt, F_sigm, F_sigm, F_ps, F_enold, F_buoy, &
         F_rig, zero3d, zero, F_turbreg, F_z0m, F_hpbl, F_lh, F_hpar, F_vcoef, &
         F_mrk2, F_dxdy, F_tau, F_kount)
    if (stat /= PHY_OK) then
       call physeterror('pbl_ri', 'error returned by mixing length calculation')
       return
    endif
    
    ! Change from gradient to flux form of buoyancy frequency and flux Richardson number
    do k=1,F_nkm1
       do i=1,F_ni          
          F_buoy(i,k) = F_pri(i,k) * F_buoy(i,k)
          F_rif(i,k) = F_pri(i,k) * F_rig(i,k)
       enddo
    enddo

    ! Solve the algebraic part of the TKE equation
    do k=1,F_nkm1
       do i=1,F_ni          
          b_term(i,k) = PBL_RI_CK * F_zn(i,k) * (F_shr2(i,k) - F_buoy(i,k))
          c_term(i,k) = PBL_RI_CE / F_zd(i,k)
       enddo
    enddo
    call tkestep(e_star, F_en, b_term, c_term, F_tau, F_ni, F_nkm1)
    if (phy_error_L) return

    ! Diagnose TKE budget terms on request
    if (ISREQSTEPL((/'BUEN','DSEN','SHEN'/))) then
       call tkebudget(F_en, e_star, b_term, c_term, F_zn, F_buoy, F_shr2, &
            F_buoyen, F_shren, F_dissen, F_tau, F_ni, F_nkm1)
       if (phy_error_L) return
    endif

    ! Prepare for diffusion term of the TKE equation
    call kcalc(ke, F_zn, F_enold, F_pri, F_vcoef, K_M, F_ni, F_nkm1)
    call ktosig(ke, ke, F_tt, F_hu, F_sigm, F_sigt, F_vcoef, K_M, F_ni, F_nkm1)
    do k=1,F_nkm1
       do i=1,F_ni
          ke(i,k) = pbl_ae * min(1. + max(F_rig(i,k), 0.), 3.) * ke(i,k)
       enddo
    enddo
    
    ! Compute surface boundary condition
    tke_sfc = PBL_RI_CU*F_frv**2 + PBL_RI_CW*F_wstar**2

    ! Diffuse TKE
    call diffuse(F_diffen, e_star, ke, tke_sfc, zero_sfc, &
         F_sigm, F_sigt, F_tau, FIELD_IS_TKE, F_ni, F_nkm1)
    if (phy_error_L) return

    ! Update TKE for the timestep
    if (F_kount > 0) F_en = max(etrmin2, e_star + F_tau * F_diffen)

    ! Compute eddy diffusivities and final buoyancy flux
    call kcalc(F_km, F_zn, F_en, F_pri, F_vcoef, K_M, F_ni, F_nkm1)
    call kcalc(F_kt, F_zn, F_en, F_pri, F_vcoef, K_T, F_ni, F_nkm1)
    F_buoy(:,:) = -F_kt(:,:) * F_buoy(:,:)

    ! Compute subgrid-scale moisture variance
    SGSVAR: if (pbl_progvar) then
       call thlqw_compute(thl, qw, F_tt, F_hu, F_lwc, F_iwc, F_sigt, F_ni, F_nkm1)
       if (F_kount == 0) then
          call covareq(F_qwvar, qw, qw, F_en, F_kt, F_zd, F_gzt, &
               F_vcoef, F_ni, F_nkm1)
       else
          call covarstep(qwvarstar, F_qwvar, qw, qw, &
               F_en, F_kt, F_zd, F_gzt, F_vcoef, F_tau, F_ni, F_nkm1)
          call ktosig(ktsg, F_kt, F_tt, F_hu, F_sigm, F_sigt, F_vcoef, K_T, F_ni, F_nkm1)
          qwvarstar(:,F_nkm1) = 4.3**2 * F_qstar(:)**2 !Eq. 63 of Mellado et al. (2017; QJRMS)
          call diffuse(qwvardiff, qwvarstar, ktsg, zero_sfc, zero_sfc, &
               F_sigm, F_sigt, F_tau, FIELD_ON_TLEV, F_ni, F_nkm1)
          F_qwvar(:,:) = qwvarstar(:,:) + F_tau * qwvardiff(:,:)
       endif
    endif SGSVAR
    
    ! End of subprogram
    return
  end subroutine rpnint

end module pbl_ri
