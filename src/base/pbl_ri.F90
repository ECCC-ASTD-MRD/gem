module pbl_ri
  implicit none
  private
  
  ! Public API functions
  public :: rpnint

contains

  subroutine rpnint(F_en, F_km, F_kt, F_pri, F_rif, F_rig, F_buoy, F_shr2, &
       F_diffen, F_dissen, F_zn, F_zd, F_enold, F_u, F_v, F_t, F_tve, F_qv, &
       F_qce, F_fbl,  F_hpbl, F_lh, F_ps, F_z0, F_frv, F_wstar, F_hpar, &
       F_sigt, F_sige, F_ze, F_dxdy, F_mrk2, F_vcoef, F_tau, F_kount, F_ni, F_nkm1)
    
    use, intrinsic :: iso_fortran_env, only: INT64
    use tdpack, only: GRAV, KARMAN, RGASD
    use phy_options, only: ilongmel, etrmin2, phyoutlist_s, nphyoutlist
    use phy_status, only: phy_error_L, PHY_OK
    use vintphy, only: vint_thermo2mom2
    use mixing_length, only: ml_compute, ML_LMDA
    use ens_perturb, only: ens_nc2d
    use pbl_ri_utils, only: buoyflux, tkestep, tkebudget, PBL_RI_CK, PBL_RI_CU, PBL_RI_CW, PBL_RI_CE, &
         PBL_RI_CAB, PBL_RI_AE
    use pbl_ri_diffuse, only: diffuse, FIELD_ON_ELEV
    implicit none
#include <rmnlib_basics.hf>
    
    !Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    integer, intent(in) :: F_kount                              !time step number
    real, intent(in) :: F_tau                                   !time step length (s)
    real, dimension(F_ni), intent(in) :: F_hpbl                 !height of the PBL (m)
    real, dimension(F_ni), intent(in) :: F_lh                   !launching height (m)
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(F_ni), intent(in) :: F_z0                   !roughness length (m)
    real, dimension(F_ni), intent(in) :: F_frv                  !friction velocity (m/s)
    real, dimension(F_ni), intent(in) :: F_wstar                !convective velocity scale (m/s)
    real, dimension(F_ni), intent(in) :: F_dxdy                 !horizontal grid area (m^2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_enold         !TKE of previous time step (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qv            !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_u             !u-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_v             !v-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_t             !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tve           !virt. temperature on e-lev (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sige          !sigma for energy levels
    real, dimension(F_ni,ens_nc2d), intent(in) :: F_mrk2        !Markov chains for stochastic parameters
    real, dimension(*), intent(in) :: F_vcoef                   !coefficients for vertical interpolation
    real, dimension(F_ni,F_nkm1), intent(in) :: F_ze            !height of energy levs (m)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qce           !PBL cloud water content on e-lev (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_fbl           !PBL cloud fraction
    real, dimension(F_ni), intent(inout) :: F_hpar              !height of parcel ascent (m)
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_zn         !mixing length (m)
    real, dimension(F_ni,F_nkm1), intent(inout) :: F_zd         !dissipation length (m)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_en           !updated TKE (m2/s2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_km           !eddy diffusivity for momentum (m2/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_kt           !eddy diffusivity for heat/moisture (m2/s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_pri          !inverse Prandtl number (ratio of KT/KM diffusion coefficients)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_rif          !flux Richardson number
    real, dimension(F_ni,F_nkm1), intent(out) :: F_rig          !gradient Richardson number
    real, dimension(F_ni,F_nkm1), intent(out) :: F_buoy         !buoyancy flux (s^-2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_shr2         !square of vertical wind shear (s^-2)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_diffen       !diffusion of TKE
    real, dimension(F_ni,F_nkm1), intent(out) :: F_dissen       !dissipation of TKE
    
    !@Object
    !          Update TKE and diffusion coefficients for the RPN Integrated PBL scheme
    include "phyinput.inc"

    ! Local variable declarations
    integer :: stat, i, k
    integer, dimension(F_ni) :: mlen
    real, dimension(F_ni) :: tke_sfc, beta_sfc
    real, dimension(F_ni,F_nkm1) :: dudz, dvdz, e_star, asig, ke, shr_term, &
         buoy_term, qe, te, b_term, c_term, zero, one
    real, dimension(F_ni,F_nkm1,3) :: zero3d

    ! Initialization
    if (F_kount == 0) then
       if (.not.any('zn' == phyinread_list_s(1:phyinread_n))) then
          do k=1,F_nkm1
             F_zn(:,k) = min(KARMAN*(F_ze(:,k) + F_z0(:)), ML_LMDA)
          enddo
       endif
    endif
    
    ! Interpolate input profiles to correct set of levels
    call vint_thermo2mom2(te,  qe, &
         &                F_t, F_qv, F_vcoef, F_ni, F_nkm1)
    
    ! Estimate buoyancy and shear gradients
    call buoyflux(F_u, F_v, F_t, F_tve, F_qv, F_qce, F_fbl, dudz, dvdz, &
         F_sigt, F_ps, F_shr2, F_rig, F_buoy, F_vcoef, F_ni, F_nkm1)

    ! Set inverse Prandtl number
    F_pri = 1.

    ! Compute mixing and dissipation length scales
    zero = 0.
    zero3d = 0.
    one = 1.
    mlen = ilongmel
    stat = ml_compute(F_zn, F_zd, F_pri, mlen, te, qe, F_qce, F_fbl, F_ze, F_ze, &
         F_sige, F_sige, F_ps, F_enold, F_buoy, F_rig, zero3d, zero, one, zero, &
         F_z0, F_hpbl, F_lh, F_hpar, F_mrk2, F_dxdy, F_tau, F_kount)
    if (stat /= PHY_OK) then
       call physeterror('pbl_ri', 'error returned by mixing length calculation')
       return
    endif

    ! Change from gradient to flux form of buoyancy flux and flux Richardson number
    do k=1,F_nkm1
       do i=1,F_ni
          F_buoy(i,k) = F_pri(i,k) * F_buoy(i,k)
          F_rif(i,k) = F_pri(i,k) * F_rig(i,k)
       enddo
    enddo

    ! Update turbulent kinetic energy
    UPDATE_TKE: if (F_kount > 0) then

       ! Solve the algabraic part of the TKE equation
       do k=1,F_nkm1
          do i=1,F_ni
             b_term(i,k) = PBL_RI_CK * F_zn(i,k) * (F_shr2(i,k) - F_buoy(i,k))
             c_term(i,k) = PBL_RI_CE / F_zd(i,k)
          enddo
       enddo
       call tkestep(e_star, F_en, b_term, c_term, F_tau, F_ni, F_nkm1)
       if (phy_error_L) return
       
       ! Diagnose TKE budget terms on request
       if (any([(any((/'dsen'/) == phyoutlist_S(i)), i=1,nphyoutlist)])) then
          call tkebudget(F_en, e_star, b_term, c_term, F_zn, F_shr2, F_buoy, &
               shr_term, buoy_term, F_dissen, F_tau, F_ni, F_nkm1)
          if (phy_error_L) return
       endif

       ! Prepare for diffusion term of the TKE equation
       asig = (GRAV/RGASD) * F_sige/F_tve
       do k=1,F_nkm1
          ke(:,k) = PBL_RI_CK * F_zn(:,k) * sqrt(F_enold(:,k)) * asig(:,k)**2
       enddo

       ! Compute surface boundary condition
       tke_sfc = PBL_RI_CU*F_frv**2 + PBL_RI_CW*F_wstar**2
       beta_sfc = 0.

       ! Diffuse TKE
       call diffuse(F_diffen, e_star, PBL_RI_AE*ke, tke_sfc, beta_sfc, &
            F_sige, F_sigt, F_tau, FIELD_ON_ELEV, F_ni, F_nkm1)
       if (phy_error_L) return       
       
       ! update TKE for the timestep
       F_en = max(etrmin2, e_star + F_tau * F_diffen)
       
    endif UPDATE_TKE

    ! Compute eddy diffusivities
    do k=1,F_nkm1
       do i=1,F_ni
          F_km(i,k) = PBL_RI_CK * F_zn(i,k) * sqrt(F_en(i,k))
          F_kt(i,k) = F_km(i,k) * F_pri(i,k)
          F_buoy(i,k) = -F_kt(i,k) * F_buoy(i,k)
       enddo
    enddo

    ! End of subprogram
    return
  end subroutine rpnint

end module pbl_ri
