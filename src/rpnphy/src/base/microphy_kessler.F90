!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module microphy_kessler
  use, intrinsic :: iso_fortran_env
  implicit none
  private

  ! Public procedures
  public :: kessler             !Condensation scheme
  public :: kessler_phybusinit  !Define bus requirements
  public :: kessler_lwc         !Compute liquid water content
  public :: kessler_iwc         !Compute ice water content

contains

  subroutine kessler(F_tttnd, F_qvtnd, F_qctnd, F_qrtnd, F_rainflux, &
       F_tt, F_qv, F_qc, F_qr, F_gz, F_sigma, F_ps, F_dt)
    use tdpack, only: TTNS3W, CHLC, TCDK, RGASD, DELTA, EPS1, CPD, foewa
    implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

    ! Arguments
    real, dimension(:,:), intent(out) :: F_tttnd          !Temperature tendency (K/s)
    real, dimension(:,:), intent(out) :: F_qvtnd          !Specific humidity tendency (kg/kg/s)
    real, dimension(:,:), intent(out) :: F_qctnd          !Cloud water tendency (kg/kg/s)
    real, dimension(:,:), intent(out) :: F_qrtnd          !Rain water tendency (kg/kg/s)
    real, dimension(:), intent(out) :: F_rainflux         !Surface rain flux (kg/m2/s)
    real, dimension(:,:), intent(in) :: F_tt              !Dry air temperature (K)
    real, dimension(:,:), intent(in) :: F_qv              !Specific humidity (kg/kg)
    real, dimension(:,:), intent(in) :: F_qc              !Specific cloud water mass (kg/kg)
    real, dimension(:,:), intent(in) :: F_qr              !Specific rain water mass (kg/kg)
    real, dimension(:,:), intent(in) :: F_gz              !Height of thermo levels (m)
    real, dimension(:,:), intent(in) :: F_sigma           !Eta coordinate of thermo levels
    real, dimension(:), intent(in) :: F_ps                !Surface pressure (Pa)
    real, intent(in) :: F_dt                              !Time step (s)

    !@Author   R. McTaggart-Cowan (2022)

    !@Revision

    !@Object
    !          Implementation of a Kessler-type simplified condensation
    !          scheme as used in the DCMIP-2016 project and described by
    !            Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
    !            Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
    !            Radius Sphere. Journal of Advances in Modeling Earth Systems.
    !            doi:10.1002/2015MS000435

    !@Notes
    !          Original code by Paul Ullrich, Kevin Reed and Joe Klemp
    !*@/

    ! Internal parameters
    real :: F5 = 237.3 * TTNS3W * CHLC / CPD

    ! Internal variables
    integer :: i, k, ni, nk, rainsplit, nt
    real :: precl, ern, qrprod, prod, qvs, dt_max, dt0, es
    real, dimension(size(F_tt,dim=2)) :: qv, qc, qr, &
         rho, velqr, sed, tt
    
    ! Set dimensions
    ni = size(F_tt, dim=1)
    nk = size(F_tt, dim=2)

    ! Process all columns in the slab
    COLUMNS: do i=1,ni

       ! Compute inputs and invert column (model top at nk)
       do k=nk,1,-1
          qv(k) = F_qv(i,k) / (1. - F_qv(i,k))
          qc(k) = F_qc(i,k) / (1. - F_qv(i,k))
          qr(k) = F_qr(i,k) / (1. - F_qv(i,k))
          rho(k) = F_sigma(i,k) * F_ps(i) / &
               (RGASD * F_tt(i,k) * (1. + DELTA * F_qv(i,k)) * (1. + qv(k)))
          tt(k) = F_tt(i,k)
       enddo

       
       ! Maximum time step size in accordance with CFL condition
       velqr(:) = fallspeed(qr(:), rho(:))
       dt_max = F_dt
       do k=2,nk
          if (velqr(k) /= 0.) then
             dt_max = min(dt_max, 0.8 * (F_gz(i,k-1) - F_gz(i,k)) / velqr(k))
          end if
       end do
       
       ! Number of subcycles
       rainsplit = ceiling(F_dt / dt_max)
       dt0 = F_dt / real(rainsplit, 8)
       
       ! Subcycle through rain process
       precl = 0.
       SUBSTEP: do nt=1,rainsplit
          
          ! Surface precipitation flux kg/m2/s)
          precl = precl + rho(nk) * qr(nk) * velqr(nk)
          
          ! Sedimentation term using upstream differencing
          do k=2,nk
             sed(k) = dt0 * (rho(k-1) * qr(k-1) * velqr(k-1) - rho(k) * qr(k) * velqr(k)) / &
                  (rho(k) * (F_gz(i,k-1) - F_gz(i,k)))
          enddo
          sed(1)  = dt0 * qr(1) * velqr(1) / (.5 * (F_gz(i,1) - F_gz(i,2)))
          
          ! Adjustment terms
          do k=1,nk
             
             ! Autoconversion and accretion rates following KW eq. 2.13a,b
             qrprod = qc(k) - (qc(k) - dt0 * max(.001 * (qc(k) - .001), 0.0)) / &
                  (1. + dt0 * 2.2 * qr(k)**.875)
             qc(k) = max(qc(k) - qrprod, 0.)
             qr(k) = max(qr(k) + qrprod + sed(k), 0.)
             
             ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
             es = foewa(tt(k))
             qvs = EPS1 * es / (F_sigma(i,k) * F_ps(i) - es)
             prod = (qv(k) - qvs) / (1. + qvs * F5 / (tt(k) - 36.)**2)
             
             ! Evaporation rate following KW eq. 2.14a,b
             ern = min(dt0 * (((1.6 + 124.9 * (0.001*rho(k) * qr(k))**.2046) * &
                  (0.001*rho(k) * qr(k))**.525) / &
                  ( 2.55e8 / (F_sigma(i,k) * F_ps(i) * qvs) + 5.4e5 )) * &
                  (dim(qvs, qv(k)) / (0.001*rho(k) * qvs)), &
                  max(-prod-qc(k), 0.), qr(k))
             
             ! Saturation adjustment following KW eq. 3.10
             tt(k) = tt(k) + CHLC/CPD * (max(prod, -qc(k)) - ern)
             qv(k) = max(qv(k) - max(prod, -qc(k)) + ern, 0.)
             qc(k) = qc(k) + max(prod, -qc(k))
             qr(k) = qr(k) - ern
          enddo
          
          ! Recalculate liquid water terminal velocity
          if (nt /= rainsplit) velqr(:) = fallspeed(qr(:), rho(:))
       enddo SUBSTEP

       ! Average precipitation rate
       precl = precl / real(rainsplit)

       ! Compute and invert tendencies and derived fields
       F_rainflux(i) = precl
       do k=1,nk
          F_tttnd(i,k) = (tt(k) - F_tt(i,k)) / F_dt
          F_qvtnd(i,k) = (qv(k)/(1. + qv(k)) - F_qv(i,k)) / F_dt
          F_qctnd(i,k) = (qc(k)/(1. + qv(k)) - F_qc(i,k)) / F_dt
          F_qrtnd(i,k) = (qr(k)/(1. + qv(k)) - F_qr(i,k)) / F_dt
       enddo
       
    enddo COLUMNS

  end subroutine kessler

  ! Define bus requirements
  function kessler_phybusinit() result(F_istat)
    use phy_status, only: PHY_OK, PHY_ERROR
    use bus_builder, only: bb_request
    implicit none
    integer :: F_istat                          !Function return status
    F_istat = PHY_ERROR
    if (bb_request((/ &
         'CLOUD_WATER_MASS', &
         'RAIN_MASS       ' &
         /)) /= PHY_OK) then
       call physeterror('microphy_kessler::kessler_phybusinit', &
            'Cannot construct bus request list')
       return
    endif
    F_istat = PHY_OK
  end function kessler_phybusinit

  ! Compute total water mass
  function kessler_lwc(F_qltot, F_dbus, F_pbus, F_vbus) result(F_istat)
    use phybus
    use phy_status, only: PHY_OK, PHY_ERROR
    implicit none
    real, dimension(:,:), intent(out) :: F_qltot        !Total water mass (kg/kg)
    real, dimension(:), pointer, contiguous :: F_dbus   !Dynamics bus
    real, dimension(:), pointer, contiguous :: F_pbus   !Permanent bus
    real, dimension(:), pointer, contiguous :: F_vbus   !Volatile bus
    integer :: F_istat                                  !Return status
#include "phymkptr.hf"
    integer :: ni, nkm1
    real, dimension(:,:), pointer :: zqcp
    F_istat = PHY_ERROR
    ni = size(F_qltot, dim=1); nkm1 = size(F_qltot, dim=2)
    MKPTR2Dm1(zqcp, qcplus, F_dbus)
    F_qltot(:,:) = zqcp(:,:)
    F_istat = PHY_OK
    return
  end function kessler_lwc

  ! Compute total ice mass
  function kessler_iwc(F_qitot, F_dbus, F_pbus, F_vbus) result(F_istat)
    use phybus
    use phy_status, only: PHY_OK, PHY_ERROR
    implicit none
    real, dimension(:,:), intent(out) :: F_qitot        !Total ice mass (kg/kg)
    real, dimension(:), pointer, contiguous :: F_dbus   !Dynamics bus
    real, dimension(:), pointer, contiguous :: F_pbus   !Permanent bus
    real, dimension(:), pointer, contiguous :: F_vbus   !Volatile bus
    integer :: F_istat                                  !Return status
    F_istat = PHY_ERROR
    F_qitot(:,:) = 0.
    F_istat = PHY_OK
    return
  end function kessler_iwc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Internal utilities
  function fallspeed(qr, rho) result(velqr)
    real, dimension(:), intent(in) :: qr, rho
    real, dimension(size(qr)) :: velqr
    integer :: k,nk
    nk = size(qr)
    do k=1,nk
       velqr(k) = 36.34*(qr(k)*(0.001*rho(k)))**0.1364*sqrt(rho(nk)/rho(k))
    enddo
  end function fallspeed
  
end MODULE microphy_kessler
