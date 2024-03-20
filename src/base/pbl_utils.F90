module pbl_utils
  implicit none
  private

  ! API functions
  public :: blheight                            !Diagnose height of the PBL
  public :: blweight                            !Set vertical tapering of PBL scheme influence
  public :: dissheat                            !Compute dissipative heating
  public :: dvrtdf                              !Compute vertical derivative
  public :: ficemxp                             !Diagnose ice fraction
  public :: sfcflux                             !Diagnose implicit surface fluxes

  ! API parameters
  integer, parameter, public :: LAMINAR = 1     !Code for laminar flow (hysteresis)
  integer, parameter, public :: TURBULENT = -1  !Code for turbulent flow (hysteresis)
  real, parameter, public :: WSTAR_MIN = 0.01   !Minimum convective velocity scale

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function blheight(hpbl,temp,hum,uwind,vwind,hghtm,hghtt,sigt,temp_s,hum_s,psfc, &
       z0m,z0t,lat,fcor,n,nk) result(status)
    use tdpack_const, only: CAPPA, CPD, EPS1, GRAV
    use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
    implicit none
#include <rmnlib_basics.hf>

    !@Arguments
    integer, intent(in) :: n                          !horizontal dimension
    integer, intent(in) :: nk                         !vertical dimension
    real, dimension(n), intent(in) :: temp_s          !surface temperature (K)
    real, dimension(n), intent(in) :: hum_s           !surface specific humidity (kg/kg)
    real, dimension(n), intent(in) :: psfc            !surface pressure (Pa)
    real, dimension(n), intent(in) :: z0m             !momentum mixing length (m)
    real, dimension(n), intent(in) :: z0t             !thermodynamic mixing length (m)
    real, dimension(n), intent(in) :: lat             !latitude (rad)
    real, dimension(n), intent(in) :: fcor            !Coriolis factor (/s)
    real, dimension(n,nk), intent(in) :: temp         !dry bulb temperature (K)
    real, dimension(n,nk), intent(in) :: hum          !specific humidity (kg/kg)
    real, dimension(n,nk), intent(in) :: uwind        !x-component wind (m/s)
    real, dimension(n,nk), intent(in) :: vwind        !y-component wind (m/s)
    real, dimension(n,nk), intent(in) :: hghtm        !momentum-level level heights (m)
    real, dimension(n,nk), intent(in) :: hghtt        !thermodynamic-level level heights (m)
    real, dimension(n,nk), intent(in) :: sigt         !sigma of thermodynamic levels
    real, dimension(n), intent(out) :: hpbl           !PBL depth (m)
    integer :: status                                 !completion status (RMN_OK or RMN_ERR)

    !@Author R. McTaggart-Cowan (Fall 2016)
    !@Object
    !          Calculate the depth of the PBL using the approach proposed
    !          by Seidel et al. (2012; JGR-Atm) based on the work of 
    !          Vogelezang and Holtslag (1996; BLM).  Variable names follow
    !          page 50 of IFS documentation (Chapter 3: Turbulent transport
    !          and interactions with the surface) for Cy40r1.
    !*@/

    ! Local parameter definitions
    real, parameter :: HGHT_NEARSFC=2.,RIB_CRIT=0.25,BETA=100.

    ! Local variable declarations
    integer :: i,k
    real :: svhbl,rib,rib_below,shrsq,hght,hght_below,zp
    real, dimension(n) :: theta_air,spd_air,dir_air,temp_n,hum_n,svn,uwind_n,vwind_n,ustar

    ! Initialize status
    status = RMN_ERR

    ! Compute basic near-surface thermodynamic properties
    if (sl_prelim(temp(:,nk),hum(:,nk),uwind(:,nk),vwind(:,nk),psfc,hghtm(:,nk), &
         spd_air=spd_air,dir_air=dir_air) /= SL_OK) then
       print*, 'Error in pbl_height for surface layer preparation'
       return
    endif
    theta_air = temp(:,nk)*sigt(:,nk)**(-CAPPA)
    if (sl_sfclayer(theta_air,hum(:,nk),spd_air,dir_air,hghtm(:,nk),hghtt(:,nk),temp_s,hum_s, &
         z0m,z0t,lat,fcor,hghtt_diag=HGHT_NEARSFC,hghtm_diag=HGHT_NEARSFC, &
         t_diag=temp_n,q_diag=hum_n,u_diag=uwind_n,v_diag=vwind_n,ue=ustar) /= SL_OK) then
       print*, 'Error in pbl_height for surface layer calculations'
       return
    endif

    ! Compute near-surface dry static energy
    svn = CPD*temp_n*(1+EPS1*hum_n) + GRAV*HGHT_NEARSFC

    ! Scan bulk Richardson numbers until the critical value is reached
    do i=1,n
       rib = 0.
       hght = HGHT_NEARSFC
       k = nk
       do while (rib < RIB_CRIT .and. k > 0)
          if (hghtt(i,k) <= HGHT_NEARSFC) then
             k = k-1
             cycle
          endif
          rib_below = rib
          hght_below = hght
          svhbl = CPD*temp(i,k)*(1+EPS1*hum(i,k)) + GRAV*hghtt(i,k)
          shrsq = (uwind(i,k)-uwind_n(i))**2+(vwind(i,k)-vwind_n(i))**2
          hght = 0.5*(hghtm(i,k)+hghtt(i,k))
          zp = (hghtm(i,k)-HGHT_NEARSFC)**2/(hghtt(i,k)-HGHT_NEARSFC)
          rib = zp*2*GRAV*(svhbl-svn(i)) / &
               ((svhbl+svn(i)-GRAV*hghtt(i,k)-GRAV*HGHT_NEARSFC)*(shrsq + BETA*ustar(i)**2))
          k = k-1
       enddo
       if (k == 0) then
          hpbl(i) = hght
       else
          hpbl(i) = hght - (rib-RIB_CRIT)/(rib-rib_below)*(hght-hght_below)
       endif
    enddo

    ! Successful completion
    status = RMN_OK
    return
  end function blheight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blweight(w,s,ps,n,nk)
    implicit none
!!!#include <arch_specific.hf>

    ! Arguments
    integer, intent(in) :: n                             !horizontal dimension
    integer, intent(in) :: nk                            !vertical dimension
    real, dimension(n,nk), intent(in) :: s               !sigma levels
    real, dimension(n), intent(in) :: ps                 !surface pressure (Pa)
    real, dimension(n,nk), intent(out) :: w              !weighting function [0,1]

    !Author
    !          J. Mailhot and B. Bilodeau (Dec 2001)

    !Revision
    !001       A-M. Leduc and S. Belair (Jun 2003) - ps as argument
    !                   blweight ---> blweight2. Change weighting
    !                   profile from sigma to pressure dependent.

    !Object
    !          Compute a weighting profile to be used with the moist
    !          turbulence scheme.

    !Notes
    !          The profile is set to:
    !            1 in the lower part of the atmosphere (S .ge. SMAX) if pres .ge. pmax
    !            0 in the upper part of the atmosphere (S .le. SMIN) if pres .le. pmin
    !            (with a linear interpolation in-between)

    ! Local parameters
    real, parameter :: PMIN=45000,PMAX=55000

    ! Local variables
    integer :: j,k
    real :: pres

    ! Compute profile of weighting functions
    do k=1,nk
       do j=1,n
          pres = s(j,k)*ps(j)
          w(j,k) = (1.-(pmax-pres)/(pmax-pmin))
          w(j,k) = min(max(w(j,k),0.),1.)
       end do
    end do

    return
  end subroutine blweight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine dissheat(F_ttendd, F_u, F_v, F_utend, F_vtend, F_km, F_sigm, F_sigt, &
       F_vcoef, F_tau, F_ni, F_nkm1)
    use phy_options, only: pbl_dissheat_opt => pbl_dissheat
    use vintphy, only: vint_mom2thermo1
    implicit none

    ! Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_u             !u-wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_v             !v-wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_utend         !u-wind tendency (m/s^2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_vtend         !v-wind tendency (m/s^2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_km            !diffusion coefficients for momentum
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigm          !sigma for momentum levels
    real, dimension(F_ni,F_nkm1), intent(in) :: F_sigt          !sigma for thermo levels
    real, dimension(*), intent(in) :: F_vcoef                   !coefficients for vertical interpolation
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_ttendd       !resultant temperature tendency (K/s)

    ! Local declarations
    integer :: i, k
    real :: dsig, du, dv
    real, dimension(F_ni,F_nkm1) :: dkem, tu2m, tu2t

    ! Dissipative heating
    select case (pbl_dissheat_opt)
    case ('LOCAL_TEND')
       dkem(:,:) = (F_u(:,1:F_nkm1) + 0.5*F_tau*F_utend(:,1:F_nkm1)) * F_utend(:,1:F_nkm1) + &
            (F_v(:,1:F_nkm1) + 0.5*F_tau*F_vtend(:,1:F_nkm1)) * F_vtend(:,1:F_nkm1)
       call vint_mom2thermo1(F_ttendd, dkem, F_vcoef, F_ni, F_nkm1)
       do i=1,F_ni
          dsig = (F_sigt(i,F_nkm1) - 1.) / (F_sigt(i,F_nkm1-1) - 1.)
          F_ttendd(i,F_nkm1) = dsig * F_ttendd(i,F_nkm1-1)
       end do
    case ('LOCAL_K')
       tu2m(:,:) = F_utend(:,1:F_nkm1)**2 + F_vtend(:,1:F_nkm1)**2
       call vint_mom2thermo1(tu2t, tu2m, F_vcoef, F_ni, F_nkm1)
       do k=1,(F_nkm1-1)
          do i=1,F_ni
             du = F_u(i,k) - F_u(i,k+1) + F_tau*(F_utend(i,k) - F_utend(i,k+1))
             dv = F_v(i,k) - F_v(i,k+1) + F_tau*(F_vtend(i,k) - F_vtend(i,k+1))
             dsig = F_sigm(i,k) - F_sigm(i,k+1)
             F_ttendd(i,k) = - F_km(i,k)*( (du/dsig)**2 + (dv/dsig)**2 ) &
                  - 0.5*F_tau*tu2t(i,k)
          end do
       end do
       do i=1,F_ni
          dsig = (F_sigt(i,F_nkm1) - 1.)/(F_sigt(i,F_nkm1-1) - 1.)
          F_ttendd(i,F_nkm1) = dsig * F_ttendd(i,F_nkm1-1)
       end do
    case DEFAULT
       F_ttendd(:,:) = 0.
    end select

  end subroutine dissheat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine dvrtdf(R , X , DS, N , MR , MX , NK)
    implicit none
!!!#include <arch_specific.hf>
    integer N, MR, MX, NK
    real R(MR,NK),X(MX,NK),DS(n,NK)

    !@Author R. Benoit RPN(Mar 1989)

    !@Revisions
    ! 001      R. Benoit (Aug 93) - DS(2D) for Local sigma

    !@Object calculate the vertical derivative by centred finite differences

    !@Arguments
    !          - Output -
    ! R        result
    !          - Input -
    ! X        variable to derive
    ! DS       distance between sigma levels 'U'
    ! N        horizontal dimensions
    ! MR       1st dimension of R
    ! MX       1st dimension of X
    ! NK       vertical dimension

    !@Notes
    !          R and X can share the same space, R(*,NK)=0

    integer J,K

    do K=1,NK-1
       do J=1,N
          R(J,K)=(X(J,K+1)-X(J,K))/DS(j,K)
       enddo
    enddo

    do J=1,N
       R(J,NK)=0
    enddo

    return
  end subroutine dvrtdf
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine ficemxp(fice,tf,df,t,n,nk)

    implicit none
!!!#include <arch_specific.hf>

    ! Arguments
    integer, intent(in) :: n                             !horizontal dimension
    integer, intent(in) :: nk                            !vertical dimension
    real, dimension(n,nk), intent(in) :: t               !dry air temperature (K)
    real, dimension(n,nk), intent(out) :: fice           !fraction of ice
    real, dimension(n,nk), intent(out) :: tf             !threshold for saturation wrt ice or liquid
    real, dimension(n,nk), intent(out) :: df             !derivative of fraction wrt T for ice or liquid

    !Author
    !          J. Mailhot (Jan 2000)

    !Object
    !          Calculate the fraction of ice and set the values of threshold
    !          and derivative w/r to T for computation of saturation values
    !          in the presence of mixed phases.

    !Notes
    !          Based on the definition in:
    !          - Burk et al. 1997, JGR 102, 16529-16544
    !          and the observations of:
    !          - Curry et al. 1990, Int. J. Climatol. 10, 749-764.
    !          - Curry et al. 1997, JGR 102, 13851-13860.
    !
    !          For F (fraction of ice), linear variation between Tmin and Tmax
    !          For TF and DF, values are set such that saturation is w/r to liquid for T > Tmin
    !                 "         "             "        saturation is w/r to ice    for T < Tmin

    ! Local parameter declarations
    real, parameter :: TMIN=248.16,TMAX=258.16

    ! Local variable declarations
    integer :: j,k
    real :: dt

    ! Compute ice fraction and saturation measures
    dt = 1.0/(tmax-tmin)
    do k=1,nk
       do j=1,n
          fice(j,k) = max(0.0,min(1.0,(tmax-t(j,k))*dt))
          tf(j,k) = fice(j,k)
          df(j,k) = -dt
          if( t(j,k).lt.tmin .or. t(j,k).gt.tmax) df(j,k) = 0.0
       end do
    end do

    return
  end subroutine ficemxp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine sfcflux(F_utau, F_vtau, F_fq, F_ustar, F_fsh, F_flw, &
       F_u, F_v, F_t, F_q, F_utend, F_vtend, F_ttend, F_qtend, &
       F_bmsg, F_btsg, F_alphaq, F_rhosfc, F_ps, F_tau, F_ni, F_nkm1)
    use tdpack, only: GRAV

    ! Argument declarations
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, intent(in) :: F_tau                                   !time step (s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_u             !u-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_v             !v-component wind (m/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_t             !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_q             !specific humidity (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_utend         !u-component wind tendency (m/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_vtend         !v-component wind tendency (m/s2)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_ttend         !dry air temperature tendency (K/s)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_qtend         !specific humidity tendency (kg/kg/s)
    real, dimension(F_ni), intent(in) :: F_bmsg                 !surface component of momentum flux
    real, dimension(F_ni), intent(in) :: F_btsg                 !surface component of heat flux
    real, dimension(F_ni), intent(in) :: F_alphaq               !inhomogeneous part of sfc heat flux
    real, dimension(F_ni), intent(in) :: F_rhosfc               !near-surface density (kg/m3) 
    real, dimension(F_ni), intent(in) :: F_ps                   !surface pressure (Pa)
    real, dimension(F_ni), intent(out) :: F_utau                !u-component surface stress
    real, dimension(F_ni), intent(out) :: F_vtau                !v-component surface stress
    real, dimension(F_ni), intent(out) :: F_fq                  !total surface stress
    real, dimension(F_ni), intent(out) :: F_ustar               !friction velocity (m/s)
    real, dimension(F_ni), intent(out) :: F_fsh                 !specific humidity flux (m/s)
    real, dimension(F_ni), intent(out) :: F_flw                 !water density flux (kg m-2 s-1)

    !@Object   Diagnose (implicit) surface fluxes
    
    ! Local variable declarations
    integer :: i
    real :: rhortvsg, mrhocmu, tplusnk, qplusnk, uplusnk, vplusnk, fsh

    ! Compute surface fluxes
    do i=1,F_ni

       ! Compute lowest-level time-plus state
       rhortvsg = F_ps(i) / GRAV
       mrhocmu  = rhortvsg * F_bmsg(i)
       tplusnk  = F_t(i,F_nkm1) + F_tau * F_ttend(i,F_nkm1)
       qplusnk  = F_q(i,F_nkm1) + F_tau * F_qtend(i,F_nkm1)
       uplusnk  = F_u(i,F_nkm1) + F_tau * F_utend(i,F_nkm1)
       vplusnk  = F_v(i,F_nkm1) + F_tau * F_vtend(i,F_nkm1)

       ! Diagnose surface momentum fluxes
       F_utau(i) = -mrhocmu * uplusnk
       F_vtau(i) = -mrhocmu * vplusnk
       F_fq(i) = -mrhocmu * sqrt(uplusnk**2 + vplusnk**2)
       F_ustar(i) = sqrt(F_fq(i) / F_rhosfc(i))

       ! Diagnose surface moisture fluxes
       fsh = F_alphaq(i) + F_btsg(i)*qplusnk
       F_flw(i) = rhortvsg * fsh
       F_fsh(i) = F_flw(i) / F_rhosfc(i)

    enddo

    ! End of subprogram
    return
  end subroutine sfcflux
  
end module pbl_utils
