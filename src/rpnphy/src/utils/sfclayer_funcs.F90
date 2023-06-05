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

module sfclayer_funcs
  implicit none
  private

  public :: stability_function

#include <rmnlib_basics.hf>

  ! Stability function interface prototype
  abstract interface
     subroutine stability_function(F_fm, F_fh, F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, &
          F_zm, F_zh, F_dfm, F_dfh)
       real, intent(out) :: F_fm                        !integrated stability function for momentum
       real, intent(out) :: F_fh                        !integrated stability function for heat
       real, intent(in) :: F_lzz0                       !neutral (log) stability function for momentum
       real, intent(in) :: F_lzz0t                      !neutral (log) stability function for heat
       real, intent(in) :: F_ilmo                       !inverse of Obukhov length (1/m)
       real, intent(in) :: F_h                          !height of the PBL (m)
       real, intent(in) :: F_beta                       !Prandtl number for a neutral profile
       real, intent(in), dimension(2) :: F_zm           !heights for momentum function calculation as (/high,low/)
       real, intent(in), dimension(2) :: F_zh           !heights for heat function calculation as (/high,low/)
       real, intent(out), optional :: F_dfm             !derivative of momentum stability function wrt ilmo
       real, intent(out), optional :: F_dfh             !derivative of heat stability function wrt ilmo
     end subroutine stability_function
  end interface

  ! Utility to acquire stability functions
  public :: sf_get

  ! Stability functions
  public :: sf_stable_delage97, &
       sf_stable_lock07, &
       sf_stable_beljaars91, &
       sf_unstable_delage92, &
       sf_unstable_dyer74

  ! Public parameters
  integer, parameter, public :: SF_OK    = RMN_OK       !return code for success
  integer, parameter, public :: SF_ERROR = RMN_ERR      !return code for error

  ! Parameters for selecting class of stability function
  integer, parameter, public :: SF_MOMENTUM = 0         !momentum stability function id
  integer, parameter, public :: SF_HEAT = 1             !heat stability function id

  ! Parameters for the Beljaars and Holtslag (1991) stability functions (stable)
  real, parameter, public :: BH91_A = 1.                !A coefficient
  real, parameter, public :: BH91_B = 2./3.             !B coefficient
  real, parameter, public :: BH91_C = 5.                !C coefficient
  real, parameter, public :: BH91_D = 0.35              !D coefficient

  ! Parameters for the Lock (2007) stability functions (stable)
  real, parameter, public :: L07_AH = 1.                !AH coefficient [was 4 in Lock (2007)]
  real, parameter, public :: L07_AM = 1.                !AM coefficient [was 4 in Lock (2007)]

  ! Parameters for the Delage (1997) stability functions (stable)
  real, parameter, public :: D97_AS = 12.               !AS coefficient

  ! Parameters for the Delage and Girard (1991) stability functions (unstable)
  real, parameter, public :: DG92_CI = 40.              !CI coefficient
  real, parameter, public :: DG92_RAC3 = 1.732050776    !RAC3 coefficient

  ! Package-local variables

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sf_get(func, stability, val) result(status)
    ! Obtain the stability function pointer as requested
    procedure(stability_function), pointer :: func      !retrieved stability function
    character(len=*), intent(in) :: stability           !stability branch to get (stable or unstable)
    character(len=*), intent(in) :: val                 !value of the function
    integer :: status                                   !return status of the operation
    status = SF_ERROR
    if (stability == 'stable') then
       select case (val)
       case ('DELAGE97')
          func => sf_stable_delage97
       case ('LOCK07')
          func => sf_stable_lock07
       case ('BELJAARS91')
          func => sf_stable_beljaars91
       case DEFAULT
          nullify(func)
          return
       end select
    else if (stability == 'unstable') then
       select case (val)
       case ('DELAGE92')
          func => sf_unstable_delage92
       case ('DYER74')
          func => sf_unstable_dyer74
       case DEFAULT
          nullify(func)
          return
       end select
    else
       return
    endif
    status = SF_OK
  end function sf_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sf_stable_delage97(F_fm, F_fh, F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, &
       F_zm, F_zh, F_dfm, F_dfh)
    ! Stability functions from Delage (1997; BLM)
    real, intent(in) :: F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, F_zm(2), F_zh(2)
    real, intent(out) :: F_fm, F_fh
    real, intent(out), optional :: F_dfm, F_dfh
    integer :: j
    real :: hi,a,b,c,d
    real, dimension(2) :: sf, dsf

    hi = 1./F_h
    d  = 4*D97_AS*F_beta*F_ilmo
    c  = d*hi - hi**2
    b  = d - 2*hi

    ! Stability function for momentum
    do j=1,2
       a  = SQRT(1 + b*F_zm(j) - c*F_zm(j)**2)
       sf(j) = 0.5*( a - F_zm(j)*hi - LOG(1+b*F_zm(j)*0.5+a) &
            - b/(2*SQRT(c))*ASIN((b-2*c*F_zm(j))/d) )
       if (present(F_dfm)) then
          dsf(j) = 0.5*( 1/(2*a)*((d*F_zm(j)-b*hi/c)*(1-F_zm(j)*hi) &
               - d*F_zm(j)*(1-F_zm(j)*hi+a)/(1+0.5*b*F_zm(j)+a)) &
               -(d**2*hi/(4*c**1.5))*ASIN((b-2*c*F_zm(j))/d) )
       endif
    enddo
    F_fm = sfadjust(F_lzz0, sf, 1.)
    if (present(F_dfm)) F_dfm = sfderiv(dsf, 1.)

    ! Stability function for heat
    do j=1,2
       a  = SQRT(1 + b*F_zh(j) - c*F_zh(j)**2)
       sf(j) = 0.5*( a - F_zh(j)*hi - LOG(1+b*F_zh(j)*0.5+a) &
            - b/(2*SQRT(c))*ASIN((b-2*c*F_zh(j))/d) )
       if (present(F_dfh)) then
          dsf(j) = 0.5*( 1/(2*a)*((d*F_zh(j)-b*hi/c)*(1-F_zh(j)*hi) &
               - d*F_zh(j)*(1-F_zh(j)*hi+a)/(1+0.5*b*F_zh(j)+a)) &
               -(d**2*hi/(4*c**1.5))*ASIN((b-2*c*F_zh(j))/d) )
       endif
    enddo
    F_fh = sfadjust(F_lzz0t, sf, F_beta)
    if (present(F_dfh)) F_dfh = sfderiv(dsf, F_beta)
  end subroutine sf_stable_delage97

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sf_stable_lock07(F_fm, F_fh, F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, &
       F_zm, F_zh, F_dfm, F_dfh)
    ! Stability functions from Lock (2007; UKMO Tech Report)
    real, intent(in) :: F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, F_zm(2), F_zh(2)
    real, intent(out) :: F_fm, F_fh
    real, intent(out), optional :: F_dfm, F_dfh
    integer :: j
    real :: a
    real, dimension(2) :: sf, dsf

    ! Stability function for momentum
    do j=1,2
       a = L07_AM
       sf(j) = a*F_zm(j)*F_ilmo
       if (present(F_dfm)) dsf(j) = a*F_zm(j)*F_ilmo
    enddo
    F_fm = sfadjust(F_lzz0, sf, 1.)
    if (present(F_dfm)) F_dfm = sfderiv(dsf, 1.)

    ! Stability function for heat
    do j=1,2
       a = L07_AH
       sf(j) = (1.+ 0.5*a*F_zh(j)*F_ilmo)**2
       if (present(F_dfh)) dsf(j) = a*F_zh(j)*F_ilmo*(1.+ 0.5*a*F_zh(j)*F_ilmo)
    enddo
    F_fh = sfadjust(F_lzz0t, sf, F_beta)
    if (present(F_dfh)) F_dfh = sfderiv(dsf, F_beta)
  end subroutine sf_stable_lock07
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sf_stable_beljaars91(F_fm, F_fh, F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, &
       F_zm, F_zh, F_dfm, F_dfh)
    ! Stability functions from Beljaars and Holtslag (1991; JAM)
    real, intent(in) :: F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, F_zm(2), F_zh(2)
    real, intent(out) :: F_fm, F_fh
    real, intent(out), optional :: F_dfm, F_dfh
    integer :: j
    real :: a,b,c,d
    real, dimension(2) :: sf, dsf
 
    a = BH91_A
    b = BH91_B
    c = BH91_C
    d = BH91_D

    ! Stability function for momentum
    do j=1,2
       sf(j) = a*F_zm(j)*F_ilmo &
            + b*(F_zm(j)*F_ilmo-c/d)*EXP(-d*F_zm(j)*F_ilmo) + b*c/d
       if (present(F_dfm)) dsf = a*F_zm(j)*F_ilmo &
            + b*F_zm(j)*F_ilmo*EXP(-d*F_zm(j)*F_ilmo) &
            - d*F_zm(j)*F_ilmo*b*(F_zm(j)*F_ilmo-c/d)*EXP(-d*F_zm(j)*F_ilmo)
    enddo
    F_fm = sfadjust(F_lzz0, sf, 1.)
    if (present(F_dfm)) F_dfm = sfderiv(dsf, 1.)
 
    ! Stability function for heat
    do j=1,2
       sf(j) = (1.+b*a*F_zh(j)*F_ilmo)**1.5 - 1 &
            + b*(F_zh(j)*F_ilmo-c/d)*EXP(-d*F_zh(j)*F_ilmo) + b*c/d
       if (present(F_dfh)) dsf = a*F_zh(j)*F_ilmo*(1.+b*a*F_zh(j)*F_ilmo)**0.5 &
            + b*F_zh(j)*F_ilmo*EXP(-d*F_zh(j)*F_ilmo) &
            - d*F_zh(j)*F_ilmo*b*(F_zh(j)*F_ilmo-c/d)*EXP(-d*F_zh(j)*F_ilmo)
    enddo
    F_fh = sfadjust(F_lzz0t, sf, F_beta)
    if (present(F_dfh)) F_dfh = sfderiv(dsf, F_beta)
  end subroutine sf_stable_beljaars91

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sf_unstable_delage92(F_fm, F_fh, F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, &
       F_zm, F_zh, F_dfm, F_dfh)
    ! Stability functions from Delage and Girard (1992; BLM)
    real, intent(in) :: F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, F_zm(2), F_zh(2)
    real, intent(out) :: F_fm, F_fh
    real, intent(out), optional :: F_dfm, F_dfh
    integer :: j
    real :: a
    real, dimension(2) :: sf, dsf

    ! Stability function for momentum
    do j=1,2
       a =(1-DG92_CI*F_zm(j)*F_beta*F_ilmo)**(0.16666666)
       sf(j) = - LOG( (a+1)**2*SQRT(a**2-a+1)*(a**2+a+1)**1.5 ) &
            +DG92_RAC3*ATAN((a**2-1)/(DG92_RAC3*a))
       if (present(F_dfm)) dsf(j) = (1./a) - 1
    enddo
    F_fm = sfadjust(F_lzz0, sf, 1.)
    if (present(F_dfm)) F_dfm = sfderiv(dsf, 1.)
 
    ! Stability function for heat
    do j=1,2
       a = (1-DG92_CI*F_zh(j)*F_beta*F_ilmo)**(0.33333333)
       sf(j) = -1.5*LOG(a**2+a+1)+DG92_RAC3*ATAN((2*a+1)/DG92_RAC3)
       if (present(F_dfh)) dsf(j) = (1./a) - 1
    enddo
    F_fh = sfadjust(F_lzz0t, sf, F_beta)
    if (present(F_dfh)) F_dfh = sfderiv(dsf, F_beta)
  end subroutine sf_unstable_delage92

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sf_unstable_dyer74(F_fm, F_fh, F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, &
       F_zm, F_zh, F_dfm, F_dfh)
    ! Stability functions from Dyer (1974; BLM)
    real, intent(in) :: F_lzz0, F_lzz0t, F_ilmo, F_h, F_beta, F_zm(2), F_zh(2)
    real, intent(out) :: F_fm, F_fh
    real, intent(out), optional :: F_dfm, F_dfh
    integer :: j
    real :: a
    real, dimension(2) :: sf, dsf

    ! Stability function for momentum
    do j=1,2
       a = (1-16*F_zm(j)*F_ilmo)**(0.25)
       sf(j) = - 2.*LOG(a+1) - LOG(a**2+1) + 2.*ATAN(a)
       if (present(F_dfm)) dsf = (1./a) - 1
    enddo
    F_fm = sfadjust(F_lzz0, sf, 1.)
    if (present(F_dfm)) F_dfm = sfderiv(dsf, 1.)

    ! Stability function for heat
    do j=1,2
       a = (1-16*F_zh(j)*F_ilmo)**(0.5)
       sf(j) = - 2.*LOG(a+1)
       if (present(F_dfh)) dsf = (1./a) - 1
    enddo
    F_fh = sfadjust(F_lzz0t, sf, F_beta)
    if (present(F_dfh)) F_dfh = sfderiv(dsf, F_beta)
  end subroutine sf_unstable_dyer74

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sfadjust(lzz0, sf, beta) result(sfadj)
    ! Adjust neutral stability function
    real, intent(in) :: lzz0              !neutral (log) stability function
    real, dimension(2), intent(in) :: sf  !integrated function at boundaries
    real, intent(in) :: beta              !scale factor to apply
    real :: sfadj                         !adjusted stability function
    sfadj = beta * (lzz0 + sf(1) - sf(2))
    return
  end function sfadjust

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sfderiv(dsf, beta) result(sfder)
    ! Compute stability function derivative (primitive)
    real, dimension(2), intent(in) :: dsf !function derivative at boundaries
    real, intent(in) :: beta              !scale factor to apply
    real :: sfder                         !stability function derivative
    sfder = beta * (dsf(1) - dsf(2))
    return
  end function sfderiv

end module sfclayer_funcs
