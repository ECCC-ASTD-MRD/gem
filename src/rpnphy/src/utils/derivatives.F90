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

module derivatives
   use clib_itf_mod, only: clib_toupper
!!!#include <arch_specific.hf>
   implicit none
   private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Internal parameters
   integer, parameter :: LONG_CHAR=16                           !Long character string length
   real, parameter :: NOT_DEFINED=huge(1.)                      !Undefined value
   character(len=LONG_CHAR), parameter :: DEFAULT_TYPE='pchip'  !Default derivation method

   ! Define API
   integer, parameter, public :: DERIV_OK=0,DERIV_ERR=1         !Return statuses of functions
   public :: deriv_profile                                      !Compute derivatives for a profile

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function derivative_pchip(dydz,yi,zi,zo) result(status)
      ! Compute derivatives for a profile using PCHIP
      use pchip_utils, only: PCHIP_OK,PCHIP_EXT,pchip_extend,pchip_layers,pchip_coef
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Profile values
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic decreasing)
      real, dimension(:,:), intent(in) :: zo                    !Coordinate values for derivative (monotonic decreasing)
      real, dimension(:,:), intent(out) :: dydz                 !Derivative of profile
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (Mar 2018)

      !Object
      !          to calculate vertical derivatives of a given profile
      !          using PCHIP interpolation

      ! Local variable declaration
      integer :: n,nko,nkext,i,k,ko,istat
      real :: s
      real, dimension(size(yi,dim=1),size(yi,dim=2)+PCHIP_EXT) :: y,z,h,del,b,c,d

      ! Set error return status
      status = DERIV_ERR

      ! Initialize bounds
      n = size(yi,dim=1)
      nko = size(zo,dim=2)
      nkext = size(yi,dim=2)+PCHIP_EXT

      ! Invert input arrays if necessary
      if (pchip_extend(yi,zi,y,z) /= PCHIP_OK) then
         call msg(MSG_WARNING,'(integral_pchip) Error returned by pchip_extend')
         return
      endif

      ! Compute layer thickness and deltas
      istat = pchip_layers(y,z,h,del)

      ! Compute interpolating polynomial coefficients
      istat = pchip_coef(h,del,b,c,d)

      ! Compute derivatives
      do i=1,n
         k = 1
         do ko=1,nko
            do while (zo(i,ko) <= z(i,k+1) .and. k < nkext)
               k = k+1
            enddo
            s = z(i,k) - zo(i,ko)
            dydz(i,ko) = d(i,k) + 2.*s*c(i,k) + 3.*s**2*b(i,k)
         enddo
      enddo

      ! Successful completion
      status = DERIV_OK
      return
   end function derivative_pchip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function derivative_spline(dydz,yi,zi,zo,deriv_min,deriv_max) result(status)
      ! Compute derivatives for a profile using cublic splines
      use spline_utils, only: SPLINE_OK,spline_layers,spline_coef
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Profile values
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic decreasing)
      real, dimension(:,:), intent(in) :: zo                    !Coordinate values for derivative (monotonic decreasing)
      real, dimension(:), intent(in) :: deriv_min               !Boundary condition derivative at coordinate minimum
      real, dimension(:), intent(in) :: deriv_max               !Boundary condition  derivative at coordinate maximum
      real, dimension(:,:), intent(out) :: dydz                 !Derivative of profile
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (Mar 2018)

      !Object
      !          to calculate vertical derivatives of a given profile
      !          using PCHIP interpolation

      ! Local variable declaration
      integer :: n,nk,nko,i,k,ko,istat
      real :: s
      real, dimension(size(yi,dim=1),size(yi,dim=2)) :: h,del,b,c,d

      ! Set error return status
      status = DERIV_ERR

      ! Initialize bounds
      n = size(yi,dim=1)
      nk = size(yi,dim=2)
      nko = size(zo,dim=2)

      ! Compute layer thickness and deltas
      istat = spline_layers(yi,zi,h,del)

      ! Compute interpolating polynomial coefficients
      if (any(deriv_min == NOT_DEFINED) .and. any(deriv_max == NOT_DEFINED)) then
         istat = spline_coef(h,del,b,c,d)
      else if (any(deriv_min == NOT_DEFINED)) then
         istat = spline_coef(h,del,b,c,d,slope_max=deriv_max)
      else if (any(deriv_max == NOT_DEFINED)) then
         istat = spline_coef(h,del,b,c,d,slope_min=deriv_min)
      else
         istat = spline_coef(h,del,b,c,d,slope_min=deriv_min,slope_max=deriv_max)
      endif

      ! Compute derivatives
      do i=1,n
         k = 1
         do ko=1,nko
            do while (zo(i,ko) <= zi(i,k+1) .and. k < nk-1)
               k = k+1
            enddo
            if (zo(i,ko) <= zi(i,k+1)) k = k+1
            s = zi(i,k) - zo(i,ko)
            dydz(i,ko) = d(i,k) + 2.*s*c(i,k) + 3.*s**2*b(i,k)
         enddo
      enddo

      ! Successful completion
      status = DERIV_OK
      return
   end function derivative_spline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function deriv_profile(dydz,yi,zi,zo,type,deriv_min,deriv_max) result(status)
      ! Derivative of a quantity over a profile

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Profile values
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic decreasing)
      real, dimension(:,:), intent(in) :: zo                    !Coordinate values for derivative (monotonic decreasing)
      character(len=*), intent(in), optional :: type            !Type of derivation strategy to use (['pchip'],spline)
      real, dimension(:), intent(in), optional :: deriv_min     !Boundary condition derivative at coordinate minimum
      real, dimension(:), intent(in), optional :: deriv_max     !Boundary condition  derivative at coordinate maximum
      real, dimension(:,:), intent(out) :: dydz                 !Derivative of profile
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (Mar 2018)

      !Object
      !          Wrap call to derivative calculations.

      ! Local variables
      integer :: istat
      real, dimension(size(zi,dim=1)) :: myDeriv_min,myDeriv_max
      character(len=LONG_CHAR) :: myType

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type
      myDeriv_min = NOT_DEFINED
      if (present(deriv_min)) myDeriv_min = deriv_min
      myDeriv_max = NOT_DEFINED
      if (present(deriv_max)) myDeriv_max = deriv_max

      ! Derivative calculation dispatcher
      istat = clib_toupper(myType)
      select case (myType)
      case ('PCHIP')
         status = derivative_pchip(dydz,yi,zi,zo)
      case ('SPLINE')
         status = derivative_spline(dydz,yi,zi,zo,myDeriv_min,myDeriv_max)
      case DEFAULT
         call msg(MSG_WARNING,'(deriv_profile) invalid type of derivative ('//trim(myType)//') requested')
      end select

      ! End of subprogram
      return
   end function deriv_profile

end module derivatives
