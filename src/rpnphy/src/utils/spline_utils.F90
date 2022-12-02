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

module spline_utils
!!!#include <arch_specific.hf>
   implicit none
   private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Internal parameters
   integer, parameter :: LONG_CHAR=16                           !Long character string length

   ! Define API
   integer, parameter, public :: SPLINE_OK=0,SPLINE_ERR=1       !Return statuses of functions
   integer, parameter, public :: SPLINE_EXT=2                   !Number of points added to profile for boundaries
   public :: spline_layers                                      !Compute layer thicknesses and deltas
   public :: spline_coef                                        !Compute coefficients for interpolating polynomial

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function spline_layers(y,z,h,del) result(status)
      ! Compute layer thicknesses and deltas
      implicit none
      real, dimension(:,:), intent(in) :: y                     !Value of extended input profile (with vertical "halos")
      real, dimension(:,:), intent(in) :: z                     !Vertical coordinate (with vertical "halos")
      real, dimension(:,:), intent(out) :: h                    !Thickness of layers in profile
      real, dimension(:,:), intent(out) :: del                  !Linear first derivatives of profile
      integer :: status                                         !Return status

      ! Local variables
      integer :: k,nk

      ! Set return status
      status = SPLINE_ERR

      ! Compute layer properties
      nk = size(y,dim=2)
      do k=1,nk-1
         h(:,k)   =  z(:,k+1)-z(:,k) 
         del(:,k) = (y(:,k+1)-y(:,k))/h(:,k) 
      enddo
      h(:,nk)   = h(:,nk-1)
      del(:,nk) = del(:,nk-1)

      ! Successful completion
      status = SPLINE_OK
      return
   end function spline_layers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function spline_coef(h,del,b,c,d,slope_min,slope_max) result(status)
      ! Compute coefficients for the interpolating polynomial
      implicit none
      real, dimension(:,:), intent(in) :: h                     !Thickness of layers in profile
      real, dimension(:,:), intent(in) :: del                   !Linear first derivatives of profile       
      real, dimension(:), intent(in), optional :: slope_min     !Derivative at coordinate minimum (use clamped BC)
      real, dimension(:), intent(in), optional :: slope_max     !Derivative at coordinate maximum (use clamped BC)
      real, dimension(:,:), intent(out) :: b                    !Coefficient for cubic term
      real, dimension(:,:), intent(out) :: c                    !Coefficient for square term
      real, dimension(:,:), intent(out) :: d                    !Coefficient for linear term
      integer :: status                                         !Return status

      ! The interpolating polynomial has the form :
      !           P(x) = y + d*s + c*s^2 + b*s^3
      ! where s = x-x(k) is the distance between x(k) and x(k+1) for the interpolation.

      ! Local variables
      integer :: k,n,nk
      real, dimension(size(h,dim=1),size(h,dim=2)) :: ta,tb,tc,tr,work

      ! Set return status
      status = SPLINE_ERR

      ! Set dimensions
      n = size(h,dim=1)
      nk = size(h,dim=2)

      ! Prepare tridiagonal matrix
      ta(:,1) = 0.
      tb(:,1) = h(:,2)
      tc(:,1) = h(:,1)+h(:,2)
      do k=2,nk-1
         ta(:,k) = h(:,k)
         tb(:,k) = 2*(h(:,k-1)+h(:,k))
         tc(:,k) = h(:,k-1)
      enddo
      ta(:,nk) = h(:,nk-1)+h(:,nk-2)
      tb(:,nk) = h(:,nk-2)
      tc(:,nk) = 0.

      ! Prepare right-hand side
      tr(:,1) = ((h(:,1)+2*tc(:,1))*h(:,2)*del(:,1)+h(:,1)**2*del(:,2))/tc(:,1)
      do k=2,nk-1
         tr(:,k) = 3*(h(:,k)*del(:,k-1) + h(:,k-1)*del(:,k))
      enddo
      tr(:,nk) = (h(:,nk-1)**2*del(:,nk-2)+(2*ta(:,nk)+h(:,nk-1))*h(:,nk-2)*del(:,nk-1))/ta(:,nk)

      ! Make clamping adjustments if boundary conditions are supplied
      if (present(slope_max)) then
         tb(:,1) = 1.
         tc(:,1) = 0.
         tr(:,1) = slope_max(:)
      endif
      if (present(slope_min)) then
         ta(:,nk) = 0.
         tb(:,nk) = 1.
         tr(:,nk) = slope_min(:)
      endif

      ! Solve the matrix problem for the slopes (d)
      call difuvd2(d,ta,tb,tc,tr,work,n,n,nk)

      ! Compute coefficients for higher order terms
      nk = size(h,dim=2)
      do k=1,nk-1
         c(:,k) = (3.*del(:,k) - 2.*d(:,k) - d(:,k+1)) / h(:,k)
         b(:,k) = (dble(d(:,k)) - dble(2.*del(:,k)) + dble(d(:,k+1)) ) / dble(h(:,k)**2)
      enddo
      c(:,nk) = 0.; b(:,nk) = 0.

      ! Successful completion
      status = SPLINE_OK
      return
   end function spline_coef

end module spline_utils
