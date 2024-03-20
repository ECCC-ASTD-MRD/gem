!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------

module pchip_utils
!!!#include <arch_specific.hf>
   implicit none
   private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Define API
   integer, parameter, public :: PCHIP_OK=0,PCHIP_ERR=1         !Return statuses of functions
   integer, parameter, public :: PCHIP_EXT=2                    !Number of points added to profile for boundaries
   integer, parameter, public :: PCHIP_DIR_UP=3                 !Direction Up
   integer, parameter, public :: PCHIP_DIR_DOWN=4               !Direction Down
   public :: pchip_extend                                       !Extend input column and invert if necessary
   public :: pchip_layers                                       !Compute layer thicknesses and deltas
   public :: pchip_coef                                         !Compute coefficients for interpolating polynomial

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function pchip_extend(yi,zi,y,z,ni,nki,nko,z1i,z2i,z1,z2,direc) result(status)
      ! Extend input column and invert if necessary
      integer, intent(in) :: ni,nki,nko                         !Arrays Dimensions
      real, dimension(ni,nki), intent(in) :: yi                 !Value of input profile
      real, dimension(ni,nki), intent(in) :: zi                 !Vertical coordinate (monotonic decreasing)
      real, dimension(ni,nko), intent(out) :: y                 !Value of extended input profile (with vertical "halos")
      real, dimension(ni,nko), intent(out) :: z                 !Vertical coordinate (with vertical "halos")
      real, dimension(ni), intent(in), optional :: z1i,z2i      !Specified input bounds for integration
      real, dimension(ni), intent(out), optional :: z1,z2       !Integration bounds appropriate for extended profile
      integer, intent(in), optional :: direc                    !Direction of profile for integration (flips profile)
      integer :: status                                         !Return status

      ! Local variables
      integer :: k,i,ko
      real, dimension(ni) :: myZ1i,myZ2i
      real, dimension(ni,nko) :: yext,zext
      integer :: myDirec

      ! Set return status
      status = PCHIP_ERR

      ! Set default values
      myDirec = PCHIP_DIR_UP
      if (present(direc)) myDirec = direc
      do i=1,ni
         myZ1i(i) = zi(i,1)
         if (present(z1i)) myZ1i(i) = z1i(i)
         myZ2i(i) = zi(i,nki)
         if (present(z2i)) myZ2i(i) = z2i(i)
      enddo
      
#ifdef DEBUG
      ! Check direction of array
      if (any(zi(:,1)<zi(:,nki))) then
         call msg(MSG_WARNING,'(pchip_extend) Inverted coordinate detected')
         y = 0.; z = 0.
         if (present(z1)) z1 = 0.
         if (present(z2)) z2 = 0.
         return
      endif

      ! Check size of array
      if (nko /= nki + PCHIP_EXT) then
         call msg(MSG_WARNING,'(pchip_extend) Output profile size should be nk_input + PCHIP_EXT')
         y = 0.; z = 0.
         if (present(z1)) z1 = 0.
         if (present(z2)) z2 = 0.
         return
      endif
#endif

      ! Extend arrays to integration bounds
      do k=1,nki
         ko = k + PCHIP_EXT/2
         do i=1,ni
            y(i,ko) = yi(i,k)
            z(i,ko) = zi(i,k)
         enddo
      enddo
      do i=1,ni
         y(i,1) = yi(i,1)
         y(i,nko) = yi(i,nki)
         z(i,1) = max(myZ2i(i), 2*zi(i,1)-zi(i,2))
         z(i,nko) = min(myZ1i(i), 2*zi(i,nki)-zi(i,nki-1))
      enddo

      ! Invert arrays if necessary
      select case (myDirec)
      case (PCHIP_DIR_DOWN)
         zext(:,:) = z(:,:)
         yext(:,:) = y(:,:)
         do k=1,nko
            do i=1,ni
               z(i,nko-k+1) = zi(i,1) - zext(i,k)
               y(i,nko-k+1) = yext(i,k)
            enddo
         enddo
         if (present(z1)) z1 = zi(:,1) - myZ2i
         if (present(z2)) z2 = zi(:,1) - myZ1i
      case (PCHIP_DIR_UP)
         if (present(z1)) z1(:) = myZ1i(:)
         if (present(z2)) z2(:) = myZ2i(:)
      case DEFAULT
         call msg(MSG_WARNING,'(pchip_extend) Invalid direc provided (use PCHIP_DIR_UP or PCHIP_DIR_DOWN)')
         return
      end select

      ! Successful completion
      status = PCHIP_OK
   end function pchip_extend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function pchip_layers(y,z,h,del,ni,nk) result(status)
      ! Compute layer thicknesses and deltas
      integer, intent(in) :: ni,nk
      real, dimension(ni,nk), intent(in) :: y      !Value of extended input profile (with vertical "halos")
      real, dimension(ni,nk), intent(in) :: z      !Vertical coordinate (with vertical "halos")
      real, dimension(ni,nk), intent(out) :: h     !Thickness of layers in profile
      real, dimension(ni,nk), intent(out) :: del   !Linear first derivatives of profile
      integer :: status                            !Return status

      ! Local variables
      integer :: k,i

      ! Set return status
      status = PCHIP_ERR
      
      ! Compute layer properties
      do k=2,nk
         do i=1,ni
            h(i,k) = z(i,k-1)-z(i,k)
            del(i,k) = (y(i,k-1)-y(i,k))/h(i,k)
         enddo
      enddo
      do i=1,ni
         h(i,1)   = h(i,2)
         del(i,1) = del(i,2)
      enddo

      ! Successful completion
      status = PCHIP_OK
      return
   end function pchip_layers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function pchip_coef(h,del,b,c,d,ni,nk) result(status)
      ! Compute coefficients for the interpolating polynomial
      integer, intent(in) :: ni,nk                              !Arrays Dimensions
      real, dimension(ni,nk), intent(in) :: h                   !Thickness of layers in profile
      real, dimension(ni,nk), intent(in) :: del                 !Linear first derivatives of profile
      real, dimension(ni,nk), intent(out) :: b                  !Coefficient for cubic term
      real, dimension(ni,nk), intent(out) :: c                  !Coefficient for square term
      real, dimension(ni,nk), intent(out) :: d                  !Coefficient for linear term
      integer :: status                                         !Return status

      ! The interpolating polynomial has the form :
      !           P(x) = y + d*s + c*s^2 + b*s^3
      ! where s = x-x(k) is the distance between x(k) and x(k+1) for the interpolation.

      ! Local variables
      integer :: k,i
      real :: w1,w2,w3,delk,delkp,delkpk,invh(ni,nk)

      real, parameter :: eps4 = tiny(1.)

      ! Set return status
      status = PCHIP_ERR
      
      do i=1,ni
         d(i,nk) = ((2.*h(i,nk) + h(i,nk-1))*del(i,nk) - h(i,nk)*del(i,nk-1)) / (h(i,nk) + h(i,nk-1))
         if (sign(1.,d(i,nk)) /= sign(1.,del(i,nk))) d(i,nk) = 0.
         if (sign(1.,del(i,nk)) /= sign(1.,del(i,nk-1))) then
            if (abs(d(i,nk)) < abs(3*del(i,nk))) d(i,nk) = 3.*del(i,nk)
         endif
      enddo

      invh = 1./h
      do k=nk-1,2,-1
         do i=1,ni
            delk = del(i,k)
            delkp = del(i,k+1)
            delkpk = delk*delkp
            w1 = 2.*h(i,k) + h(i,k+1)
            w2 = 2.*h(i,k+1) + h(i,k)
            w3 = delk*w1+delkp*w2
            d(i,k) = (w1+w2)*max(0., delkpk) / (sign(1.,w3)*max(eps4,abs(w3)))
            c(i,k+1) = (3.*del(i,k+1) - 2.*d(i,k+1) - d(i,k)) * invh(i,k+1)
            b(i,k+1) = (d(i,k+1) - 2.*del(i,k+1) + d(i,k)) * (invh(i,k+1)*invh(i,k+1))
         enddo
      enddo

      do i=1,ni
         d(i,1) = ((2.*h(i,2) + h(i,3))*del(i,2) - h(i,2)*del(i,3)) / (h(i,2) + h(i,3))
         if (sign(1.,d(i,2)) /= sign(1.,del(i,2))) d(i,1) = 0.
         if (sign(1.,del(i,2)) /= sign(1.,del(i,3))) then
            if (abs(d(i,1)) > abs(3.*del(i,2))) d(i,1) = 3.*del(i,1)
         endif
         c(i,1) = 0.
         b(i,1) = 0.
         c(i,2) = (3.*del(i,2) - 2.*d(i,2) - d(i,1)) * invh(i,2)
         b(i,2) = (d(i,2) - 2.*del(i,2) + d(i,1) ) * (invh(i,2)*invh(i,2))
      enddo

      ! Successful completion
      status = PCHIP_OK
      return
   end function pchip_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function pchip_slopes(h,del,d,ni,nk) result(status)
      ! Compute profile slopes
      integer, intent(in) :: ni,nk
      real, dimension(ni,nk), intent(in) :: h,del
      real, dimension(ni,nk), intent(out) :: d
      integer :: status                                         !Return status

      ! Local variables
      integer :: i,k
      real :: w1,w2,delk,delkp

      ! Set return status
      status = PCHIP_ERR

      ! Compute slopes at interior points
      do k=2,nk-1
         do i=1,ni
            delk = del(i,k)
            delkp = del(i,k+1)
            if (delk == 0. .or. delkp == 0. .or. sign(1.,delk)*sign(1.,delkp) < 0.) then
               d(i,k) = 0.
            else
               w1 = 2.*h(i,k) + h(i,k+1)
               w2 = 2.*h(i,k+1) + h(i,k)
               d(i,k) = (w1+w2)*delkp*delk / (delk*w1+delkp*w2)
            endif
         enddo
      enddo

      ! Compute slopes at end points
      do i=1,ni
         d(i,1) = ((2.*h(i,2) + h(i,3))*del(i,2) - h(i,2)*del(i,3)) / (h(i,2) + h(i,3))
         if (sign(1.,d(i,2)) /= sign(1.,del(i,2))) d(i,1) = 0.
         if (sign(1.,del(i,2)) /= sign(1.,del(i,3))) then
            if (abs(d(i,1)) > abs(3*del(i,2))) d(i,1) = 3.*del(i,1)
         endif
      enddo
      do i=1,ni
         d(i,nk) = ((2.*h(i,nk) + h(i,nk-1))*del(i,nk) - h(i,nk)*del(i,nk-1)) / (h(i,nk) + h(i,nk-1))
         if (sign(1.,d(i,nk)) /= sign(1.,del(i,nk))) d(i,nk) = 0.
         if (sign(1.,del(i,nk)) /= sign(1.,del(i,nk-1))) then
            if (abs(d(i,nk)) < abs(3*del(i,nk))) d(i,nk) = 3.*del(i,nk)
         endif
      enddo
 
      ! Successful completion
      status = PCHIP_OK
      return
   end function pchip_slopes

end module pchip_utils
