!-------------------------------------- LICENCE BEGIN --------------------------
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
!-------------------------------------- LICENCE END ----------------------------

module integrals
   use, intrinsic :: iso_fortran_env, only: REAL64, REAL32
   use pchip_utils, only: PCHIP_OK,PCHIP_EXT,PCHIP_DIR_UP,PCHIP_DIR_DOWN, &
        pchip_coef
   use rootfind, only: RF_OK,rf_nrbnd
!!!#include <arch_specific.hf>
   implicit none
   private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Generic functions
   interface int_profile
      module procedure integral_profile_vec
      module procedure integral_profile_vec_dble
      module procedure integral_profile_sum
      module procedure integral_profile_sum_dble
      module procedure integral_profile_sum_scalbnd
      module procedure integral_profile_sum_scalbnd_dble
   end interface int_profile
   interface int_solve
      module procedure integral_solve_vdep
      module procedure integral_solve_sdep
      module procedure integral_solve_sdepsa
   end interface int_solve
 
   ! Define API
   integer, parameter, public :: INT_TYPE_PCHIP  = 1
   integer, parameter, public :: INT_TYPE_STEP   = 2
   integer, parameter, public :: INT_TYPE_LINEAR = 3
   integer, parameter, public :: INT_DIR_UP      = PCHIP_DIR_UP
   integer, parameter, public :: INT_DIR_DOWN    = PCHIP_DIR_DOWN
   integer, parameter, public :: INT_OK=0,INT_ERR=1             !Return statuses of functions
   public :: int_profile                                        !Compute integrals for a profile
   public :: int_solve                                          !Solve an integral equation

   ! Internal parameters
   integer, parameter :: DEFAULT_TYPE=INT_TYPE_PCHIP            !Default integration method
   real, parameter :: TRD = 1./3.

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_pchip(zi,z1i,z2i,n,nk, &
        yi,Aid,Ais,Sid,Sis) result(status)
      implicit none

      !@Arguments
      integer, intent(in) :: n, nk
      real, intent(in) :: zi(n,nk)        !Vertical coordinate (monotonic, increasing from nk to 1)
      real, intent(in) :: z1i(n)          !Lower limit of integration
      real, intent(in) :: z2i(n)          !Upper limit of integration
      real, intent(in) :: yi(n,nk)        !Integrand
      
      real(REAL64), intent(out), optional :: Aid(n,nk)   !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                         !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      real, intent(out), optional :: Ais(n,nk)           !Integrated profile (as Aid but real32)
      real(REAL64), intent(out), optional :: Sid(n)      !Integrated sum (Aid(:,1))
      real, intent(out), optional :: Sis(n)              !Integrated sum (as Sid but real32)

      integer :: status                                  !Return status of function

      !@Author A.Zadra (Aug 2016)
      !@Object calculate vertical integrals of a given profile
      !        using PCHIP interpolation
      !        Integrate a quantity over a profile using PCHIP,
      !        returning integrated quantity along the profile

      ! Local variable declaration
      real(REAL64), dimension(n) :: delta1, delta2
      real :: w0,w1,w2
      integer, dimension(n) :: k1,k2
      integer :: nko,i,k,istat,ko,k1max,k2min
      real(REAL64) :: hk,yk,dk,ck,bk,r,s,hk2,hk3,hk4,r2,r3,r4,s2,s3,s4
      real(REAL64) :: delta
      real, dimension(n) :: z1,z2
      real, dimension(n,nk+PCHIP_EXT) :: z,y,h,del,b,c,d,h2,h4
      real(REAL64), dimension(n,nk+PCHIP_EXT) :: A   !# REAL64 takes 130% more time
      logical :: bnd_invert

      ! Set error return status
      status = INT_ERR

      ! Determine direction of integration
      ! Extend input arrrays
      ! Invert input arrays if necessary
      nko = nk+2
      bnd_invert = any(z1i > z2i)
      do i=1,n
         if (bnd_invert) then
            z1(i) = z2i(i)
            z2(i) = z1i(i)
         else
            z1(i) = z1i(i)
            z2(i) = z2i(i)
         endif
         y(i,1) = yi(i,1)
         y(i,nko) = yi(i,nk)
         z(i,1) = max(z2(i), 2*zi(i,1)-zi(i,2))
         z(i,nko) = min(z1(i), 2*zi(i,nk)-zi(i,nk-1))
      enddo
      do k=1,nk
         ko = k + 1
         do i=1,n
            y(i,ko) = yi(i,k)
            z(i,ko) = zi(i,k)
            h(i,ko) = z(i,ko-1)-z(i,ko)
            del(i,ko) = (y(i,ko-1)-y(i,ko))/h(i,ko)
         enddo
      enddo
      do i=1,n
         h(i,1) = h(i,2) 
         del(i,1) = del(i,2)
         h(i,nko) = z(i,nko-1)-z(i,nko)
         del(i,nko) = (y(i,nko-1)-y(i,nko))/h(i,nko)
      enddo

      ! Compute interpolating polynomial coefficients
      istat = pchip_coef(h,del,b,c,d,n,nko)
      
      ! Find subdomain
      istat = subdomain(z,z1,z2,k1,k2,n,nko)
      
      ! Integrate profile - upper and lower bounds values 
      do i=1,n
         k = min(max(1, k1(i)), nko)
         yk = y(i,k)
         dk = d(i,k)*0.5
         ck = c(i,k)*TRD
         bk = b(i,k)*0.25
         r = z1(i) - z(i,k)
         r2 = r*r
         r3 = r2*r
         r4 = r3*r
         delta1(i) =  -1. * (&
              r*yk + &
              r2*dk + &
              r3*ck + &
              r4*bk)
      enddo
      do i=1,n
         k = min(max(1, k2(i)), nko)
         hk  = h(i,k)
         hk2 = hk*hk
         hk3 = hk2*hk
         hk4 = hk3*hk
         yk = y(i,k)
         dk = d(i,k)*0.5
         ck = c(i,k)*TRD
         bk = b(i,k)*0.25
         s = z2(i) - z(i,k)
         s2 = s*s
         s3 = s2*s
         s4 = s3*s
         delta2(i) =  &
              (s - hk)*yk + &
              (s2 - hk2)*dk + &
              (s3 - hk3)*ck + &
              (s4 - hk4)*bk
      enddo

      ! Integral profile by interpolation
      A(:,nko) = 0.
      do k=nko,3,-1
         do i=1,n
            w1 = float(1 - min(abs(k-k1(i)), 1))
            w2 = float(1 - min(abs(k-k2(i)), 1))
            w0 = float(min( &
                 min(max(k1(i)-k, 0), 1), &
                 min(max(k-k2(i), 0), 1) &
                 ))
            w0 = max(w0,w1,w2)
            hk  = h(i,k)
            hk2 = hk*hk
            hk3 = hk2*hk
            hk4 = hk3*hk
            yk = y(i,k)
            dk = d(i,k)*0.5
            ck = c(i,k)*TRD
            bk = b(i,k)*0.25
            delta =  &
                 hk*yk + &
                 hk2*dk + &
                 hk3*ck + &
                 hk4*bk
            A(i,k-1) = A(i,k) + w0*(delta &
                 + w1 * delta1(i) &
                 + w2 * delta2(i) &
                 )
         enddo
      enddo
      
      ! Return interior of domain only
      ! Negate result if bounds were inverted
      if (bnd_invert) then
         if (present(Aid)) Aid = -A(:,2:nk+1)
         if (present(Ais)) Ais = -A(:,2:nk+1)
         if (present(Sid)) Sid = -A(:,2)
         if (present(Sis)) Sis = -A(:,2)
      else
         if (present(Aid)) Aid = A(:,2:nk+1)
         if (present(Ais)) Ais = A(:,2:nk+1)
         if (present(Sid)) Sid = A(:,2)
         if (present(Sis)) Sis = A(:,2)
      endif
   
      ! Successful completion
      status = INT_OK
      return
   end function integral_pchip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_linear(zi,z1i,z2i,type,n,nk, &
        yi,Aid,Ais,Sid,Sis) result(status)
      ! Integrate a quantity over a profile using linear interpolation, returning integrated quantity along the profile
      implicit none

      ! Arguments
      integer, intent(in) :: n, nk
      real, intent(in) :: zi(n,nk)        !Vertical coordinate (monotonic, increasing from nk to 1)
      real, intent(in) :: z1i(n)          !Lower limit of integration
      real, intent(in) :: z2i(n)          !Upper limit of integration
      integer, intent(in) :: type         !Linear subtype (INT_TYPE_LINEAR or INT_TYPE_STEP)
      real, intent(in) :: yi(n,nk)        !Integrand
      
      real(REAL64), intent(out), optional :: Aid(n,nk)   !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                         !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      real, intent(out), optional :: Ais(n,nk)           !Integrated profile (as Aid but real32)
      real(REAL64), intent(out), optional :: Sid(n)      !Integrated sum (Aid(:,1))
      real, intent(out), optional :: Sis(n)              !Integrated sum (as Sid but real32)

      integer :: status                                  !Return status of function

      !Author
      !          R. McTaggart-Cowan and A.Zadra (July 2017)

      !Object
      !          to calculate vertical integrals of a given profile
      !          using linear interpolation

      ! Local variable declaration
      integer, dimension(n) :: k1, k2
      real(REAL64), dimension(n) :: delta1, delta2
      integer :: istat,k1i,k2i,kp
      real :: invzz1,invzz2,w0,w1,w2
      real :: zmid1,zmid2
      integer :: nko,i,k,ko
      real :: zmid_step_mask
      real, dimension(n) :: z1,z2
      real, dimension(n,nk+2) :: z,y,zh
      real(REAL64), dimension(n,nk+2) :: A  !# REAL64 takes 160% more time
      logical :: bnd_invert

      ! Set error return status
      status = INT_ERR

      ! Determine direction of integration
      ! Extend input arrrays
      ! Define layer interfaces
      nko = nk+2
      bnd_invert = any(z1i > z2i)
      do i=1,n
         if (bnd_invert) then
            z1(i) = z2i(i)
            z2(i) = z1i(i)
         else
            z1(i) = z1i(i)
            z2(i) = z2i(i)
         endif
         y(i,1) = yi(i,1)
         y(i,nko) = yi(i,nk)
         z(i,1) = max(z2(i), 2*zi(i,1)-zi(i,2))
         z(i,nko) = min(z1(i), 2*zi(i,nk)-zi(i,nk-1))
         zh(i,1) = 0.
         zh(i,nko) = 0.
      enddo
      do k=1,nk
         ko = k + 1
         do i=1,n
            y(i,ko) = yi(i,k)
            z(i,ko) = zi(i,k)
            zh(i,ko) = 0.5*(z(i,ko) + z(i,ko-1))
         enddo
      enddo
      
      ! Find subdomain
      istat = subdomain2(zh,z1,z2,k1,k2,n,nko)

      ! Integrate profile - upper and lower bounds values 
      zmid_step_mask = 1.
      if (type == INT_TYPE_STEP) zmid_step_mask = 0.
      do i=1,n
         k1i = min(max(2, k1(i)), nko)
         zmid1 = 0.5*(zh(i,k1i)+z1(i))
         zmid1 = zmid_step_mask * zmid1 + (1.-zmid_step_mask) * z(i,k1i)
         invzz1 = 1./(z(i,k1i-1)-z(i,k1i))
         delta1(i) = (y(i,k1i) + (y(i,k1i-1)-y(i,k1i))*((zmid1-z(i,k1i))*invzz1)) * (zh(i,k1i)-z1(i))
      enddo
      do i=1,n
         k2i = min(max(1, k2(i)), nko-1)
         zmid2 = 0.5*(zh(i,k2i+1)+z2(i))
         zmid2 = zmid_step_mask * zmid2 + (1.-zmid_step_mask) * z(i,k2i)
         invzz2 = 1./(z(i,k2i)-z(i,k2i+1))
         delta2(i) = (y(i,k2i) - (y(i,k2i)-y(i,k2i+1))*((z(i,k2i)-zmid2)*invzz2)) * (z2(i)-zh(i,k2i+1))
      enddo
      
      ! Integrate profile
      A(:,nko) = 0.
      do k=nko,2,-1
         kp = min(k+1, nko)
         do i=1,n
            w1 = float(1 - min(abs(k-k1(i)), 1))
            w2 = float(1 - min(abs(k-k2(i)), 1))
            w0 = float(min( &
                 min(max(k1(i)-k, 0), 1), &
                 min(max(k-k2(i), 0), 1) &
                 ))
            w0 = max(w0,w1,w2)
            A(i,k) = A(i,kp) + w0*( &
                 w1 * delta1(i) + &
                 w2 * delta2(i) + &
                 (1. - max(w1,w2)) * y(i,k) * (zh(i,k)-zh(i,kp)) &
                 )
         enddo
      enddo
      
      ! Return interior of domain only
      ! Negate result if bounds were inverted
      if (bnd_invert) then
         if (present(Aid)) Aid = -A(:,2:nk+1)
         if (present(Ais)) Ais = -A(:,2:nk+1)
         if (present(Sid)) Sid = -A(:,2)
         if (present(Sis)) Sis = -A(:,2)
      else
         if (present(Aid)) Aid = A(:,2:nk+1)
         if (present(Ais)) Ais = A(:,2:nk+1)
         if (present(Sid)) Sid = A(:,2)
         if (present(Sis)) Sis = A(:,2)
      endif

      ! Successful completion
      status = INT_OK
      return
   end function integral_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_profile_vec_dble(Ai,yi,zi,z1i,z2i,type) result(status)
      ! Integrate a quantity over a profile, returning integrated quantity along the profile

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, dimension(:), target, intent(in) :: z1i             !Lower limit of integration
      real, dimension(:), target, intent(in) :: z2i             !Upper limit of integration
      integer, intent(in), optional :: type                     !Type of integration strategy to use (INT_TYPE_PCHIP, INT_TYPE_LINEAR or INT_TYPE_STEP)
      real(REAL64), dimension(:,:), intent(out) :: Ai           !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                                !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (July 2017)

      !Object
      !          Dispatch requests for vertical integration.

      ! Local variables
      integer :: myType, ni, nk
 
      ! Initialize return value
      status = INT_ERR

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ni = size(yi,dim=1)
      nk = size(yi,dim=2)
      ! Dispatch to integration function
      select case (myType)
      case (INT_TYPE_PCHIP)
         status = integral_pchip(zi,z1i,z2i,ni,nk,yi,aid=ai)
      case (INT_TYPE_LINEAR, INT_TYPE_STEP)
         status = integral_linear(zi,z1i,z2i,myType,ni,nk,yi,aid=ai)
      case DEFAULT
         call msg(MSG_WARNING,'(integral_profile_vec_dble) invalid type of integral requested')
         return
      end select

      ! Error handling
      if (status /= INT_OK) then
         call msg(MSG_WARNING,'(integral_profile_vec_dble) error returned from integration subprogram')
         return
      end if

      ! End of subprogram
      return
   end function integral_profile_vec_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_profile_vec(Ai,yi,zi,z1i,z2i,type) result(status)
      ! Integrate a quantity over a profile, returning integrated quantity along the profile

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, dimension(:), target, intent(in) :: z1i             !Lower limit of integration
      real, dimension(:), target, intent(in) :: z2i             !Upper limit of integration
      integer, intent(in), optional :: type                     !Type of integration strategy to use (INT_TYPE_PCHIP, INT_TYPE_LINEAR or INT_TYPE_STEP)
      real, dimension(:,:), intent(out) :: Ai                   !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                                !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (May 2017)

      !Object
      !          Wrap call to double-precision calculation.

      ! Local variables
      integer :: myType, ni, nk

      ! Initialize return value
      status = INT_ERR

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ni = size(yi,dim=1)
      nk = size(yi,dim=2)
      ! Dispatch to integration function
      select case (myType)
      case (INT_TYPE_PCHIP)
         status = integral_pchip(zi,z1i,z2i,ni,nk,yi,ais=ai)
      case (INT_TYPE_LINEAR, INT_TYPE_STEP)
         status = integral_linear(zi,z1i,z2i,myType,ni,nk,yi,ais=ai)
      case DEFAULT
         call msg(MSG_WARNING,'(integral_profile_vec) invalid type of integral requested')
         return
      end select

      ! Error handling
      if (status /= INT_OK) then
         call msg(MSG_WARNING,'(integral_profile_vec) error returned from integration subprogram')
         return
      end if
      
      ! End of subprogram
      return
   end function integral_profile_vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_profile_sum_dble(Si,yi,zi,z1i,z2i,type) result(status)
      ! Integrate a quantity over a profile, returning integrated value

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, dimension(:), intent(in) :: z1i                     !Lower limit of integration
      real, dimension(:), intent(in) :: z2i                     !Upper limit of integration
      integer, intent(in), optional :: type                     !Type of integration strategy to use (INT_TYPE_PCHIP, INT_TYPE_LINEAR or INT_TYPE_STEP)
      real(REAL64), dimension(:), intent(out) :: Si             !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          A.Zadra (Aug 2016)

      !Object
      !          to calculate vertical integrals of a given profile
      !          using PCHIP interpolation

      ! Local variables
      integer :: myType, ni, nk

      ! Set error return status
      status = INT_ERR

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Compute integral profile
      ni = size(yi,dim=1)
      nk = size(yi,dim=2)
      ! Dispatch to integration function
      select case (myType)
      case (INT_TYPE_PCHIP)
         status = integral_pchip(zi,z1i,z2i,ni,nk,yi,sid=Si)
      case (INT_TYPE_LINEAR, INT_TYPE_STEP)
         status = integral_linear(zi,z1i,z2i,myType,ni,nk,yi,sid=Si)
      case DEFAULT
         call msg(MSG_WARNING,'(integral_profile_sum_dble) invalid type of integral requested')
         return
      end select

      ! Error handling
      if (status /= INT_OK) then
         call msg(MSG_WARNING,'(integral_profile_sum_dble) error returned from integration subprogram')
         return
      end if
      
      ! End of subprogram
      return
   end function integral_profile_sum_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_profile_sum(Si,yi,zi,z1i,z2i,type) result(status)
      ! Integrate a quantity over a profile, returning integrated value

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, dimension(:), intent(in) :: z1i                     !Lower limit of integration
      real, dimension(:), intent(in) :: z2i                     !Upper limit of integration
      integer, intent(in), optional :: type                     !Type of integration strategy to use (INT_TYPE_PCHIP, INT_TYPE_LINEAR or INT_TYPE_STEP)
      real, dimension(:), intent(out) :: Si                     !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (May 2017)

      !Object
      !          Wrap call to double-precision calculations.

      ! Local variables
      integer :: myType, ni, nk

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Compute integral profile
      ni = size(yi,dim=1)
      nk = size(yi,dim=2)
      ! Dispatch to integration function
      select case (myType)
      case (INT_TYPE_PCHIP)
         status = integral_pchip(zi,z1i,z2i,ni,nk,yi,sis=Si)
      case (INT_TYPE_LINEAR, INT_TYPE_STEP)
         status = integral_linear(zi,z1i,z2i,myType,ni,nk,yi,sis=Si)
      case DEFAULT
         call msg(MSG_WARNING,'(integral_profile_sum) invalid type of integral requested')
         return
      end select

      ! Error handling
      if (status /= INT_OK) then
         call msg(MSG_WARNING,'(integral_profile_sum) error returned from integration subprogram')
         return
      end if

      ! End of subprogram
      return
   end function integral_profile_sum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_profile_sum_scalbnd_dble(Si,yi,zi,z1i,z2i,type) result(status)
      ! Integrate a quantity over a profile, returning integrated value

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, intent(in) :: z1i                                   !Lower limit of integration
      real, intent(in) :: z2i                                   !Upper limit of integration
      integer, intent(in), optional :: type                     !Type of integration strategy to use (INT_TYPE_PCHIP, INT_TYPE_LINEAR or INT_TYPE_STEP)
      real(REAL64), dimension(:), intent(out) :: Si             !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          A.Zadra (Aug 2016)

      !Object
      !          to calculate vertical integrals of a given profile
      !          using PCHIP interpolation

      ! Local variables
      real, dimension(size(yi,dim=1)) :: vz1i,vz2i
      integer :: myType

      ! Set error return status
      status = INT_ERR

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Compute integral profile
      vz1i = z1i
      vz2i = z2i
      if (integral_profile_sum_dble(Si,yi,zi,vz1i,vz2i,type=myType) /= INT_OK) then
         call msg(MSG_WARNING,'(integral_profile_sum_scalbnd_dble) Error returned by integration subprogram')
         return
      endif

      ! Successful completion
      status = INT_OK
      return
   end function integral_profile_sum_scalbnd_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_profile_sum_scalbnd(Si,yi,zi,z1i,z2i,type) result(status)
      ! Integrate a quantity over a profile, returning integrated value

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, intent(in) :: z1i                                   !Lower limit of integration
      real, intent(in) :: z2i                                   !Upper limit of integration
      integer, intent(in), optional :: type                     !Type of integration strategy to use (INT_TYPE_PCHIP, INT_TYPE_LINEAR or INT_TYPE_STEP)
      real, dimension(:), intent(out) :: Si                     !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (May 2017)

      !Object
      !          Wrap call to double-precision calculations.

      ! Local variables
      real, dimension(size(yi,dim=1)) :: vz1i,vz2i
      integer :: myType

      ! Set error return status
      status = INT_ERR

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Compute integral profile
      vz1i = z1i
      vz2i = z2i
      if (integral_profile_sum(Si,yi,zi,vz1i,vz2i,type=myType) /= INT_OK) then
         call msg(MSG_WARNING,'(integral_profile_sum_scalbnd) Error returned by integration subprogram')
         return
      endif

      ! End of subprogram
      status = INT_OK
      return
   end function integral_profile_sum_scalbnd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_solve_vdep(zci,yi,zi,zidep,a,direc,n,nki, &
#ifdef DEBUG
        fd_unittest, &
#endif
        found) result(status)
      ! Solve an integral equation using PCHIP interpolation
      implicit none
      !@Arguments
      integer, intent(in) :: n,nki
      real, dimension(n,nki), intent(in) :: yi                  !Integrand
      real, dimension(n,nki), intent(in) :: zi                  !Vertical coordinate (monotonic, increasing from nki to 1)
      real, dimension(n), intent(in) :: zidep                   !Departure level for integral
      real, dimension(n), intent(in) :: a                       !Right hand side of integral equation
      integer, intent(in) :: direc                              !Direction of integration (INT_DIR_UP or INT_DIR_DOWN)
      real, dimension(n), intent(out) :: zci                    !Value of (z-zdep) for int^z (g(z')*y(z')*dz') = a
#ifdef DEBUG
      integer, intent(in), optional :: fd_unittest              !Open file descriptor for unit test results
#endif
      logical, dimension(n), intent(out), optional :: found     !Found a solution to the integral equation within the column
      integer :: status                                         !Return status of function
      !@Author A.Zadra (Aug 2016)
      !@Object solve an integral equation using PCHIP interpolation
      !@Notes
      !          The result is zci: the value of z at which the equation
      !                int^z (g(z')*y(z')*dz') = a
      !          is solved, with weight function g given by
      !                     | 1.  , between zidep and zi(1) or zi(nki) depending on direc
      !              g(z) = |
      !                     | 0.  , otherwise

      ! Local variable declaration
      integer :: i,k,istat,wi
      integer, dimension(n) :: k1,k2,k3
      real :: delF,root,dfflow,dffhigh,ffzlow,ffzhigh,w0,w1,w2,hlow,hhigh,fflow,ffhigh,da
      real, dimension(n) :: z1,z2,zc,hlow1,dfflow1,fflow1,hhigh2
      real, dimension(n,nki+PCHIP_EXT) :: F2d
      real, dimension(n,nki+PCHIP_EXT) :: z,y,h,del,b,c,d,d2c
      integer :: k1max,k2min,nko,ko
#ifdef DEBUG
      integer :: ku
      real :: dz,floc,zloc
#endif

      nko = nki+PCHIP_EXT
      
      ! Set error return status
      status = INT_ERR
      
      ! Extend arrays to integration bounds
      ! Set bounds of integral
      ! Invert input arrays if necessary
      IFDIREC: if (direc == PCHIP_DIR_DOWN) then
         do i=1,n
            y(i,nko) = sign(1., a(i))*yi(i,1)
            y(i,1)   = sign(1., a(i))*yi(i,nki)
            z(i,nko) = zi(i,1) - max(zidep(i),  2*zi(i,1)   - zi(i,2))
            z(i,1)   = zi(i,1) - min(zi(i,nki), 2*zi(i,nki) - zi(i,nki-1))
            z1(i)    = zi(i,1) - zidep(i)
            z2(i)    = zi(i,1) - zi(i,nki)
         enddo
         do k=nki,1,-1
            ko = nko-k !# ko=k+PCHIP_EXT/2 ; ko=nko-ko+1
            do i=1,n
               y(i,ko)   = sign(1., a(i))*yi(i,k)
               z(i,ko)   = zi(i,1) - zi(i,k)
               h(i,ko)   = z(i,ko-1) - z(i,ko)
               del(i,ko) = (y(i,ko-1) - y(i,ko)) / h(i,ko)
            enddo
         enddo
      else  !#IFDIREC
         do i=1,n
            y(i,1)   = sign(1., a(i))*yi(i,1)
            y(i,nko) = sign(1., a(i))*yi(i,nki)
            z(i,1)   = max(zi(i,1),  2*zi(i,1)   - zi(i,2))
            z(i,nko) = min(zidep(i), 2*zi(i,nki) - zi(i,nki-1))
            z1(i)    = zidep(i)
            z2(i)    = zi(i,1)
         enddo
         do k=1,nki
            ko = k + 1 !#PCHIP_EXT/2
            do i=1,n
               y(i,ko)  = sign(1., a(i))*yi(i,k)
               z(i,ko)  = zi(i,k)
               h(i,ko)   = z(i,ko-1) - z(i,ko)
               del(i,ko) = (y(i,ko-1) - y(i,ko)) / h(i,ko)
            enddo
         enddo
      endif IFDIREC
           
      do i=1,n
         h(i,1)   = h(i,2)
         del(i,1) = del(i,2)
         h(i,nko) = z(i,nko-1) - z(i,nko)
         del(i,nko) = (y(i,nko-1) - y(i,nko)) / h(i,nko)
      enddo

      ! Compute interpolating polynomial coefficients
      istat = pchip_coef(h,del,b,c,d,n,nko)

      ! Find subdomain
      istat = subdomain(z,z1,z2,k1,k2,n,nko)

      ! Solve equation
      do i=1,n
         k = k1(i)
         hlow1(i)   = z1(i)-z(i,k)
         dfflow1(i) = dff(hlow1(i))
         fflow1(i)  = ff(hlow1(i))
         hhigh2(i)  = z2(i)-z(i,k2(i))
      enddo

      zc = z2
      k3 = -1
      F2d = 0.
      k1max = min(maxval(k1), nko)
      k2min = max(2, minval(k2))
      KLOOP: do k=k1max,k2min,-1
         do i=1,n
            if (k > k1(i) .or. k < k2(i) .or. k3(i) /= -1) cycle
            w1 = float(1 - min(abs(k-k1(i)), 1))
            w2 = float(1 - min(abs(k-k2(i)), 1))
            w0 = 1. - max(w1,w2)
            
            hlow  = w1*hlow1(i)  !# + w2*hlow2(i) + w0*0.
            hhigh = w2*hhigh2(i) + (w0+w1)*h(i,k)
            dfflow   = w1*dfflow1(i) + (w0+w2)*y(i,k)
            fflow = w1*fflow1(i)  !# + (w0+w2)*0.
            F2d(i,k) = F2d(i,k) - w1*fflow1(i)
            
            ffhigh  = ff(hhigh)
            dffhigh = dff(hhigh)
            
            ! Remaining RHS value for this layer
            da  = abs(a(i)) - F2d(i,k)

            ! Check for a possible maximum within the layer
            LOCAL_MAX: if (sign(1.,dfflow) /= sign(1.,dffhigh)) then
               d2c(i,k) = d(i,k) + 2.*c(i,k)
               istat = rf_nrbnd(root,dff,d2ff,hlow,hhigh,dfflow,dffhigh)
#ifdef DEBUG
               if (istat /= RF_OK) then
                  call msg(MSG_WARNING,'(integral_solve_vdep) Cannot find local root')
                  return
               endif
#endif
               delF = ff(root)
               w0 = min(max(sign(1.,da - delF), 0.), 1.)
               hhigh  = w0*hhigh  + (1.-w0)*root
               ffhigh = w0*ffhigh + (1.-w0)*delF
            endif LOCAL_MAX

            ! Integrate over layer
            F2d(i,k-1) = F2d(i,k) + ffhigh

            ! Numerical solution to compute position within the layer
            if (ffhigh >= da) then
               k3(i) = k
               istat = rf_nrbnd(root,ffz,dff,hlow,hhigh,fflow-da,ffhigh-da)
#ifdef DEBUG
               if (istat /= RF_OK) then
                  call msg(MSG_WARNING,'(integral_solve_vdep) Cannot solve integral equation')
                  return
               endif
#endif
               zc(i) = root + z(i,k)
            endif

         enddo
      enddo KLOOP

#ifdef DEBUG
      ! Final check that the solution height remains within bracketing interval
      if (any(zc < z1) .or. any(zc > z2)) then
         call msg(MSG_WARNING,'(integral_solve_vdep) Integral solution out of bounds')
         return
      endif
#endif
      
      ! Re-invert level if necessary
      select case (direc)
      case (PCHIP_DIR_DOWN)
         zci = zidep - (zi(:,1) - zc)
      case (PCHIP_DIR_UP)
         zci = zc - zidep
      end select
 
      ! Generate solution report on request
      if (present(found)) found = (k3 /= -1)

      ! Successful completion
      status = INT_OK
      return

    contains

      ! Internal functions

      real function ff(x)
        ! Interploating polynomial
        real, intent(in) :: x
        ff = x*(y(i,k) + x*(0.5*d(i,k) + x*(TRD*c(i,k) + x*0.25*b(i,k))))
      end function ff

      real function ffz(x)
        ! Interpolating polynomial with 0 on the RHS for root-finding
        real, intent(in) :: x
        ffz = ff(x) - da
       end function ffz

      real function dff(x)
        ! First derivative of interpolating polynomial
        real, intent(in) :: x
        dff = (y(i,k) + x*(d(i,k) + x*(c(i,k) + x*b(i,k))))
      end function dff
     
      real function d2ff(x)
         ! Second derivative of interpolating polynomial
         real, intent(in) :: x
         d2ff = (d2c(i,k) + x*x*3.*b(i,k))
      end function d2ff
 
   end function integral_solve_vdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_solve_sdep(zci,yi,zi,zidep,a,direc) result(status)
      implicit none
      !@Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, intent(in) :: zidep                                 !Departure level for integral
      real, dimension(:), intent(in) :: a                       !Right hand side of integral equation
      integer, intent(in) :: direc                              !Direction of integration (INT_DIR_UP or INT_DIR_DOWN)
      real, dimension(:), intent(out) :: zci                    !Value of (z-zdep) for int^z (g(z')*y(z')*dz') = a
      integer :: status                                         !Return status of function
      !@Author A.Zadra (Aug 2016)
      !@Object solve an integral equation using PCHIP interpolation

      integer :: n,nk
      real, dimension(size(yi,dim=1)) :: vzidep
      
      n = size(yi,dim=1)
      nk = size(yi,dim=2)
      vzidep = zidep
      status = integral_solve_vdep(zci,yi,zi,vzidep,a,direc,n,nk)
      if (status /= INT_OK) then
         call msg(MSG_WARNING,'(integral_solve_sdep) Error returned by integral solver subprogram')
      endif

      return
   end function integral_solve_sdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_solve_sdepsa(zci,yi,zi,zidep,a,direc) result(status)
      implicit none
      !@Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, intent(in) :: zidep                                 !Departure level for integral
      real, intent(in) :: a                                     !Right hand side of integral equation
      integer, intent(in) :: direc                              !Direction of integration (INT_DIR_UP or INT_DIR_DOWN)
      real, dimension(:), intent(out) :: zci                    !Value of (z-zdep) for int^z (g(z')*y(z')*dz') = a
      integer :: status                                         !Return status of function
      !@Author A.Zadra (Aug 2016)
      !@Object solve an integral equation using PCHIP interpolation

      real, dimension(size(yi,dim=1)) :: va

      va = a
      status = integral_solve_sdep(zci,yi,zi,zidep,va,direc)
      if (status /= INT_OK) then
         call msg(MSG_WARNING,'(integral_solve_sdepsa) Error returned by integral solver subprogram')
      endif

      return
   end function integral_solve_sdepsa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dir$ attributes forceinline :: subdomain
   function subdomain(z,z1,z2,k1,k2,n,nk) result(status)
      !@Object Find subdomain
      !@Arguments
      integer, intent(in) :: n,nk
      real, dimension(n,nk), intent(in) :: z
      real, dimension(n), intent(in) :: z1,z2
      integer, dimension(n), intent(out) :: k1,k2
      !@Author A.Zadra (Aug 2016)

      integer :: status

      ! Local variables
      integer :: i,k

      ! Set return status
      status = INT_ERR

      ! Find subdomain
      k1 = 2
      k2 = 2
      do k = nk,1,-1
         do i = 1,n
            if (z(i,k) > z1(i)) then
               k1(i) = min(max(k+1, k1(i)), nk)
            endif
            if (z(i,k) >= z2(i)) then
               k2(i) = max(k, k2(i))
            endif
         enddo
      enddo

      ! Successful completion
      status = INT_OK
      return
   end function subdomain

   
!dir$ attributes forceinline :: subdomain2
   function subdomain2(z,z1,z2,k1,k2,n,nk) result(status)
      !@Object Find subdomain
      !@Arguments
      integer, intent(in) :: n,nk
      real, dimension(n,nk), intent(in) :: z
      real, dimension(n), intent(in) :: z1,z2
      integer, dimension(n), intent(out) :: k1,k2

      integer :: status

      ! Local variables
      integer :: i,k

      ! Set return status
      status = INT_ERR

      ! Find subdomain
      k1 = 1
      k2 = 1
      do k = nk,2,-1
         do i = 1,n
            if (z(i,k) >= z2(i) .and. k1(i) /= 1) then
               k2(i) = max(k2(i), k)
            endif
            if (z(i,k) >= z1(i)) then
               k1(i) = max(k1(i), k)
            endif
         enddo
      enddo

      ! Successful completion
      status = INT_OK
      return
   end function subdomain2
   
end module integrals
