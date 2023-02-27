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
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_toupper
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
 
   ! Internal parameters
   integer, parameter :: LONG_CHAR=16                           !Long character string length
   character(len=LONG_CHAR), parameter :: DEFAULT_TYPE='pchip'  !Default integration method

   ! Define API
   integer, parameter, public :: INT_OK=0,INT_ERR=1             !Return statuses of functions
   public :: int_profile                                        !Compute integrals for a profile
   public :: int_solve                                          !Solve an integral equation

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_pchip(Ai,yi,zi,z1i,z2i) result(status)
      use pchip_utils, only: PCHIP_OK,PCHIP_EXT,pchip_extend,pchip_layers,pchip_coef
      implicit none

      !@Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, dimension(:), target, intent(in) :: z1i             !Lower limit of integration
      real, dimension(:), target, intent(in) :: z2i             !Upper limit of integration
      real(REAL64), dimension(:,:), intent(out) :: Ai     !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                                !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      integer :: status                                         !Return status of function

      !@Author A.Zadra (Aug 2016)
      !@Object calculate vertical integrals of a given profile
      !        using PCHIP interpolation
      !        Integrate a quantity over a profile using PCHIP,
      !        returning integrated quantity along the profile

      ! Local variable declaration
      integer :: n,nk,i,k,istat
      integer, dimension(size(yi,dim=1)) :: k1,k2
      real :: hk,yk,dk,ck,bk,r,s
      real, dimension(:), pointer :: z1iloc,z2iloc
      real, dimension(size(yi,dim=1)) :: z1,z2
      real, dimension(size(yi,dim=1),size(yi,dim=2)+PCHIP_EXT) :: z,y,h,del,b,c,d
      real(REAL64), dimension(size(yi,dim=1),size(yi,dim=2)+PCHIP_EXT) :: A
      logical :: bnd_invert

      ! Set error return status
      status = INT_ERR

      ! Initialize bounds
      n = size(yi,dim=1)
      nk = size(yi,dim=2)

      ! Determine direction of integration
      if (any(z1i<=z2i)) then
         z1iloc => z1i
         z2iloc => z2i
         bnd_invert = .false.
      else
         z1iloc => z2i
         z2iloc => z1i
         bnd_invert = .true.
      endif
 
      ! Invert input arrays if necessary
      if (pchip_extend(yi,zi,y,z,z1i=z1iloc,z2i=z2iloc,z1=z1,z2=z2,direc='UP') /= PCHIP_OK) then
         call msg(MSG_WARNING,'(integral_pchip) Error returned by pchip_extend')
         return
      endif

      ! Compute layer thickness and deltas
      istat = pchip_layers(y,z,h,del)

      ! Compute interpolating polynomial coefficients
      istat = pchip_coef(h,del,b,c,d)

      ! Find subdomain
      istat = subdomain(z,z1,z2,k1,k2)

      ! Integral profile by interpolation
      A = 0.
      do i=1,n
         do k=k1(i),k2(i),-1
            hk = h(i,k)
            yk = y(i,k)
            dk = d(i,k)
            ck = c(i,k)
            bk = b(i,k)

            A(i,k-1) = A(i,k) + hk*yk &
                 + (hk**2)*dk/2. &
                 + (hk**3)*ck/3. &
                 + (hk**4)*bk/4.

            if (k == k1(i)) then
               r = z1(i) - z(i,k)
               A(i,k-1) = A(i,k-1) - ( r*yk &
                    + (r**2)*dk/2. &
                    + (r**3)*ck/3. &
                    + (r**4)*bk/4. )
            endif
            if (k == k2(i)) then
               s = z2(i) - z(i,k)
               A(i,k-1) = A(i,k-1) + (s - hk)*yk &
                    + (s**2 - hk**2)*dk/2. &
                    + (s**3 - hk**3)*ck/3. &
                    + (s**4 - hk**4)*bk/4.
            endif
         enddo
         do k=k1(i)+1,nk
            A(i,k) = A(i,k1(i))
         enddo
         do k=k2(i)-2,1,-1
            A(i,k) = A(i,k2(i)-1)
         enddo
      enddo

      ! Return interior of domain only
      ! Negate result if bounds were inverted
      if (bnd_invert) then
         Ai = -A(:,2:nk+1)
      else
         Ai = A(:,2:nk+1)
      endif

      ! Successful completion
      status = INT_OK
      return
   end function integral_pchip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_linear(Ai,yi,zi,z1i,z2i,type) result(status)
      ! Integrate a quantity over a profile using linear interpolation, returning integrated quantity along the profile
      use pchip_utils, only: PCHIP_OK,PCHIP_EXT,pchip_extend
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, dimension(:), target, intent(in) :: z1i             !Lower limit of integration
      real, dimension(:), target, intent(in) :: z2i             !Upper limit of integration
      character(len=*), intent(in) :: type                      !Linear subtype ('linear' or 'step')
      real(REAL64), dimension(:,:), intent(out) :: Ai     !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                                !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      integer :: status                                         !Return status of function

      !Author
      !          R. McTaggart-Cowan and A.Zadra (July 2017)

      !Object
      !          to calculate vertical integrals of a given profile
      !          using linear interpolation

      ! Local variable declaration
      integer :: n,nk,nk2,i,k
      real :: zmid, zmid_step_mask
      real, dimension(:), pointer :: z1iloc,z2iloc
      real, dimension(size(yi,dim=1)) :: z1,z2
      real, dimension(size(yi,dim=1),size(yi,dim=2)+2) :: z,y,zh
      real(REAL64), dimension(size(yi,dim=1),size(yi,dim=2)+2) :: A
      logical :: bnd_invert

      ! Set error return status
      status = INT_ERR

      ! Initialize bounds
      n = size(yi,dim=1)
      nk = size(yi,dim=2)

      ! Determine direction of integration
      if (any(z1i<=z2i)) then
         z1iloc => z1i
         z2iloc => z2i
         bnd_invert = .false.
      else
         z1iloc => z2i
         z2iloc => z1i
         bnd_invert = .true.
      endif

      ! Invert input arrays if necessary
      if (pchip_extend(yi,zi,y,z,z1i=z1iloc,z2i=z2iloc,z1=z1,z2=z2,direc='UP') /= PCHIP_OK) then
         call msg(MSG_WARNING,'(integral_linear) Error returned by pchip_extend')
         return
      endif
 
      ! Define layer interfaces
      zh = 0.
      nk2 = size(z,dim=2)
      do k=nk2,2,-1
         zh(:,k) = (z(:,k) + z(:,k-1)) / 2.
      enddo

      ! Integrate profile
      A = 0.
      zmid_step_mask = 1.
      if (type == 'STEP') zmid_step_mask = 0.
      do i=1,n
         k = nk2
         do while (zh(i,k) < z1(i) .and. k > 1)
            k = k-1
         enddo
         zmid = 0.5*(zh(i,k)+z1(i))
         zmid = zmid_step_mask * zmid + (1.-zmid_step_mask) * z(i,k)
         A(i,k) = (y(i,k) + (y(i,k-1)-y(i,k))*((zmid-z(i,k))/(z(i,k-1)-z(i,k)))) * (zh(i,k)-z1(i))
         do while (zh(i,k) < z2(i) .and. k > 1)
            k = k-1
            A(i,k) = A(i,k+1) + y(i,k) * (zh(i,k)-zh(i,k+1))
         enddo
         zmid = 0.5*(zh(i,k+1)+z2(i))
         zmid = zmid_step_mask * zmid + (1.-zmid_step_mask) * z(i,k)
         A(i,k) = A(i,k+1) + (y(i,k) - (y(i,k)-y(i,k+1))*((z(i,k)-zmid)/(z(i,k)-z(i,k+1)))) * (z2(i)-zh(i,k+1))
         do while (k > 1)
            k = k-1
            A(i,k) = A(i,k+1)
         enddo
      enddo

      ! Return interior of domain only
      if (bnd_invert) then
         Ai = -A(:,2:nk+1)
      else
         Ai = A(:,2:nk+1)
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
      character(len=*), intent(in), optional :: type            !Type of integration strategy to use (['pchip'],'linear','step')
      real(REAL64), dimension(:,:), intent(out) :: Ai     !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                                !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (July 2017)

      !Object
      !          Dispatch requests for vertical integration.

      ! Local variables
      character(len=LONG_CHAR) :: myType
      integer :: istat
 
      ! Initialize return value
      status = INT_ERR

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Dispatch to integration function
      istat = clib_toupper(myType)
      select case (myType)
      case ('PCHIP')
         status = integral_pchip(Ai,yi,zi,z1i,z2i)
      case ('LINEAR','STEP')
         status = integral_linear(Ai,yi,zi,z1i,z2i,myType)
      case DEFAULT
         call msg(MSG_WARNING,'(integral_profile_vec_dble) invalid type of integral ('//trim(myType)//') requested')
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
      character(len=*), intent(in), optional :: type            !Type of integration strategy to use (['pchip'],'linear','step')
      real, dimension(:,:), intent(out) :: Ai                   !Integrated profile A(z) = int^z (g(z')*y(z')*dz')
                                                                !   where g(z)=1 for z1i<=z<=z2i (0 elsewhere)
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (May 2017)

      !Object
      !          Wrap call to double-precision calculation.

      ! Local variables
      real(REAL64), dimension(size(Ai,dim=1),size(Ai,dim=2)) :: Ai_dble
      character(len=LONG_CHAR) :: myType

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Call full-precision calculation
      status = integral_profile_vec_dble(Ai_dble,yi,zi,z1i,z2i,type=myType)
      Ai = real(Ai_dble)

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
      character(len=*), intent(in), optional :: type            !Type of integration strategy to use (['pchip'],'linear','step')
      real(REAL64), dimension(:), intent(out) :: Si       !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          A.Zadra (Aug 2016)

      !Object
      !          to calculate vertical integrals of a given profile
      !          using PCHIP interpolation

      ! Local variables
      real(REAL64), dimension(size(yi,dim=1),size(yi,dim=2)) :: Ai
      character(len=LONG_CHAR) :: myType

      ! Set error return status
      status = INT_ERR

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Compute integral profile
      if (integral_profile_vec_dble(Ai,yi,zi,z1i,z2i,type=myType) /= INT_OK) then
         call msg(MSG_WARNING,'(integral_profile_sum_dble) error returned by integration subprogram')
         return
      endif

      ! Select end-point for integrated value
      Si = Ai(:,1)

      ! Successful completion
      status = INT_OK
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
      character(len=*), intent(in), optional :: type            !Type of integration strategy to use (['pchip'],'linear','step')
      real, dimension(:), intent(out) :: Si                     !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (May 2017)

      !Object
      !          Wrap call to double-precision calculations.

      ! Local variables
      real(REAL64), dimension(size(Si)) :: Si_dble
      character(len=LONG_CHAR) :: myType

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Call full-precision calculation
      status = integral_profile_sum_dble(Si_dble,yi,zi,z1i,z2i,type=myType)
      Si = real(Si_dble)

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
      character(len=*), intent(in), optional :: type            !Type of integration strategy to use (['pchip'],'linear','step')
      real(REAL64), dimension(:), intent(out) :: Si       !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          A.Zadra (Aug 2016)

      !Object
      !          to calculate vertical integrals of a given profile
      !          using PCHIP interpolation

      ! Local variables
      real, dimension(size(yi,dim=1)) :: vz1i,vz2i
      character(len=LONG_CHAR) :: myType

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
      character(len=*), intent(in), optional :: type            !Type of integration strategy to use (['pchip'],'linear','step')
      real, dimension(:), intent(out) :: Si                     !Integrated value
      integer :: status                                         !Return status of function

      !Author
      !          R.McTaggart-Cowan (May 2017)

      !Object
      !          Wrap call to double-precision calculations.

      ! Local variables
      real(REAL64), dimension(size(Si)) :: Si_dble
      character(len=LONG_CHAR) :: myType

      ! Handle optional arguments
      myType = DEFAULT_TYPE
      if (present(type)) myType = type

      ! Call full-precision calculation
      status = integral_profile_sum_scalbnd_dble(Si_dble,yi,zi,z1i,z2i,type=myType)
      Si = real(Si_dble)

      ! End of subprogram
      return
   end function integral_profile_sum_scalbnd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integral_solve_vdep(zci,yi,zi,zidep,a,direc,fd_unittest,found) result(status)
      ! Solve an integral equation using PCHIP interpolation
      use pchip_utils, only: PCHIP_OK,PCHIP_EXT,pchip_extend,pchip_layers,pchip_coef
      use rootfind, only: RF_OK,rf_nrbnd
      implicit none
      !@Arguments
      real, dimension(:,:), intent(in) :: yi                    !Integrand
      real, dimension(:,:), intent(in) :: zi                    !Vertical coordinate (monotonic, increasing from nk to 1)
      real, dimension(:), intent(in) :: zidep                   !Departure level for integral
      real, dimension(:), intent(in) :: a                       !Right hand side of integral equation
      character(len=*), intent(in) :: direc                     !Direction of integration ('up', 'down')
      real, dimension(:), intent(out) :: zci                    !Value of (z-zdep) for int^z (g(z')*y(z')*dz') = a
      integer, intent(in), optional :: fd_unittest              !Open file descriptor for unit test results
      logical, dimension(:), intent(out), optional :: found     !Found a solution to the integral equation within the column
      integer :: status                                         !Return status of function
      !@Author A.Zadra (Aug 2016)
      !@Object solve an integral equation using PCHIP interpolation
      !@Notes
      !          The result is zci: the value of z at which the equation
      !                int^z (g(z')*y(z')*dz') = a
      !          is solved, with weight function g given by
      !                     | 1.  , between zidep and zi(1) or zi(nk) depending on direc
      !              g(z) = |
      !                     | 0.  , otherwise

      ! Local variable declaration
      integer :: n,nk,i,k,istat
      integer, dimension(size(yi,dim=1)) :: k1,k2
      real :: da,delF,hlow,hhigh
      real, dimension(1) :: root
      real, dimension(size(yi,dim=1)) :: z1,z2,zc,z1i,z2i,signa
      real, dimension(size(yi,dim=2)+PCHIP_EXT) :: F
      real, dimension(size(yi,dim=1),size(yi,dim=2)) :: yyi
      real, dimension(size(yi,dim=1),size(yi,dim=2)+PCHIP_EXT) :: z,y,h,del,b,c,d
      logical, dimension(size(yi,dim=1)) :: myFound
      character(len=LONG_CHAR) :: myDirec
#ifdef DEBUG
      integer :: ku
      real :: dz,floc,zloc
#endif
      
      ! Set error return status
      status = INT_ERR

      ! Initialize array bounds
      n = size(yi,dim=1)
      nk = size(yi,dim=2)
 
      ! For the integral solution a() should always be positive, so negate yi() if necessary
      signa(:) = sign(1.,a(:))
      do k=1,nk
         do i=1,n
            yyi(i,k) = signa(i) * yi(i,k)
         enddo
      enddo

      ! Set bounds of integral
      myDirec = direc
      istat = clib_toupper(myDirec)
      select case (myDirec)
      case ('DOWN')
         z1i = zi(:,nk)
         z2i = zidep
      case ('UP')
         z1i = zidep
         z2i = zi(:,1)
      case DEFAULT
         z1i = -1.
         z2i = -1.
      end select

      ! Invert input arrays if necessary
      if (pchip_extend(yyi,zi,y,z,z1i=z1i,z2i=z2i,z1=z1,z2=z2,direc=direc) /= PCHIP_OK) then
         call msg(MSG_WARNING,'(integral_solve_vdep) Error returned by pchip_extend')
         return
      endif

      ! Compute layer thickness and deltas
      istat = pchip_layers(y,z,h,del)

      ! Compute interpolating polynomial coefficients
      istat = pchip_coef(h,del,b,c,d)

      ! Find subdomain
      istat = subdomain(z,z1,z2,k1,k2)

      ! Solve equation
      zc = z2
      do i=1,n
         myFound(i) = .false.
         F = 0.
         KLOOP: do k=k1(i),k2(i),-1

            ! Set up integration bounds for polyomial in this layer
            if (k == k1(i)) then
               hlow = z1(i)-z(i,k)
               hhigh = h(i,k)
               F(k) = -ff(hlow)
            elseif (k == k2(i)) then
               hlow = 0.
               hhigh = z2(i)-z(i,k)
            else
               hlow = 0.
               hhigh = h(i,k)
            endif

            ! Remaining RHS value for this layer
            da  = abs(a(i)) - F(k)

#ifdef DEBUG
            ! Numerical integration for unit testing
            if (present(fd_unittest)) then
               do ku=1,50
                  dz = (ku-1)/50.*(hhigh-hlow)
                  zloc = z(i,k)
                  floc = F(k)+ff(hlow+dz)
                  dz = hlow + dz
                  if (myDirec == 'DOWN') then
                     zloc = z(i,1)-zloc
                     dz = -dz
                  endif
                  if (a(i) < 0.) floc = -floc
                  write(fd_unittest, '(2(e,x))') zloc+dz, floc
               enddo
            endif
#endif

            ! Check for a possible maximum within the layer
            LOCAL_MAX: if (dff(hlow)*dff(hhigh) < 0.) then
               if (rf_nrbnd(root,dff,d2ff,(/hlow/),(/hhigh/)) /= RF_OK) then
                  call msg(MSG_WARNING,'(integral_solve_vdep) Cannot find local root')
                  return
               endif
               delF = ff(root(1))
               if (delF >= da) hhigh = root(1)
            endif LOCAL_MAX

            ! Integrate over layer
            delF = ff(hhigh)
            F(k-1) = F(k) + delF

            ! Numerical solution to compute position within the layer
            if (delF >= da) then
               myFound(i) = .true.
               if (rf_nrbnd(root,ffz,dff,(/hlow/),(/hhigh/)) /= RF_OK) then
                  call msg(MSG_WARNING,'(integral_solve_vdep) Cannot solve integral equation')
                  return
               endif
               zc(i) = root(1) + z(i,k)
            endif

            if (myFound(i)) exit KLOOP
         enddo KLOOP

         ! Final check that the solution height remains within bracketing interval
         if (zc(i) < z1(i) .or. zc(i) > z2(i)) then
            call msg(MSG_WARNING,'(integral_solve_vdep) Integral solution out of bounds')
            return
         endif
      enddo

      ! Re-invert level if necessary
      select case (myDirec)
      case ('DOWN')
         zci = zidep - (zi(:,1) - zc)
      case ('UP')
         zci = zc - zidep
      end select
 
      ! Generate solution report on request
      if (present(found)) found(:) = myFound(:)

      ! Successful completion
      status = INT_OK
      return

    contains

      ! Internal functions

      real function ff(x)
        ! Interploating polynomial
        real, intent(in) :: x
        ff = x*y(i,k) + (d(i,k)/2.)*x**2. + (c(i,k)/3.)*x**3. + (b(i,k)/4.)*x**4
      end function ff

      real function ffz(x)
        ! Interpolating polynomial with 0 on the RHS for root-finding
        real, intent(in) :: x
        ffz = ff(x) - da
      end function ffz

      real function dff(x)
        ! First derivative of interpolating polynomial
        real, intent(in) :: x
        dff = y(i,k) + d(i,k)*x + c(i,k)*x**2 + b(i,k)*x**3
      end function dff

      real function d2ff(x)
         ! Second derivative of interpolating polynomial
         real, intent(in) :: x
         d2ff = d(i,k) + 2.*c(i,k) + 3*b(i,k)*x**2
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
      character(len=*), intent(in) :: direc                     !Direction of integration ('up', 'down')
      real, dimension(:), intent(out) :: zci                    !Value of (z-zdep) for int^z (g(z')*y(z')*dz') = a
      integer :: status                                         !Return status of function
      !@Author A.Zadra (Aug 2016)
      !@Object solve an integral equation using PCHIP interpolation

      real, dimension(size(yi,dim=1)) :: vzidep

      vzidep = zidep
      status = integral_solve_vdep(zci,yi,zi,vzidep,a,direc)
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
      character(len=*), intent(in) :: direc                     !Direction of integration ('up', 'down')
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
  function subdomain(z,z1,z2,k1,k2) result(status)
      !@Object Find subdomain
      !@Arguments
      real, dimension(:,:), intent(in) :: z
      real, dimension(:), intent(in) :: z1,z2
      integer, dimension(:), intent(out) :: k1,k2
      !@Author A.Zadra (Aug 2016)

      integer :: status

      ! Local variables
      integer :: n,nk,i,k

      ! Set return status
      status = INT_ERR

      ! Find subdomain
      n = size(z,dim=1)
      nk = size(z,dim=2)
      do i=1,n
         k = nk
         do while (z1(i) >= z(i,k) .and. k > 1)
            k = k-1
         enddo
         k1(i) = min(nk, k+1)
         do while (z2(i) >= z(i,k) .and. k > 1)
            k = k-1
         enddo
         k2(i) = min(nk, k+1)
      enddo

      ! Successful completion
      status = INT_OK
      return
   end function subdomain

end module integrals
