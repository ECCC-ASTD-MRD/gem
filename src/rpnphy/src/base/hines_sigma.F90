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

module hines_sigma
   implicit none
   private
   public :: hines_sigma1

contains

subroutine hines_sigma1(sigma_t, sigma_alpha, sigsqh_alpha,  &
     &                  lev, ni, nig, nkm1, naz)
   use, intrinsic :: iso_fortran_env, only: REAL32
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: lev, ni, nig, nkm1, naz
   real(REAL32), intent(out) :: sigma_t(ni,nkm1)
   real(REAL32), intent(out) :: sigma_alpha(ni,naz,nkm1)
   real(REAL32), intent(in)  :: sigsqh_alpha(ni,naz)

   !@Author
   !  aug. 7/95 - c. mclandress
   !@Object
   !  This routine calculates the total rms and azimuthal rms horizontal
   !  velocities at a given level on a longitude by altitude grid for
   !  the hines' doppler spread gwd parameterization scheme.
   !  note: only four or eight azimuths can be used.
   !  output arguements:
   !@Arguments
   !               - Output -
   ! sigma_t       total rms horizontal wind (m/s).
   ! sigma_alpha   total rms wind in each azimuth (m/s).
   !               - Input -
   ! sigsqh_alpha  portion of wind variance from waves having wave
   !               normals in the alpha azimuth (m/s).
   ! lev           altitude level to process.
   ! ni            number of longitudes.
   ! nig           horizontal operator scope
   ! nkm1          number of vertical levels.
   ! naz           azimuthal array dimension (naz >= naz).

   !  internal variables.
   real(REAL32), parameter :: ZERO = 0.
   
   integer :: i, n
   real(REAL32) :: sum_even, sum_odd, above_zero
   !-----------------------------------------------------------------------

   !  calculate azimuthal rms velocity for the 8 azimuth case.
   
   do i = 1,nig
      sum_odd  = (sigsqh_alpha(i,1) + sigsqh_alpha(i,3)   &
           &             + sigsqh_alpha(i,5) + sigsqh_alpha(i,7)) / 2.
      sum_even = (sigsqh_alpha(i,2) + sigsqh_alpha(i,4)   &
           &             + sigsqh_alpha(i,6) + sigsqh_alpha(i,8)) / 2.
      above_zero= max(ZERO, sigsqh_alpha(i,1) + sigsqh_alpha(i,5) + sum_even)
      sigma_alpha(i,1,lev) = sqrt(above_zero)
      ! this next sqrt was indicated to produce FP invalid operation by xlf13
      ! on at least one point on processor oe-0000-0385 out of 480 in a GY25 grid
      above_zero= max(ZERO, sigsqh_alpha(i,2) + sigsqh_alpha(i,6) + sum_odd)
      sigma_alpha(i,2,lev) = sqrt(above_zero)
      above_zero= max(ZERO, sigsqh_alpha(i,3) + sigsqh_alpha(i,7) + sum_even)
      sigma_alpha(i,3,lev) = sqrt(above_zero)
      above_zero= max(ZERO, sigsqh_alpha(i,4) + sigsqh_alpha(i,8) + sum_odd)
      sigma_alpha(i,4,lev) = sqrt(above_zero)
      sigma_alpha(i,5,lev) = sigma_alpha(i,1,lev)
      sigma_alpha(i,6,lev) = sigma_alpha(i,2,lev)
      sigma_alpha(i,7,lev) = sigma_alpha(i,3,lev)
      sigma_alpha(i,8,lev) = sigma_alpha(i,4,lev)
   end do

   !  calculate total rms velocity.
   do i = 1,nig
      sigma_t(i,lev) = 0.
   end do
   do n = 1,naz
      do i = 1,nig
         sigma_t(i,lev) = sigma_t(i,lev) + sigsqh_alpha(i,n)
      end do
   end do
   do i = 1,nig
      above_zero = max(ZERO, sigma_t(i,lev))
      sigma_t(i,lev) = sqrt(above_zero)
   end do

   !-----------------------------------------------------------------------
   return
end subroutine hines_sigma1

end module hines_sigma
