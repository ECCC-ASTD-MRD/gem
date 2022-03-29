!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------

subroutine hines_wind(v_alpha,vel_u,vel_v,  &
     &                naz,il1,il2,lev1,lev2,nlons,nlevs,nazmth)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer :: naz, il1, il2, lev1, lev2
   integer :: nlons, nlevs, nazmth
   real(REAL64) :: v_alpha(nlons,nlevs,nazmth)
   real(REAL64) ::  vel_u(nlons,nlevs), vel_v(nlons,nlevs)

   !@Author
   !  aug. 7/95 - c. mclandress
   !@Object
   !  This routine calculates the azimuthal horizontal background wind component
   !  on a longitude by altitude grid for the case of 4 or 8 azimuths for
   !  the hines' doppler spread gwd parameterization scheme.
   !@Arguments
   !          - Output -
   ! v_alpha  background wind component at each azimuth (m/s).
   !             (note: first azimuth is in eastward direction
   !              and rotate in counterclockwise direction.)
   !          - Input -
   ! vel_u    background zonal wind component (m/s).
   ! vel_v    background meridional wind component (m/s).
   ! naz      actual number of horizontal azimuths used (must be 4 or 8).
   ! il1      first longitudinal index to use (il1 >= 1).
   ! il2      last longitudinal index to use (il1 <= il2 <= nlons).
   ! lev1     first altitude level to use (lev1 >=1).
   ! lev2     last altitude level to use (lev1 < lev2 <= nlevs).
   ! nlons    number of longitudes.
   ! nlevs    number of vertical levels.
   ! nazmth   azimuthal array dimension (nazmth >= naz).

   !          - constants in data statements -
   ! cos45    cosine of 45 degrees.
   ! umin     minimum allowable value for zonal or meridional
   !          wind component (m/s).


   !  internal variables.

   integer  i, l
   real(REAL64) :: u, v, cos45, umin
   !-----------------------------------------------------------------------

   cos45 = 0.7071068
   umin  = 0.001

   select case (naz)

   case(4)  !  case with 4 azimuths.

      do l = lev1,lev2
         do i = il1,il2
            u = vel_u(i,l)
            v = vel_v(i,l)
            if (abs(u) .lt. umin)  u = umin
            if (abs(v) .lt. umin)  v = umin
            v_alpha(i,l,1) = u
            v_alpha(i,l,2) = v
            v_alpha(i,l,3) = - u
            v_alpha(i,l,4) = - v
         end do
      end do

   case (8)   !  case with 8 azimuths.

      do l = lev1,lev2
         do i = il1,il2
            u = vel_u(i,l)
            v = vel_v(i,l)
            if (abs(u) .lt. umin)  u = umin
            if (abs(v) .lt. umin)  v = umin
            v_alpha(i,l,1) = u
            v_alpha(i,l,2) = cos45 * ( v + u )
            v_alpha(i,l,3) = v
            v_alpha(i,l,4) = cos45 * ( v - u )
            v_alpha(i,l,5) = - u
            v_alpha(i,l,6) = - v_alpha(i,l,2)
            v_alpha(i,l,7) = - v
            v_alpha(i,l,8) = - v_alpha(i,l,4)
         end do
      end do

   end select

   !-----------------------------------------------------------------------
   return
end subroutine hines_wind
