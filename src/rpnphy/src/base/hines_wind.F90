
module hines_wind
   implicit none
   private
   public :: hines_wind1

contains

subroutine hines_wind1(v_alpha, vel_u, vel_v, levbot, ni, nig, nkm1, naz)
   use, intrinsic :: iso_fortran_env, only: REAL32
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: levbot, ni, nig, nkm1, naz
   real(REAL32), intent(out) :: v_alpha(ni,naz,nkm1)
   real(REAL32), intent(in) :: vel_u(ni,nkm1), vel_v(ni,nkm1)

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
   ! levbot   last altitude level to use (1 < levbot <= nkm1).
   ! ni       number of longitudes.
   ! nig      horizontal operator scope
   ! nkm1     number of vertical levels.
   ! naz      azimuthal array dimension (naz >= naz).

   !          - constants in data statements -
   ! cos45    cosine of 45 degrees.
   ! umin     minimum allowable value for zonal or meridional
   !          wind component (m/s).


   !  internal variables.

   real(REAL32), parameter :: cos45 = 0.7071068 !# cos(45) computed by compiler?
   real(REAL32), parameter :: umin = 0.001

   integer :: i, l
   real(REAL32) :: u, v
   !-----------------------------------------------------------------------
   
   do l = 1,levbot
      do i = 1,nig
         u = vel_u(i,l)
         v = vel_v(i,l)
         if (abs(u) < umin)  u = umin  !# u = sign(max(umin, abs(u)), u)
         if (abs(v) < umin)  v = umin  !# v = sign(max(umin, abs(v)), v)
         v_alpha(i,1,l) = u
         v_alpha(i,2,l) = cos45 * ( v + u )
         v_alpha(i,3,l) = v
         v_alpha(i,4,l) = cos45 * ( v - u )
         v_alpha(i,5,l) = - u
         v_alpha(i,6,l) = - v_alpha(i,2,l)
         v_alpha(i,7,l) = - v
         v_alpha(i,8,l) = - v_alpha(i,4,l)
      end do
   end do

   !-----------------------------------------------------------------------
   return
end subroutine hines_wind1

end module hines_wind
