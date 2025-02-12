
module hines_flux
   implicit none
   private
   public :: hines_flux4

contains

subroutine hines_flux4(drag_u, drag_v,        &
     &                 alt, density, densb,   &
     &                 m_alpha, ak_alpha,     &
     &                 levbot, ni, nig, nkm1, naz, &
     &                 hflt)
   use, intrinsic :: iso_fortran_env, only: REAL32
   use mo_gwspectrum, only: m_min, kstar
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: levbot, hflt, ni, nig, nkm1, naz
   real, intent(out) :: drag_u(ni,nkm1), drag_v(ni,nkm1)
   real(REAL32), intent(in) :: alt(ni,nkm1), density(ni,nkm1), densb(ni)
   real(REAL32), intent(in) :: m_alpha(ni,naz,nkm1)
   real(REAL32), intent(in) :: ak_alpha(ni,naz)

   !@Authors
   !  aug. 6/95 - c. mclandress
   !       2001 - m. charron
   !@Revision
   ! 001   L. Spacek (Sep 2008) - calculate drag using centered differences
   ! 002   M. Charron & A. Zadra (June 2011) - use closest levels for centered differences
   !@Object
   !  Calculate zonal and meridional components of the vertical flux
   !  of horizontal momentum and corresponding wave drag (force per unit mass)
   !  on a longitude by altitude grid for the hines' doppler spread
   !  gwd parameterization scheme.
   !  note: only 4 or 8 azimuths can be used.
   !            - Output -
   ! drag_u     zonal component of drag (m/s^2).
   ! drag_v     meridional component of drag (m/s^2).
   !            - Input -
   ! alt        altitudes (m).
   ! density    background density (kg/m^3).
   ! densb      background density at bottom level (kg/m^3).
   ! m_alpha    cutoff vertical wavenumber (1/m).
   ! ak_alpha   spectral amplitude factor (i.e., {ajkj} in m^4/s^2).
   ! m_min      minimum allowable cutoff wavenumber (1/m)
   !            for spectral slope of one.
   ! levbot     last altitude level to use (1 < levbot <= nkm1).
   ! ni         number of longitudes.
   ! nig        horizontal operator scope
   ! nkm1       number of vertical levels.
   ! naz        azimuthal array dimension (naz >= naz).

   !  constant in data statement.

   ! cos45      cosine of 45 degrees.
   real(REAL32), parameter :: cos45 = 0.7071068 !# cos(45) computed by compiler?
   
   integer :: i, l, n, it
   real(REAL32) :: dendz, dendz2, kstar8
   real(REAL32) :: flux(ni,naz)
   real(REAL32) :: flux_u(ni,nkm1) ! zonal momentum flux (pascals)
   real(REAL32) :: flux_v(ni,nkm1) ! meridional momentum flux (pascals)
   real :: work_u(ni,nkm1), work_v(ni,nkm1)
   !-----------------------------------------------------------------------
   kstar8 = kstar
   
   !  calculate flux from sum.
   do l = 1,levbot
      do n = 1, naz
         do i = 1,nig
            flux(i,n) = ak_alpha(i,n)*kstar8*(m_alpha(i,n,l)-m_min)
         end do
      end do
      do i = 1,nig
         flux_u(i,l) = flux(i,1) - flux(i,5) + cos45 *     &
              (flux(i,2) - flux(i,4) - flux(i,6) + flux(i,8))
         flux_v(i,l) = flux(i,3) - flux(i,7) + cos45 *     &
              (flux(i,2) + flux(i,4) - flux(i,6) - flux(i,8))
         flux_u(i,l) = flux_u(i,l) * densb(i)
         flux_v(i,l) = flux_v(i,l) * densb(i)
      end do
   end do

   !  filter fluxes (hflt iterations)
   do it=1,hflt
      do i = 1,nig
         work_u(i,1) = 0.25*(2.*flux_u(i,1) + flux_u(i,2) )
         work_v(i,1) = 0.25*(2.*flux_v(i,1) + flux_v(i,2) )
         work_u(i,levbot) = 0.25*( flux_u(i,levbot-1) + 2.*flux_u(i,levbot) )
         work_v(i,levbot) = 0.25*( flux_v(i,levbot-1) + 2.*flux_v(i,levbot) )
      end do
      do l = 2,levbot-1
         do i = 1,nig
            work_u(i,l) = 0.25*(flux_u(i,l-1) + 2.*flux_u(i,l) + flux_u(i,l+1))
            work_v(i,l) = 0.25*(flux_v(i,l-1) + 2.*flux_v(i,l) + flux_v(i,l+1))
         end do
      end do
      do l = 1,levbot
         do i = 1,nig
            flux_u(i,l) = work_u(i,l)
            flux_v(i,l) = work_v(i,l)
         end do
      end do
   end do

   !  calculate drag at full levels using centered differences
   do l = 2,levbot-1
      do i = 1,nig
         dendz2 = 1. / (density(i,l) * ( alt(i,l-1) - alt(i,l) ))
         drag_u(i,l) = - ( flux_u(i,l-1) - flux_u(i,l) ) * dendz2
         drag_v(i,l) = - ( flux_v(i,l-1) - flux_v(i,l) ) * dendz2
      end do
   end do

   !  drag at first and last levels using one-side differences.
   do i = 1,nig
      dendz = 1. / (density(i,1) * ( alt(i,1) - alt(i,2) ))
      drag_u(i,1) =  flux_u(i,1) * dendz
      drag_v(i,1) =  flux_v(i,1) * dendz
   end do

   do i = 1,nig
      dendz = 1. / (density(i,levbot) * ( alt(i,levbot-1) - alt(i,levbot) ))
      drag_u(i,levbot) = - ( flux_u(i,levbot-1) - flux_u(i,levbot) ) * dendz
      drag_v(i,levbot) = - ( flux_v(i,levbot-1) - flux_v(i,levbot) ) * dendz
   end do
   if (nkm1 > levbot) then
      do i = 1,nig
         dendz = 1. / (density(i,levbot+1) * (alt(i,levbot) - alt(i,levbot+1)))
         drag_u(i,levbot+1) = - flux_u(i,levbot) * dendz
         drag_v(i,levbot+1) = - flux_v(i,levbot) * dendz
      end do
   endif

   return
   !-----------------------------------------------------------------------
end subroutine hines_flux4

end module hines_flux
