!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module hines_extro
   implicit none
   private
   public :: hines_extro5

contains

!/@*
subroutine hines_extro5(ni, nig, nkm1, naz,  &
     drag_u, drag_v,  &
     vel_u, vel_v, bvfreq, density, visc_mol, alt,  &
     rmswind,  &
     levbot, hflt)
   use, intrinsic :: iso_fortran_env, only: REAL32
   use mo_gwspectrum, only: m_min
   use hines_wind, only: hines_wind1
   use hines_wavnum, only: hines_wavnum1
   use vert_smooth, only: vert_smooth1
   use hines_flux, only: hines_flux4
   use phy_status, only: phy_error_L
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: ni, nig, nkm1, naz, levbot, hflt

   real, intent(out) :: drag_u(ni,nkm1),   drag_v(ni,nkm1)
   real(REAL32), intent(in) :: vel_u(ni,nkm1),    vel_v(ni,nkm1)
   real(REAL32), intent(in) :: bvfreq(ni,nkm1),   density(ni,nkm1)
   real(REAL32), intent(in) :: visc_mol(ni,nkm1), alt(ni,nkm1)
   real(REAL32), intent(in) :: rmswind(ni)

   !@Authors
   !  aug. 13/95 - c. mclandress
   !  sept. /95  - n. mcfarlane
   !  1995- 2002 - e. manzini
   !@Revision
   !  001   A. Zadra & R. McTaggart-Cowan (June 2011) - provide upper boundary index through interface
   !  001   A. Plante - Remove toplev (notop). Set lev1 to 1 since our model to is well below 120km.
   !@Object
   !  main routine for hines' "extrowave" gravity wave parameterization based
   !  on hines' doppler spread theory. this routine calculates zonal
   !  and meridional components of gravity wave drag, heating rates
   !  and diffusion coefficient on a longitude by altitude grid.
   !  no "mythical" lower boundary region calculation is made.
   !@Argumentsxs
   !              - Output -
   ! drag_u       zonal component of gravity wave drag (m/s^2).
   ! drag_v       meridional component of gravity wave drag (m/s^2).
   ! flux_u       zonal component of vertical momentum flux (pascals)
   ! flux_v       meridional component of vertical momentum flux (pascals)
   !              - Input -
   ! vel_u        background zonal wind component (m/s).
   ! vel_v        background meridional wind component (m/s).
   ! bvfreq       background brunt vassala frequency (radians/sec).
   ! densit       background density (kg/m^3)
   ! visc_mol     molecular viscosity (m^2/s)
   ! alt          altitude of momentum, density, buoyancy levels (m)
   !              (note: levels ordered so that alt(i,1) > alt(i,2), etc.)
   ! rmswind      root mean square gravity wave wind at lowest level (m/s).
   ! k_alpha      horizontal wavenumber of each azimuth (1/m).
   ! levbot       index of last level (eg bottom) for drag calculation
   !              (i.e., lev1 < levbot <= nkm1).
   ! ni           number of longitudes.
   ! nig          horizontal operator scope
   ! nkm1         number of vertical levels.
   ! naz          azimuthal array dimension (naz >= naz).
   !              - ouput diagnostics -
   ! m_alpha      cutoff vertical wavenumber (1/m).
   !              - work arrays -
   ! v_alpha      wind component at each azimuth (m/s) and if iheatcal=1
   !              holds vertical derivative of cutoff wavenumber.
   ! sigma_alpha  total rms wind in each azimuth (m/s).
   ! ak_alpha     spectral amplitude factor at each azimuth
   !              (i.e.,{ajkj}) in m^4/s^2.
   ! densb        background density at bottom level.
   ! bvfb         buoyancy frequency at bottom level and
   !              work array for icutoff = 1.
   ! losigma_t    .true. for total sigma not zero
   !*@/
   
   integer :: i, n, l
   
   real(REAL32) :: m_alpha(ni,naz,nkm1) ! cutoff vertical wavenumber (1/m)
   real(REAL32) :: v_alpha(ni,naz,nkm1)
   real(REAL32) :: densb(ni)
   real(REAL32) :: ak_alpha(ni,naz)
   !-----------------------------------------------------------------------

   !  buoyancy and density at bottom level.
   do i = 1,nig
      densb(i) = density(i,levbot)
   end do

   !  initialize some variables
   do l=1,levbot
      do n = 1,naz
         do i = 1,nig
            m_alpha(i,n,l) =  m_min
         end do
      end do
   end do

   !  compute azimuthal wind components from zonal and meridional winds.
   call hines_wind1(v_alpha, vel_u, vel_v, levbot, ni, nig, nkm1, naz)

   !  calculate cutoff vertical wavenumber and velocity variances.
   call hines_wavnum1(m_alpha, ak_alpha,  &
        &             v_alpha, visc_mol, density, densb,  &
        &             bvfreq, rmswind,  &
        &             levbot, ni, nig, nkm1, naz)
   if (phy_error_L) return

   !  smooth cutoff wavenumbers and total rms velocity in the vertical
   !  direction nsmax times

   call vert_smooth1(m_alpha, levbot, ni, nig, nkm1, naz)
      
   !  calculate zonal and meridional components of the
   !  momentum flux and drag.
   call hines_flux4(drag_u, drag_v,        &
        &           alt, density, densb,   &
        &           m_alpha,  ak_alpha,    &
        &           levbot, ni, nig, nkm1, naz, &
        &           hflt)

   !-----------------------------------------------------------------------
end subroutine hines_extro5

end module hines_extro
