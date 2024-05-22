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
!-------------------------------------- LICENCE END --------------------------

module hines_wavnum
   implicit none
   private
   public :: hines_wavnum1

contains

!/@*
subroutine hines_wavnum1(m_alpha, ak_alpha,  &
     &                   v_alpha, visc_mol, density, densb,  &
     &                   bvfreq, rms_wind,  &
     &                   levbot, ni, nig, nkm1, naz)
   use, intrinsic :: iso_fortran_env, only: REAL32
   use mo_gwspectrum, only: kstar, m_min, f1, f2, f3, rnaz
   use hines_sigma, only: hines_sigma1
   use hines_intgrl, only: hines_intgrl4
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: levbot, ni, nig, nkm1, naz
   real(REAL32), intent(inout) :: m_alpha(ni,naz,nkm1) ! cutoff vertical wavenumber (1/m)
   real(REAL32), intent(out) :: ak_alpha(ni,naz)
   real(REAL32), intent(in) :: v_alpha(ni,naz,nkm1)
   real(REAL32), intent(in) :: visc_mol(ni,nkm1)
   real(REAL32), intent(in) :: density(ni,nkm1),  densb(ni)
   real(REAL32), intent(in) :: bvfreq(ni,nkm1),  rms_wind(ni)

   !@Author
   !  aug. 10/95 - c. mclandress
   !  2000-2001  - m. charron
   !  2002       - e. manzini
   !@Object

   !  This routine calculates the cutoff vertical wavenumber and velocity
   !  variances on a longitude by altitude grid for the hines' doppler
   !  spread gravity wave drag parameterization scheme.
   !  note: (1) only values of four or eight can be used for # azimuths (naz).
   !        (2) only values of 1.0, 1.5 or 2.0 can be used for slope (slope).
   !        (3) if m_min not zero, only slope=1. can be used.

   !@Arguments

   !              - Output
   ! m_alpha      cutoff wavenumber at each azimuth (1/m).
   ! sigma_t      total rms horizontal wind (m/s).
   ! sigma_alpha  total rms wind in each azimuth (m/s).
   ! ak_alpha     spectral amplitude factor at each azimuth
   !              (i.e.,{ajkj}) in m^4/s^2.
   ! losigma_t    .true. for total sigma not zero
   ! mmin_alpha   minimum value of cutoff wavenumber.

   !              - Input -
   ! v_alpha      wind component at each azimuth (m/s).
   ! visc_mol     molecular viscosity (m^2/s)
   ! density      background density (kg/m^3).
   ! densb        background density at model bottom (kg/m^3).
   ! bvfreq       background brunt vassala frequency (radians/sec).
   ! bvfb         background brunt vassala frequency at model bottom.
   ! rms_wind     root mean square gravity wave wind at lowest level (m/s).
   ! levbot       index of lowest vertical level.
   ! ni           number of longitudes.
   ! nig          horizontal operator scope
   ! nkm1         number of vertical levels.
   ! naz          azimuthal array dimension

   !             - work arrays -
   ! i_alpha      hines' integral at a single level.
   ! do_alpha     .true. for the azimuths and longitudes for
   !                   which to continue to compute the drag above
   !                   the lowest level
   !*@/

   integer, parameter :: KLEV = 1
   integer, parameter :: KBELOW = KLEV+1
   
   real(REAL32), parameter :: t1 = 1.d-10
   real(REAL32), parameter :: visc_min = 1.e-10
   real(REAL32), parameter :: mmsq = m_min**2
   real(REAL32), parameter :: m_min8 = m_min

   integer :: i, l, n, lbelow

   real(REAL32) :: m_sub_m_turb, m_sub_m_mol, m_trial
   real(REAL32) :: visc, deno

   real(REAL32) :: n_over_m(ni), sigfac(ni), bvfb(ni), rbvfb(ni)
   real(REAL32) :: sigma_t(ni,nkm1)        ! total rms gw wind (m/s)
   real(REAL32) :: sigma_alpha(ni,naz,nkm1)
   real(REAL32) :: mmin_alpha(ni,naz)   ! minumum value of m_alpha

   real(REAL32) :: sigsqh_alpha(ni,naz)
   real(REAL32) :: i_alpha(ni,naz)

   logical :: do_alpha(ni,naz)
   logical :: losigma_t(ni,2)

   !-----------------------------------------------------------------------

   ! initialize logical flags and arrays
   ! calculate azimuthal variances at bottom level using anisotropy factor
   do i = 1,nig
      bvfb(i)  = bvfreq(i,levbot)
      rbvfb(i) = 1. / bvfb(i)  !#TODO: use vect pow fn
   enddo

   do n=1,naz
      do i = 1,nig
         do_alpha(i,n) = .true.
         sigsqh_alpha(i,n) = rnaz * rms_wind(i)**2
      enddo
   enddo

   !  velocity variances at bottom level.
   call hines_sigma1(sigma_t, sigma_alpha, sigsqh_alpha, levbot, ni, nig, nkm1, naz)

   !  calculate cutoff wavenumber and spectral amplitude factor
   !  at bottom level where it is assumed that background winds vanish
   !  and also initialize minimum value of cutoff wavnumber.

   do n = 1,naz
      do i = 1,nig
            m_alpha(i,n,levbot) = bvfb(i) /  &
                 &                          ( f1 * sigma_alpha(i,n,levbot)  &
                 &                          + f2 * sigma_t(i,levbot) )
            ak_alpha(i,n) = 2. * sigsqh_alpha(i,n) &
                 &                        / (m_alpha(i,n,levbot)**2 - mmsq)
            mmin_alpha(i,n) = m_alpha(i,n,levbot)
      end do
   end do

   !  calculate quantities from the bottom upwards,
   !  starting one level above bottom.

   losigma_t(1:nig,KBELOW) = .true.
   DOLEVELS: do l = levbot-1,1,-1
      lbelow = l + 1

      !  calculate n/m_m where m_m is maximum permissible value of the vertical
      !  wavenumber (i.e., m > m_m are obliterated) and n is buoyancy frequency.
      !  m_m is taken as the smaller of the instability-induced
      !  wavenumber (m_sub_m_turb) and that imposed by molecular viscosity
      !  (m_sub_m_mol). since variance at this level is not yet known
      !  use value at level below.

      do i = 1,nig
         losigma_t(i,KLEV) = .true.
         if (losigma_t(i,KBELOW))   then
            visc = max(visc_mol(i,l), visc_min)
            m_sub_m_turb = bvfreq(i,l) / (f2 * sigma_t(i,lbelow))
            m_sub_m_mol = (bvfreq(i,l)*kstar/visc)**0.33333333/f3
            if (m_sub_m_turb < m_sub_m_mol)  then
               n_over_m(i) = f2 *sigma_t(i,lbelow)
            else
               n_over_m(i) = bvfreq(i,l) / m_sub_m_mol
            end if
         endif
      end do

      !  calculate cutoff wavenumber at this level.

      DONAZ: do n = 1,naz
         do i = 1,nig
            if (do_alpha(i,n) .and. losigma_t(i,KBELOW)) then

               ! calculate trial value (variance at this level is not yet known:
               ! use value at level below). if trial value negative or larger
               ! minimum value (not permitted) then set it to minimum value.

               deno= f1 * sigma_alpha(i,n,lbelow) + n_over_m(i) + v_alpha(i,n,l)
               deno= sign(max(t1, abs(deno)), deno)
               m_trial = bvfb(i) / deno

               if (m_trial <= 0.) m_trial = mmin_alpha(i,n)
               m_alpha(i,n,l) = min(m_trial, mmin_alpha(i,n))

               ! do not permit cutoff wavenumber to be less than minimum value.
               m_alpha(i,n,l) = max(m_min8, m_alpha(i,n,l))
               if (m_alpha(i,n,l) == m_min8) do_alpha(i,n) = .false.

               ! reset minimum value of cutoff wavenumber if necessary.
               mmin_alpha(i,n) = min(m_alpha(i,n,l), mmin_alpha(i,n))
               
            else

               m_alpha(i,n,l) = m_min

            endif
         end do
      end do DONAZ

      !  calculate the hines integral at this level.
      call hines_intgrl4(i_alpha, v_alpha, m_alpha, rbvfb,  &
           &             l, ni, nig, nkm1, naz)

      !  calculate the velocity variances at this level.
      do i = 1,nig
         sigfac(i) = densb(i) / density(i,l) * bvfreq(i,l) * rbvfb(i)
      end do

      do n = 1,naz
         do i = 1,nig
            sigsqh_alpha(i,n) = sigfac(i) * ak_alpha(i,n) * i_alpha(i,n)
         end do
      end do
      call hines_sigma1(sigma_t, sigma_alpha, sigsqh_alpha, l, ni, nig, nkm1, naz)

      !  if total rms wind zero (no more drag) then set drag to false
      do i=1,nig
         if (sigma_t(i,l) < epsilon(1.)) losigma_t(i,KLEV) = .false.
         losigma_t(i,KBELOW) = losigma_t(i,KLEV)
      enddo

   end do DOLEVELS

   !-----------------------------------------------------------------------
   return
end subroutine hines_wavnum1

end module hines_wavnum
