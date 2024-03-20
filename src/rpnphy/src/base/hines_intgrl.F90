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

module hines_intgrl
   implicit none
   private
   public :: hines_intgrl4

contains

!/@*
subroutine hines_intgrl4(i_alpha, v_alpha, m_alpha, rbvfb,   &
     &                  lev, ni, nig, nkm1, naz)
   use, intrinsic :: iso_fortran_env, only: REAL32
   use mo_gwspectrum, only: m_min  !#, slope
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: lev, ni, nig, nkm1, naz
   real(REAL32), intent(out) :: i_alpha(ni,naz)
   real(REAL32), intent(in) :: v_alpha(ni,naz,nkm1)
   real(REAL32), intent(in) :: m_alpha(ni,naz,nkm1)
   real(REAL32), intent(in) :: rbvfb(ni)

   !@Authors
   !  aug. 8/95 - c. mclandress
   !  2001      - m. charron
   !  2003      - l. kornblueh
   !@Object
   !  This routine calculates the vertical wavenumber integral
   !  for a single vertical level at each azimuth on a longitude grid
   !  for the hines' doppler spread gwd parameterization scheme.
   !  note: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
   !        (2) the integral is written in terms of the product qm
   !            which by construction is always less than 1. series
   !            solutions are used for small |qm| and analytical solutions
   !            for remaining values.
   !@Arguments
   !           - Output -
   ! i_alpha   hines' integral.
   !           - Input -
   ! v_alpha   azimuthal wind component (m/s).
   ! m_alpha   azimuthal cutoff vertical wavenumber (1/m).
   ! rbvfb     1./background brunt vassala frequency at model bottom.
   ! m_min     minimum allowable cutoff vertical wavenumber (1/m)
   !           for spectral slope of one.
   ! slope     slope of initial vertical wavenumber spectrum
   !           (must use slope = 1., 1.5 or 2.)
   ! lev       altitude level to process.
   ! ni        number of longitudes.
   ! nig       horizontal operator scope
   ! nkm1      number of vertical levels.
   ! naz       azimuthal array dimension (naz >= naz).

   !  constants in data statements:

   ! qmin = minimum value of q_alpha (avoids indeterminant form of integral)
   ! qm_min = minimum value of q_alpha * m_alpha (used to avoid numerical
   !          problems).
   !*@/

   real(REAL32), parameter :: zero = 0.
   real(REAL32), parameter :: q_min = 1.0
   real(REAL32), parameter :: qm_min = 0.01

   !  internal variables.

   integer :: i, n
   real(REAL32) :: q_alpha, qm, qmm

   !  variables for sparse vector optimization
   integer :: ic, ixi(nig*naz), ix, ixnaz(nig*naz)

   !-----------------------------------------------------------------------

   ic = 0
   do n = 1,naz
      do i = 1,nig

            if (m_alpha(i,n,lev) > m_min) then
               
               q_alpha = v_alpha(i,n,lev) * rbvfb(i)
               qm      = q_alpha * m_alpha(i,n,lev)
               qmm     = q_alpha * m_min

               !  if |qm| is small then use first 4 terms series of taylor
               !  series expansion of integral in order to avoid
               !  indeterminate form of integral,
               !  otherwise use analytical form of integral.

               if (abs(q_alpha) < q_min .or. abs(qm) < qm_min)  then
                  ! taylor series expansion is a very rare event.
                  ! do sparse processing separately
                  ic = ic+1
                  ixi(ic) = i
                  ixnaz(ic) = n
               else
                  i_alpha(i,n) = - (log(1.-qm) - log(1.-qmm) + qm - qmm) / q_alpha**2
                  !  If i_alpha negative due to round off error, set it to zero
                  i_alpha(i,n) = max(i_alpha(i,n), zero)
               end if

            else
               i_alpha(i,n) = 0.
            endif

      end do
   end do
   ! taylor series expansion is a very rare event.
   ! do sparse processing here separately
   do ix = 1, ic
      n = ixnaz(ix)
      i = ixi(ix)
      q_alpha = v_alpha(i,n,lev) * rbvfb(i)
      qm      = q_alpha * m_alpha(i,n,lev)
      qmm     = q_alpha * m_min
      if (abs(q_alpha) < epsilon(1.))  then
         i_alpha(i,n) = (m_alpha(i,n,lev)**2  - m_min**2) / 2.
      else
         i_alpha(i,n) = (qm**2/2.  + qm**3/3.  + qm**4/4.  + qm**5/5.    &
              - qmm**2/2. - qmm**3/3. - qmm**4/4. - qmm**5/5.) &
              / q_alpha**2
      endif
      i_alpha(i,n) = max(i_alpha(i,n), zero)
   end do

   !-----------------------------------------------------------------------
   return
end subroutine hines_intgrl4

end module hines_intgrl
