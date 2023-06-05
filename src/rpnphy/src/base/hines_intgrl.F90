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

!/@*
subroutine hines_intgrl(i_alpha,                                     &
     &                  v_alpha, m_alpha, bvfb, m_min, slope, naz,   &
     &                  lev, il1, il2, nlons, nlevs, nazmth,         &
     &                  lorms, do_alpha)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer  lev, naz, il1, il2, nlons, nlevs, nazmth
   real(REAL64) :: i_alpha(nlons,nazmth)
   real(REAL64) :: v_alpha(nlons,nlevs,nazmth)
   real(REAL64) :: m_alpha(nlons,nlevs,nazmth)
   real(REAL64) :: bvfb(nlons), rbvfb(nlons), slope, m_min

   logical lorms(nlons), do_alpha(nlons,nazmth)
   logical lerror(nlons)

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
   ! bvfb      background brunt vassala frequency at model bottom.
   ! m_min     minimum allowable cutoff vertical wavenumber (1/m)
   !           for spectral slope of one.
   ! slope     slope of initial vertical wavenumber spectrum
   !           (must use slope = 1., 1.5 or 2.)
   ! naz       actual number of horizontal azimuths used.
   ! lev       altitude level to process.
   ! il1       first longitudinal index to use (il1 >= 1).
   ! il2       last longitudinal index to use (il1 <= il2 <= nlons).
   ! nlons     number of longitudes.
   ! nlevs     number of vertical levels.
   ! nazmth    azimuthal array dimension (nazmth >= naz).
   ! lorms     .true. for drag computation (column selector)

   !  constants in data statements:

   ! qmin = minimum value of q_alpha (avoids indeterminant form of integral)
   ! qm_min = minimum value of q_alpha * m_alpha (used to avoid numerical
   !          problems).
   !*@/


   !  internal variables.

   integer  i, n
   real(REAL64) :: q_alpha, qm, qmm, sqrtqm, q_min, qm_min
   real(REAL64) :: zero

   !  variables for sparse vector optimization
   integer ic, ixi(nlons*nazmth), ix, ixnaz(nlons*nazmth)

   !-----------------------------------------------------------------------

   !  initialize local scalar and arrays

   zero=0.
   q_min = 1.0
   qm_min = 0.01

   do i = il1,il2
      rbvfb(i)=1.0/bvfb(i)
   enddo

   !  for integer value slope = 1.

   if ( abs(slope-1.) < epsilon(1.) )  then
      ic = 0
      do n = 1,naz
         do i = il1,il2
            if (lorms(i)) then

               if (m_alpha(i,lev,n) > m_min) then
                  q_alpha = v_alpha(i,lev,n) * rbvfb(i)
                  qm      = q_alpha * m_alpha(i,lev,n)
                  qmm     = q_alpha * m_min

                  !  if |qm| is small then use first 4 terms series of taylor
                  !  series expansion of integral in order to avoid
                  !  indeterminate form of integral,
                  !  otherwise use analytical form of integral.

                  if ( abs(q_alpha) .lt. q_min .or. abs(qm).lt. qm_min)  then
                     ! taylor series expansion is a very rare event.
                     ! do sparse processing separately
                     ic = ic+1
                     ixi(ic) = i
                     ixnaz(ic) = n
                  else
                     i_alpha(i,n) = - ( log(1.-qm) - log(1.-qmm) + qm - qmm) / q_alpha**2
                  end if

                  !  If i_alpha negative due to round off error, set it to zero

                  i_alpha(i,n) = max( i_alpha(i,n) , zero )
               else
                  i_alpha(i,n) = 0.
                  do_alpha(i,n) = .false.
               endif

            endif
         end do
      end do
      ! taylor series expansion is a very rare event.
      ! do sparse processing here separately
      do ix =1, ic
         n = ixnaz(ix)
         i = ixi(ix)
         q_alpha = v_alpha(i,lev,n) * rbvfb(i)
         qm      = q_alpha * m_alpha(i,lev,n)
         qmm     = q_alpha * m_min
         if ( abs(q_alpha) < epsilon(1.) )  then
            i_alpha(i,n) = ( m_alpha(i,lev,n)**2  - m_min**2 ) / 2.
         else
            i_alpha(i,n) = ( qm**2/2.  + qm**3/3.  + qm**4/4.  + qm**5/5.    &
                 - qmm**2/2. - qmm**3/3. - qmm**4/4. - qmm**5/5. ) &
                 / q_alpha**2
         end if
         i_alpha(i,n) = max( i_alpha(i,n) , zero )
      end do
   end if

   !  for integer value slope = 2.

   if ( abs(slope-2.) < epsilon(1.) )  then
      ic = 0
      do n = 1,naz
         do i = il1,il2
            if ( lorms(i) ) then

               q_alpha = v_alpha(i,lev,n) * rbvfb(i)
               qm = q_alpha * m_alpha(i,lev,n)

               !  if |qm| is small then use first 4 terms series of taylor
               !  series expansion of integral in order to avoid
               !  indeterminate form of integral,
               !  otherwise use analytical form of integral.

               if ( abs(q_alpha) .lt. q_min .or. abs(qm) .lt. qm_min)  then
                  ! taylor series expansion is a very rare event.
                  ! do sparse processing separately
                  ic = ic+1
                  ixi(ic) = i
                  ixnaz(ic) = n
               else
                  i_alpha(i,n) = - ( log(1.-qm) + qm + qm**2/2.)    &
                       / q_alpha**3
               endif

            endif
         end do
      end do
      ! taylor series expansion is a very rare event.
      ! do sparse processing here separately
      do ix = 1, ic
         n = ixnaz(ix)
         i = ixi(ix)
         q_alpha = v_alpha(i,lev,n) * rbvfb(i)
         qm = q_alpha * m_alpha(i,lev,n)
         if ( abs(q_alpha) < epsilon(1.) )  then
            i_alpha(i,n) = m_alpha(i,lev,n)**3 / 3.
         else
            i_alpha(i,n) = ( qm**3/3. + qm**4/4. + qm**5/5.    &
                 + qm**6/6. ) / q_alpha**3
         end if
      end do
   end if

   !  for real value slope = 1.5

   if ( abs(slope-1.5) < epsilon(1.) )  then
      ic = 0
      do n = 1,naz
         do i = il1,il2
            if ( lorms(i) ) then

               q_alpha = v_alpha(i,lev,n) * rbvfb(i)
               qm = q_alpha * m_alpha(i,lev,n)

               !  if |qm| is small then use first 4 terms series of taylor
               !  series expansion of integral in order to avoid
               !  indeterminate form of integral,
               !  otherwise use analytical form of integral.

               if (abs(q_alpha) .lt. q_min .or. abs(qm) .lt. qm_min)  then
                  ! taylor series expansion is a very rare event.
                  ! do sparse processing separately
                  ic = ic+1
                  ixi(ic) = i
                  ixnaz(ic) = n
               else
                  qm     = abs(qm)
                  sqrtqm = sqrt(qm)
                  if (q_alpha .ge. 0.)  then
                     i_alpha(i,n) = ( log( (1.+sqrtqm)/(1.-sqrtqm) )  &
                          &                          -2.*sqrtqm*(1.+qm/3.) ) / q_alpha**2.5
                  else
                     i_alpha(i,n) = 2. * ( atan(sqrtqm) + sqrtqm*(qm/3.-1.) ) &
                          &                          / abs(q_alpha)**2.5
                  endif
               endif

            endif
         end do
      end do
      ! taylor series expansion is a very rare event.
      ! do sparse processing here separately
      do ix = 1, ic
         n = ixnaz(ix)
         i = ixi(ix)
         q_alpha = v_alpha(i,lev,n) * rbvfb(i)
         qm = q_alpha * m_alpha(i,lev,n)
         if ( abs(q_alpha) < epsilon(1.) )  then
            i_alpha(i,n) = m_alpha(i,lev,n)**2.5 / 2.5
         else
            i_alpha(i,n) = ( qm/2.5 + qm**2/3.5   &
                 + qm**3/4.5 + qm**4/5.5 )   &
                 * m_alpha(i,lev,n)**1.5 / q_alpha
         end if
      enddo
   end if

   !  if integral is negative (which in principal should not happen) then
   !  print a message and some info since execution will abort when calculating
   !  the variances.

   do n = 1,naz
      lerror(:) = .false.
      do i = il1, il2
         if (i_alpha(i,n) < 0.)  then
            lerror(i) = .true.
            exit
         end if
      end do

      if (any(lerror)) then
!!$          WRITE (nout,*)
!!$          WRITE (nout,*) '******************************'
!!$          WRITE (nout,*) 'hines integral i_alpha < 0 '
!!$          WRITE (nout,*) '  longitude i=',i
!!$          WRITE (nout,*) '  azimuth   n=',n
!!$          WRITE (nout,*) '  level   lev=',lev
!!$          WRITE (nout,*) '  i_alpha =',i_alpha(i,n)
!!$          WRITE (nout,*) '  v_alpha =',v_alpha(i,lev,n)
!!$          WRITE (nout,*) '  m_alpha =',m_alpha(i,lev,n)
!!$          WRITE (nout,*) '  q_alpha =',v_alpha(i,lev,n)*rbvfb(i)
!!$          WRITE (nout,*) '  qm      =',v_alpha(i,lev,n)*rbvfb(i)*m_alpha(i,lev,n)
!!$          WRITE (nout,*) '******************************'
         call physeterror('hines_intgrl', 'i_alpha integral is negative')
         return
      end if

   end do

   !-----------------------------------------------------------------------
   return
end subroutine hines_intgrl
