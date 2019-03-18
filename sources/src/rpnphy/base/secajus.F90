!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!**s/r secajus - performs a dry convective adjustment
!
      subroutine secajus(tconv, t    , s    , ps,   niter, &
                         conv , cdt1 , ni   , nk)
!
      implicit none
#include <arch_specific.hf>
!
      integer ni,nk, niter
      real cdt1, conv
!
      real ps(ni)
      real s(ni,nk), t(ni,nk), tconv(ni,nk)
!
!Author
!        Alain Patoine
!
!Revision
!
! 001    B. Bilodeau (Jan 1997) - Adaptation from 3D to 2D;
!        add tendencies calculations and dynamic allocation.
! 002    C. Girard (Mar 1997) - Conserve theta rather than T.
! 003    B. Bilodeau (Jan 2001) - Automatic arrays
!
!Object
!     to perform a dry convective adjustment
!
!Notes
!     The algorith is the same as the one used in an old version of the rfe
!     model. Examples of the original code were provided by both C. Beaudoin
!     and B. Bilodeau
!
!Arguments
!
!          - Output -
! tconv    temperature tendency due to dry convective adjustement
! niter    number of iterations
!
!          - Input -
! t        temperature field
! s        sigma levels
! conv     convergence criteria
! cdt1     timestep
! ni       field dimension in x-direction
! nk       field dimension in z-direction
!
!
!Implicits
#include "tdpack_const.hf"
!
!Modules
!     none
!
!*
      integer nitmax
      parameter ( nitmax=25 )
!
      logical adj
!
      integer i, k, nkm, nkp
!
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      real, dimension(ni)    :: wrk11,wrk12,wrk13,wrk14
      real, dimension(ni,nk) :: wrk21,wrk22,wrk23,p
!
!***********************************************************************
!
!
!     calcul de la pression
      do 100 k=1,nk
         do 100 i=1,ni
            p(i,k) = s(i,k) * ps(i)
 100  continue
!
      do 110 i=1,ni
 110  wrk23(i,1)   = (p(i,2)   - p(i,1)  )   * 0.5
!
      do 120 k=2,nk-1
      do 120 i=1,ni
 120  wrk23(i,k)   = (p(i,k+1) - p(i,k-1))   * 0.5
!
      do 130 i=1,ni
 130  wrk23(i,nk) = (p(i,nk) - p(i,nk-1)) * 0.5
!
      do 140 k=1,nk-1
      do 140 i=1,ni
      wrk21(i,k) = (p(i,k)/p(i,k+1))**cappa
      wrk23(i,k) = wrk23(i,k+1) / wrk23(i,k)
      wrk22(i,k) = 1. / (1.+wrk23(i,k))
 140  continue
!
!------------------------------------------------------------------------------
!
      do 260 niter=1,nitmax
!
      adj = .false.
!
      do 210 i=1,ni
 210  wrk11(i) = t(i,1)
!
      DO 240 k=1,nk-1
      nkm = k-1
      nkp = k+1
!
      if ( k .gt. 1 ) then
         do 220 i=1,ni
         tconv(i,nkm) = (wrk11(i) - t(i,nkm))/cdt1
 220     wrk11(i)     = wrk12(i)
      endif
!
      do 230 i=1,ni
      wrk12(i) = t(i,nkp)
      wrk14(i) = wrk11(i)-wrk12(i)*wrk21(i,k)
!
      wrk14(i) = max(0.,-wrk14(i))*wrk22(i,k)
!
      wrk11(i) = wrk11(i)+wrk14(i)*wrk23(i,k)
      wrk12(i) = wrk12(i)-wrk14(i)/wrk21(i,k)
!
      if ( abs(wrk14(i)) .gt. conv ) adj=.true.
!
 230  continue
!
 240  continue
!
      do 250 i=1,ni
      tconv(i,nk-1) = (wrk11(i) - t(i,nk-1))/cdt1
 250  tconv(i,nk  ) = (wrk12(i) - t(i,nk  ))/cdt1
!
!
 260  continue
!
!
      return
      end
